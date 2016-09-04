#include "pointmatcher/PointMatcher.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <time.h>
#include <map>
#include <octomap/octomap.h>
#include <octomap/OcTree.h>
#include <octomap/ColorOcTree.h>
#include "icp.h"
#include "dbscan/dbscan.h"
#include "ransac/ransac.h"
#include <pcl/visualization/common/common.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>
 
// #include "tracking/TrackManager.h"

using namespace std;
using namespace octomap;
using namespace Eigen;

#define MAX_RANGE 10
#define GROUND_HEIGHT -0.7
#define CEIL_HEIGHT 1

// #define MAX_RANGE 30
// #define GROUND_HEIGHT -0.2
// #define CEIL_HEIGHT 100

typedef PointMatcher<float>::TransformationParameters TransformMatrix;
typedef std::map<int, Pointcloud> MAP;
typedef std::pair<int, Pointcloud> PAIR;
vector<point3d> LidarCenter;
vector<point3d> depthpics;

// global
TransformMatrix TransAcc; // accumulated transform matrix
int total, progress; // for display progress
// TrackManager trackManager;

class CSVRow
{
    public:
        string const& operator[](size_t index) const
        {
            return m_data[index];
        }
        size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(istream& str)
        {
            string         line;
            getline(str,line);

            stringstream   lineStream(line);
            string         cell;

            m_data.clear();
            while(getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        vector<string>    m_data;
};

istream& operator>>(istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

/**
 * read pointcloud from a csvfile
 * @param  filename [description]
 * @return          [description]
 */
Pointcloud readPointCloud(const char* filename) {
  Pointcloud Q;
  ifstream in(filename);
  CSVRow row;
  while(in >> row)
  {
    float x, y, z;
    x = atof(row[0].c_str());
    y = atof(row[1].c_str());
    z = atof(row[2].c_str());
    Q.push_back(x, y, z);
  }
  return Q;
}

void extractGround(Pointcloud P, Pointcloud & ground, Pointcloud &nonGround) {
  point3d upper, lower;
  P.calcBBX(lower, upper);
  // cout << "lower: " << lower << endl;
  // cout << "upper: " << upper << endl;
  point3d ground_upper(upper.x(), upper.y(), lower.z() + 1);
  point3d ground_lower(lower.x(), lower.y(), GROUND_HEIGHT);

  // remove ceiling
  upper = point3d(upper.x(), upper.y(), CEIL_HEIGHT);

  ground.push_back(P);
  nonGround.push_back(P);
  ground.crop(lower, ground_upper);
  nonGround.crop(ground_lower, upper);
}

Pointcloud limitXY(Pointcloud P, float max) {
  point3d upper, lower;
  P.calcBBX(lower, upper);
  upper = point3d(max, max, upper.z());
  lower = point3d(-max, -max, lower.z());
  Pointcloud result;
  result.push_back(P);
  result.crop(lower, upper);
  return result;
}

void getClusterFeatures(Pointcloud cluster, point3d &centroid, point3d & boxSize, point3d &ClusterCenter) {
  int size = cluster.size();
  // calc centroid
  float x = 0, y = 0, z = 0;
  for (int i = 0; i < size; i++) {
    x += cluster[i].x();
    y += cluster[i].y();
    z += cluster[i].z();
  }
  centroid = point3d(x/size, y/size, z/size);
  // calc boxSize
  point3d upper, lower;
  cluster.calcBBX(lower, upper);
  boxSize = point3d(upper.x() - lower.x(), upper.y() - lower.y(), upper.z() - lower.z());
  ClusterCenter = point3d((upper.x()+lower.x())/2, (upper.y()+lower.y())/2, (upper.z()+lower.z())/2);
}

void getInAndOutliners(Pointcloud P, bool * cs, Pointcloud &inliners, Pointcloud &outliners) {
  inliners.clear();
  outliners.clear();
  int size = P.size();
  int count = 0;
  for (int i = 0; i < size; i++) {
    if (cs[i]) {
      inliners.push_back(P[i]);
    } else {
      outliners.push_back(P[i]);
      count++;
    }
  }
}

void octree2pcl(Pointcloud tree){
  // unsigned int maxDepth = tree.getTreeDepth();
  // cout << "tree depth is " << maxDepth << endl;
  // expand collapsed occupied nodes until all occupied leaves are at maximum depth
  vector<point3d> pcl;
  for (Pointcloud::iterator it = tree.begin(); it != tree.end(); ++it)
  {
    pcl.push_back(*it);
  }
 string outputFilename = "outliners.pcd";
  ofstream f(outputFilename.c_str(), ofstream::out);
  f << "# .PCD v0.7" << endl
    << "VERSION 0.7" << endl
    << "FIELDS x y z" << endl
    << "SIZE 4 4 4" << endl
    << "TYPE F F F" << endl
    << "COUNT 1 1 1" << endl
    << "WIDTH " << pcl.size() << endl
    << "HEIGHT 1" << endl
    << "VIEWPOINT 0 0 0 0 0 0 1" << endl
    << "POINTS " << pcl.size() << endl
    << "DATA ascii" << endl;
  for (size_t i = 0; i < pcl.size(); i++)
      f << pcl[i].x() << " " << pcl[i].y() << " " << pcl[i].z() << endl;
  f.close();
}

// pcl EduclideanClusterExtraction method
// void cluster_extraction(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){
//   pcl::PCDWriter writer;
//   // Creating the KdTree object for the search method of the extraction
//   pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
//   tree->setInputCloud (cloud);

//   std::vector<pcl::PointIndices> cluster_indices;
//   pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
//   ec.setClusterTolerance (0.07); //7cm
//   ec.setMinClusterSize (5);
//   ec.setMaxClusterSize (1000);
//   ec.setSearchMethod (tree);
//   ec.setInputCloud (cloud);
//   ec.extract (cluster_indices);
//   int j = 0;
//   for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin (); it != cluster_indices.end (); ++it)
//   {
//     pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster (new pcl::PointCloud<pcl::PointXYZ>);
//     for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
//       cloud_cluster->points.push_back (cloud->points[*pit]); //*
//     cloud_cluster->width = cloud_cluster->points.size ();
//     cloud_cluster->height = 1;
//     cloud_cluster->is_dense = true;

//     std::cout << "PointCloud representing the Cluster: " << cloud_cluster->points.size () << " data points." << std::endl;
//     std::stringstream ss;
//     ss << "cloud_cluster_" << j << ".pcd";
//     writer.write<pcl::PointXYZ> (ss.str (), *cloud_cluster, false); //*
//     j++;
//   }
// }

void deleteCrop(Pointcloud &input, point3d lowerBound, point3d upperBound) {

  Pointcloud result;

  float min_x, min_y, min_z;
  float max_x, max_y, max_z;
  float x,y,z;

  min_x = lowerBound(0); min_y = lowerBound(1); min_z = lowerBound(2);
  max_x = upperBound(0); max_y = upperBound(1); max_z = upperBound(2);

  for (Pointcloud::const_iterator it=input.begin(); it!=input.end(); it++) {
    x = (*it).x();
    y = (*it).y();
    z = (*it).z();

    if ((x >= min_x) &&
   (y >= min_y) &&
   (z >= min_z) &&
   (x <= max_x) &&
   (y <= max_y) &&
   (z <= max_z)){
    }
    else{
      result.push_back (x,y,z);
    }
  }

  input.clear();
  input.push_back(result);

}

// DBSCAN Cluster Extraction
void cluster_extraction(Pointcloud dcs, ColorOcTree tree, MAP &movingObjs, TransformMatrix TransAcc, std::vector<MAP> &DynamicObjects, Pointcloud &stationary){

  int size = dcs.size();
  int * clusters_idxs = new int[size];
 
  clusters_idxs = dbscan(dcs, 10, 0.5);//min_points and epsilon

  MAP clusterMap;
  for(int i=0; i<size; i++){
    int cluster_idx = clusters_idxs[i];
    MAP::iterator it = clusterMap.find(cluster_idx);
    float x = dcs[i].x();
    float y = dcs[i].y();
    float z = dcs[i].z();
    if (it != clusterMap.end()) {
      (it->second).push_back(point3d(x, y, z));
    } else {
      Pointcloud v;
      v.push_back(point3d(x, y, z));
      clusterMap.insert(PAIR(cluster_idx, v));
    }
  }

  cout<<"Cluster Size: "<<clusterMap.size()<<endl;
  pcl::PCDWriter writer;
  int exist = 0, movingObjs_index=0;
  float score = 0;
  for (MAP::iterator it = clusterMap.begin() ; it != clusterMap.end(); it++) {
    Pointcloud cluster = it->second;
    Pointcloud tmp(stationary);
    point3d lowerBound, upperBound;
    cluster.calcBBX(lowerBound, upperBound);
    point3d ex(0.3, 0.3, 0.3);
    // point3d ex(0.1, 0.1, 0.1);
    lowerBound -= ex;
    upperBound += ex;
    tmp.crop(lowerBound, upperBound);
    cluster.clear();
    for(Pointcloud::iterator it = tmp.begin(); it != tmp.end(); it++)
    {
      cluster.push_back(*it);
    }
    int cluster_idx = it->first;
    exist = 0;
    // compare to total map to remove some not interest cluster 
    for(Pointcloud::iterator it = cluster.begin(); it != cluster.end(); it++)
    {
      OcTreeNode* n = tree.search((*it));
      if(!n){
        exist++;
      }
    }
    score = exist/cluster.size();
    if( cluster_idx != 0 && score>0.5 && cluster.size()>50)
    {
      movingObjs.insert(PAIR(movingObjs_index,cluster));
      point3d temp;
      temp.x() = TransAcc(0, 3);
      temp.y() = TransAcc(1, 3);
      temp.z() = TransAcc(2, 3);
      LidarCenter.push_back(temp);
      pcl::PointCloud<pcl::PointXYZ> cloud;
      cloud.width = cluster.size() +1;
      cloud.height = 1;
      cloud.is_dense = false;
      cloud.points.resize (cloud.width * cloud.height); 
      // Lidar Center
      cloud.points[0].x = temp.x();
      cloud.points[0].y = temp.y();
      cloud.points[0].z = temp.z();
      int j = 1;
      for (Pointcloud::iterator it = cluster.begin(); it != cluster.end(); ++it)
      {
        // cloud.points[j].x = (*it).x() - temp.x();
        // cloud.points[j].y = -((*it).z() - temp.z());
        // cloud.points[j].z = (*it).y() - temp.y();
        cloud.points[j].x = (*it).x();
        cloud.points[j].y = (*it).y();
        cloud.points[j].z = (*it).z();
        j++;
      }
      std::stringstream ss;
      ss << "cloud_cluster_"  << DynamicObjects.size() << "_" << movingObjs_index << ".pcd";
      string outputFilename;
      ss >> outputFilename;
      pcl::io::savePCDFileBinary(outputFilename, cloud);
      movingObjs_index++; 
      deleteCrop(stationary, lowerBound, upperBound);
    }
  }
  DynamicObjects.push_back(movingObjs);
  // octree2pcl(stationary);
}

void transform2TSDF(Pointcloud dynObj, point3d LidarCenter,vector<point3d> &depthpics){
  point3d center, boxSize, centroid, normal, v, dynPoint, projection, projectionCenter, temp;
  depthpics.clear();
  getClusterFeatures(dynObj, centroid, boxSize, center) ;
  normal.x() = center.x() - LidarCenter.x();
  normal.y() = center.y() - LidarCenter.y();
  normal.z() = center.z() - LidarCenter.z(); 
  float total = sqrt(pow((normal.x()), 2) + pow((normal.y()),2) + pow((normal.z()), 2));
  normal.x() = normal.x()/total;
  normal.y() = normal.y()/total;
  normal.z() = normal.z()/total;
  for (Pointcloud::iterator it = dynObj.begin(); it != dynObj.end(); it++) {
    dynPoint.x() = (*it).x();
    dynPoint.y() = (*it).y();
    dynPoint.z() = (*it).z();
    v.x() = dynPoint.x() - LidarCenter.x();
    v.y() = dynPoint.y() - LidarCenter.y();
    v.z() = dynPoint.z() - LidarCenter.z();
    float depth = normal.dot(v);
    projection.x() = (0.2*LidarCenter.x()+depth*dynPoint.x())/(0.2+depth);
    projection.y() = (0.2*LidarCenter.y()+depth*dynPoint.y())/(0.2+depth);
    // projection.z() = (0.2*LidarCenter.z()+depth*dynPoint.z())/(0.2+depth);
    total = sqrt(pow((center.x()-LidarCenter.x()), 2) + pow((center.y()-LidarCenter.y()),2) + pow((center.z()-LidarCenter.z()), 2));
    projectionCenter.x() = 0.2*(center.x()-LidarCenter.x())/total+LidarCenter.x();
    projectionCenter.y() = 0.2*(center.y()-LidarCenter.y())/total+LidarCenter.y();
    // projectionCenter.z() = 0.2*(center.z()-LidarCenter.z())/total+LidarCenter.z();
    temp.x() = projection.x() - projectionCenter.x();
    temp.y() = projection.y() - projectionCenter.y();
    temp.z() = depth;
    depthpics.push_back(temp);
  }
}

void initMap(ColorOcTree &tree, Pointcloud P) {
  Pointcloud ground, PWithOutGround;
  extractGround(P, ground, PWithOutGround);
  P = PWithOutGround;
  P = limitXY(P, MAX_RANGE);

  for (Pointcloud::iterator it = P.begin(); it != P.end(); it++) {
    tree.updateNode((*it).x(), (*it).y(), (*it).z(), true);
    tree.setNodeColor((*it).x(), (*it).y(), (*it).z(), 0, 0, 255);
  }
}

Pointcloud updateMap(ColorOcTree &tree, Pointcloud P, Pointcloud lastP, std::vector<Pointcloud> &dcs, std::vector<MAP> &DynamicObjects) {
  long beginTime = clock();
  Pointcloud ground, PWithOutGround, P_, inliners1, outliners1, inliners2, outliners2, previousOutliners, dynObj;
  dynObj.clear();
  // icp and dynamic dectction should be implemented without ground points
  extractGround(P, ground, PWithOutGround);
  P = PWithOutGround;
  P = limitXY(P, MAX_RANGE);
  bool * cs = new bool[P.size()];
  bool *cs2= new bool[lastP.size()];
  P = icp(lastP, P, TransAcc, cs, cs2);
  getInAndOutliners(P, cs, inliners1, outliners1);
  // dynamic detection
  MAP movingObjs;
  dcs.push_back(outliners1);
  // if(dcs.size()>=1)
  // {
    previousOutliners = dcs[0];
    dcs.erase(dcs.begin());
    cluster_extraction(previousOutliners, tree, movingObjs, TransAcc, DynamicObjects, P);
  // }
  for(int i=0;i<P.size();i++)
  {
    P_.push_back(P[i]);
  }

  cout<<"P_:"<<P_.size()<<endl;
  // add points
  for (Pointcloud::iterator it = P_.begin(); it != P_.end(); it++) {
    tree.updateNode((*it), true);
    tree.setNodeColor((*it).x(), (*it).y(), (*it).z(), 0, 0, 255);
  }
   //free points update
  // mark and clear dynamic points
  for (int i = 0; i < movingObjs.size(); i++) {
    dynObj = movingObjs[i];
    // transform2TSDF(dynObj, LidarCenter, depthpics);
    for (Pointcloud::iterator it = dynObj.begin(); it != dynObj.end(); it++) {
      ColorOcTreeNode* node = tree.updateNode((*it), false);
      node->setLogOdds(-0.4);
      // ColorOcTreeNode* n = tree.updateNode((*it), true);
      // n->setColor(255,0,0); // set color to red
    }
  }

  tree.setNodeColor(TransAcc(0, 3), TransAcc(1, 3), TransAcc(2, 3), 155, 100, 255); // lidar current pos

  long endTime = clock();
  char msg[100];
  sprintf(msg, "frame %d/%d completed, consumed time: %.2f s.\n", progress, total, (float)(endTime-beginTime)/1000000);
  cout << msg;
  octree2pcl(P);
  delete[] cs;
  delete[] cs2;
  return P;
}



void viewerOneOff (pcl::visualization::PCLVisualizer& viewer)
{
    // set background to black (R = 0, G = 0, B = 0)
    viewer.setBackgroundColor (0, 0, 0);
}

void viewerPsycho (pcl::visualization::PCLVisualizer& viewer)
{
    // you can add something here, ex:  add text in viewer
}

int main(int argc, char** argv) {
  ColorOcTree tree (0.05);  // create empty tree with resolution 0.1
  int from  = atoi(argv[1]);
  int to = atoi(argv[2]);
  int step = atoi(argv[3]);
  string path = argv[4];

  std::vector<Pointcloud> dcs;//Dynamic Objects Candidates
  std::vector<MAP> DynamicObjects;

  // pcl visualization
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);

  // ptr transform
  pcl::visualization::CloudViewer viewer("Cloud Viewer");

  viewer.runOnVisualizationThreadOnce(viewerOneOff);
  viewer.runOnVisualizationThreadOnce(viewerPsycho);

  // init
  char baseFile[50];
  sprintf(baseFile, "%s (Frame %04d).csv", path.c_str(), from);
  Pointcloud base = readPointCloud(baseFile);
  initMap(tree, base);
  TransAcc.resize(4, 4);
  TransAcc.setIdentity();
  total = (int) (to - from) / step;
  progress = 1;

  Pointcloud P, lastP;
  lastP = base;
  char file[50];
  for (int i = from + 1; i <= to; i += step) {
    sprintf(file, "%s (Frame %04d).csv", path.c_str(), i);
    P = readPointCloud(file);
    lastP = updateMap(tree, P, lastP, dcs, DynamicObjects);
    progress++;
    pcl::io::loadPCDFile("outliners.pcd", *cloud);
    viewer.showCloud(cloud);
  }
  // trackManager.saveTargets();

  // string result = "map.bt";
  // tree.writeBinary(result);

  string result = "map.ot";
  tree.write(result);

  cout << "wrote example file " << result << endl;

  while(! viewer.wasStopped())
  {
  }
  return 0;
}
