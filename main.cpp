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

// global
TransformMatrix TransAcc; // accumulated transform matrix
int total, progress; // for display progress
// TrackManager trackManager;

class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str,line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str,CSVRow& data)
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
  std::ifstream in(filename);
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

// void extractGround(Pointcloud P, Pointcloud & ground, Pointcloud &nonGround) {
//   bool *isGround;
//   isGround = ransacFitPlane(P, 0.5, 6000, 500);
//   int size = P.size();
//   for (int i = 0; i < size; i++) {
//     if (isGround[i]) {
//       ground.push_back(P[i]);
//     } else {
//       nonGround.push_back(P[i]);
//     }
//   }
// }

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

void getClusterFeatures(Pointcloud cluster, point3d &centroid, point3d & boxSize) {
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
}

void getInAndOutliners(Pointcloud P, Pointcloud Q, bool * cs, bool *cs2, Pointcloud &inliners, Pointcloud &outliners, Pointcloud &inliners2, Pointcloud &outliners2) {
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
  inliners2.clear();
  outliners2.clear();
  size = Q.size();
  count = 0;
  for (int i = 0; i < size; i++) {
    if (cs2[i]) {
      inliners2.push_back(Q[i]);
    } else {
      outliners2.push_back(Q[i]);
      count++;
    }
  }
}

Pointcloud dynamicFilter(Pointcloud P_origion, Pointcloud P, Pointcloud Q, std::vector<Pointcloud> &movingObjs){
  Pointcloud inliners;
  Pointcloud outliners;
  float throttle = 0.05;
  movingObjs.clear();
  // match
  Matrix<float, 4, Dynamic> features_Q = getFeaturesMatrix(Q);
  Matrix<float, 4, Dynamic> features_P = getFeaturesMatrix(P);

  Labels featureLabels;
  featureLabels.push_back(Label("x"));
  featureLabels.push_back(Label("y"));
  featureLabels.push_back(Label("z"));

  DP data1(features_P, featureLabels);
  DP data2(features_Q, featureLabels);
  MatrixXf M = data1.features;
  MatrixXf N = data2.features;

  NNSearchF* nns = NNSearchF::createKDTreeLinearHeap(M);
  const int K = 1;
  MatrixXi indices;
  MatrixXf dists2;
  indices.resize(K, N.cols());
  dists2.resize(K, N.cols());
  nns->knn(N, indices, dists2, K, 0, throttle);

  for (int i = 0; i < N.cols(); i++) {
    if (dists2(i) < throttle && indices(i)) {
      inliners.push_back(P[i]);
    } 
    else {
      outliners.push_back(P[i]);
    }
  }
  cout<<"dynamicObjects size:"<<inliners.size()<<endl;
  for(int i=0;i<P_origion.size();i++)
  {
    outliners.push_back(P_origion[i]);
  }
  movingObjs.push_back(inliners);
  return outliners;
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

Pointcloud updateMap(ColorOcTree &tree, Pointcloud P, Pointcloud lastP, Pointcloud lastOutliners) {
  long beginTime = clock();
  Pointcloud ground, PWithOutGround, P_, inliners1, outliners1, inliners2, outliners2;
  // icp and dynamic dectction should be implemented without ground points
  extractGround(P, ground, PWithOutGround);
  P = PWithOutGround;
  P = limitXY(P, MAX_RANGE);
  bool * cs = new bool[P.size()];
  bool *cs2= new bool[lastP.size()];
  P = icp(lastP, P, TransAcc, cs, cs2);
  getInAndOutliners(P, lastP, cs, cs2, inliners1, outliners1, inliners2, outliners2);
  cout<<"inliners1: "<<inliners1.size()<<"inliners2: "<<inliners2.size()<<endl;
  cout<<"lastOutliners: "<<lastOutliners.size()<<"outliners2: "<<outliners2.size()<<endl;
  // dynamic detection
  std::vector<Pointcloud> movingObjs;

  P_ = dynamicFilter(inliners1, outliners2, lastOutliners, movingObjs);
  // trackManager.update(movingObjs);
  cout<<"P_:"<<P_.size()<<endl;
  // add points
  for (Pointcloud::iterator it = P_.begin(); it != P_.end(); it++) {
    tree.updateNode((*it), true);
    tree.setNodeColor((*it).x(), (*it).y(), (*it).z(), 0, 0, 255);
  }

  // mark and clear dynamic points
  for (int i = 0; i < movingObjs.size(); i++) {
    Pointcloud dynObj = movingObjs[i];
    for (Pointcloud::iterator it = dynObj.begin(); it != dynObj.end(); it++) {
      ColorOcTreeNode* node = tree.updateNode((*it), false);
      node->setLogOdds(-0.4);
      // ColorOcTreeNode* n = tree.updateNode((*it), true);
      // n->setColor(255,0,0); // set color to red
    }
  }

  tree.setNodeColor(TransAcc(0, 3), TransAcc(1, 3), TransAcc(2, 3), 255, 0, 255); // lidar current pos

  long endTime = clock();
  char msg[100];
  sprintf(msg, "frame %d/%d completed, consumed time: %.2f s.\n", progress, total, (float)(endTime-beginTime)/1000000);
  cout << msg;
  lastOutliners.clear();
  int size = outliners1.size();
  for(int i=0;i<size;i++)
  {
    lastOutliners.push_back(outliners1[i]);
  }
  delete[] cs;
  delete[] cs2;
  return P;
}

int main(int argc, char** argv) {
  ColorOcTree tree (0.05);  // create empty tree with resolution 0.1
  int from  = atoi(argv[1]);
  int to = atoi(argv[2]);
  int step = atoi(argv[3]);
  string path = argv[4];

  // init
  char baseFile[50];
  sprintf(baseFile, "%s (Frame %04d).csv", path.c_str(), from);
  Pointcloud base = readPointCloud(baseFile);
  initMap(tree, base);
  TransAcc.resize(4, 4);
  TransAcc.setIdentity();
  total = (int) (to - from) / step;
  progress = 1;

  Pointcloud P, lastP, lastOutliners;
  lastP = base;
  char file[50];
  for (int i = from + 1; i <= to; i += step) {
    sprintf(file, "%s (Frame %04d).csv", path.c_str(), i);
    P = readPointCloud(file);
    lastP = updateMap(tree, P, lastP, lastOutliners);
    progress++;
  }

  // showMovingObjsTrajectory(tree);
  // trackManager.saveTargets();

  // string result = "map.bt";
  // tree.writeBinary(result);

  string result = "map.ot";
  tree.write(result);

  cout << "wrote example file " << result << endl;

  return 0;
}
