#include "icp.h"
// #include <algorithm>


Matrix<float, 4, Dynamic> getFeaturesMatrix(Pointcloud P) {
  int nPoints = P.size();
  Matrix<float, 4, Dynamic> features;
  features.resize(4, nPoints);

  int col = 0;
  for (Pointcloud::iterator it = P.begin(); it != P.end(); it++) {
    features(0, col) = (*it).x();
    features(1, col) = (*it).y();
    features(2, col) = (*it).z();
    features(3, col) = 1;
    col ++;
  }
  // cout << features << endl;
  return features;
}

Pointcloud toPointCloud(Matrix<float, 4, Dynamic> features) {
  Pointcloud pc;
  int n = features.outerSize();
  for (int i = 0; i < n; i++) {
    float x = features(0, i) / features(3, i);
    float y = features(1, i) / features(3, i);
    float z = features(2, i) / features(3, i);
    pc.push_back(x, y, z);
  }
  return pc;
}

//total transformation
void accumulateTransfor(TransformMatrix & Tacc, TransformMatrix Tnew) {
  TransformMatrix R = Tacc.block(0, 0, 3, 3);
  TransformMatrix T = Tacc.block(0, 3, 3, 1);
  TransformMatrix Rt = Tnew.block(0, 0, 3, 3);
  TransformMatrix Tt = Tnew.block(0, 3, 3, 1);
  // TransformMatrix identityPart = Tacc.block(3, 0, 1, 4);
  R = Rt * R;
  T = Rt * T + Tt;
  Tacc << R, T;
}

/*code for ransacIcp*/
Matrix<float, 4, Dynamic> randomPick(Pointcloud P, int n) {
  random_shuffle(P.begin(), P.end());
  Pointcloud samples;
  for (int i = 1; i < n; i++) {
    samples.push_back(P[i]);
  }
  return getFeaturesMatrix(samples);
}

Pointcloud findInliners(DP Q, DP P_, TransformMatrix T, float throttle, bool *cs, bool *cs2) {
  Pointcloud inliners;

  // apply transformation
  DP P = P_;
  PM::Transformation* rigidTrans;
  rigidTrans = PM::get().REG(Transformation).create("RigidTransformation");
  P = rigidTrans->compute(P, T);

  // match
  MatrixXf M = Q.features;
  MatrixXf q = P.features;
  MatrixXf P_origin = P_.features;
  //find inliners1
  NNSearchF* nns = NNSearchF::createKDTreeLinearHeap(M);
  const int K = 1;
  MatrixXi indices;
  MatrixXf dists2;
  indices.resize(K, q.cols());
  dists2.resize(K, q.cols());
  nns->knn(q, indices, dists2, K, 0, throttle);

  for (int i = 0; i < q.cols(); i++) {
    Vector4f point = P_origin.col(i);
    if (dists2(i) < throttle && indices(i)) {
      inliners.push_back(point(0), point(1), point(2));
      cs[i] = 1;
    } else {
      cs[i] = 0;
    }
  }
  //find inliners2
  NNSearchF* nns2 = NNSearchF::createKDTreeLinearHeap(q);
  MatrixXi indices2;
  MatrixXf dists3;
  indices2.resize(K, M.cols());
  dists3.resize(K, M.cols());
  nns->knn(M, indices2, dists3, K, 0, throttle);

  for (int i = 0; i < M.cols(); i++) {
    if (dists3(i) < throttle && indices2(i)) {
      cs2[i] = 1;
    } else {
      cs2[i] = 0;
    }
  }  

  return inliners;
}

float calcErr(DP Q, DP P, TransformMatrix T) {
  // apply transformation
  PM::Transformation* rigidTrans;
  rigidTrans = PM::get().REG(Transformation).create("RigidTransformation");
  P = rigidTrans->compute(P, T);

  // match
  MatrixXf M = Q.features;
  MatrixXf q = P.features;
  NNSearchF* nns = NNSearchF::createKDTreeLinearHeap(M);
  const int K = 1;
  MatrixXi indices;
  MatrixXf dists2;
  indices.resize(K, q.cols());
  dists2.resize(K, q.cols());
  nns->knn(q, indices, dists2, K, 0.1);
  return dists2.mean();
}

DP matrix2DP(Matrix<float, 4, Dynamic> features) {
  Labels featureLabels;
  featureLabels.push_back(Label("x"));
  featureLabels.push_back(Label("y"));
  featureLabels.push_back(Label("z"));
  return DP(features, featureLabels);
}

TransformMatrix ransacIcp(DP ref, DP data, bool *cs,bool *cs2) {
  long beginTime = clock();
  const int MAX_ITER = 10;
  const float pIn = 0.6;
  float throttle = 0.05;
  float errMin = -1;
  float errAcc = 0.1;
  float errStop = 0.05;

  // setup icp
  PM::ICP icp;
  string configFile = "icp_config.yaml";
  ifstream ifs(configFile.c_str());
  icp.loadFromYaml(ifs);
  // icp.setDefault();

  int dataSize = data.features.cols();
  Pointcloud P = toPointCloud(data.features);
  Matrix<float, 4, Dynamic> sample;
  TransformMatrix bestT;
  int i;
  bool founded = false;
  for (i = 0; i < MAX_ITER; i++) {
    sample = randomPick(P, ceil(dataSize * pIn));

    DP dataTmp = matrix2DP(sample);
    TransformMatrix T = icp(dataTmp, ref);

    Pointcloud inliners = findInliners(ref, data, T, throttle, cs, cs2);
    cout << "inliners:" << (float)inliners.size()/dataSize << endl;

    if ((float)inliners.size()/dataSize > pIn) {

      DP inlinersDp = matrix2DP(getFeaturesMatrix(inliners));

      T = icp(inlinersDp, ref);

      float err = calcErr(ref, inlinersDp, T);
      if (err < errAcc && (errMin < 0 || err < errMin)) {
        founded = true;
        errMin = err;
        bestT = T;
        if (errMin < errStop) break;
      }
    }
  }

  if (!founded) {
    cout << "ransac icp falled back" << endl;
    bestT = icp(data, ref);
    findInliners(ref, data, bestT, throttle, cs, cs2);
  }

  // cout << "ransac icp completed in " << i << "iterations" << endl;
  // cout << "final err: " << errMin << endl;

  long endTime = clock();
  char msg[100];
  sprintf(msg, "ransac icp consume time: %.2f s.\n", (float)(endTime-beginTime)/1000000);
  cout << msg;

  return bestT;
}
/*******************/

/**
 * perform icp algorithm
 * @param  Q    [reference pointcloud]
 * @param  P    [pointcloud to be aligned]
 * @param  Tacc [accumulated transform matrix]
 * @return      [aligned P]
 */
Pointcloud icp(Pointcloud Q, Pointcloud P, TransformMatrix & Tacc, bool * cs,bool *cs2) {
  Matrix<float, 4, Dynamic> features_Q = getFeaturesMatrix(Q);
  Matrix<float, 4, Dynamic> features_P = getFeaturesMatrix(P);

  Labels featureLabels;
  featureLabels.push_back(Label("x"));
  featureLabels.push_back(Label("y"));
  featureLabels.push_back(Label("z"));

  DP data(features_P, featureLabels);
  DP ref(features_Q, featureLabels);

  PM::Transformation* rigidTrans;
  rigidTrans = PM::get().REG(Transformation).create("RigidTransformation");
  data = rigidTrans->compute(data, Tacc); // init transform

  // icp with ransac
  TransformMatrix T = ransacIcp(ref, data, cs, cs2);
  cout<< endl << "ransacT: " << endl << T;

  DP data_out(data);
  data_out = rigidTrans->compute(data_out, T);
  accumulateTransfor(Tacc, T);

  return toPointCloud(data_out.features);
}
