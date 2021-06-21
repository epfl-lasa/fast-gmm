#pragma once

#include <eigen3/Eigen/Core>

namespace fast_gmm {

struct SVRModels {
  unsigned int kernelType;
  unsigned int nbDim;
  unsigned int totalSV;
  double gamma;
  double b;
  Eigen::MatrixXd SVs;  // nbDim X totalSV
  Eigen::VectorXd alpha; // 1 X totalSV
  double mux;
};

// SVR
class SVR {
private:
  SVRModels SVRModel;

  Eigen::VectorXd diffx;
public:
  SVR(char* f_svrmodel);
  double regression(const Eigen::VectorXd& xi);
};
}// namespace fast_gmm
