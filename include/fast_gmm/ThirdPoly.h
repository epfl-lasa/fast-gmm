/*
 * ThirdPoly.h
 *
 *  Created on: Oct 30, 2012
 *      Author: seungsu
 */

#pragma once

#include <eigen3/Eigen/Core>

namespace fast_gmm {

class ThirdPoly {
private :
  int mDim;

  Eigen::VectorXd mParamA, mParamB, mParamC, mParamD;

  Eigen::VectorXd mInitPos, mInitVel;
  Eigen::VectorXd mEndPos, mEndVel;
  double mDuration;
public :
  ThirdPoly(int dim);

  void SetConstraints(Eigen::VectorXd& initPos,
                      Eigen::VectorXd& initVel,
                      Eigen::VectorXd& endPos,
                      Eigen::VectorXd& endVel,
                      double duration);

  void Get(double t, Eigen::VectorXd& pos);
  void Get(double t, Eigen::VectorXd& pos, Eigen::VectorXd& vel);
};
}
