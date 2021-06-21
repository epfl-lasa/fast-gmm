/*
 * ThridPoly.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: seungsu
 */

#include "fast_gmm/ThirdPoly.h"

namespace fast_gmm {

ThirdPoly::ThirdPoly(int dim) {
  mDim = dim;

  mInitPos.resize(mDim);
  mInitVel.resize(mDim);
  mEndPos.resize(mDim);
  mEndVel.resize(mDim);

  mParamA.resize(mDim);
  mParamB.resize(mDim);
  mParamC.resize(mDim);
  mParamD.resize(mDim);
}

void ThirdPoly::SetConstraints(Eigen::VectorXd& initPos,
                               Eigen::VectorXd& initVel,
                               Eigen::VectorXd& endPos,
                               Eigen::VectorXd& endVel,
                               double duration) {
  mInitPos = initPos;
  mInitVel = initVel;
  mEndPos = endPos;
  mEndVel = endVel;

  // calculate polynomial
  mParamA = mInitPos * (2.0) + mEndPos * (-2.0) + mInitVel + mEndVel;
  mParamB = mInitPos * (-3.0) + mEndPos * (3.0) + mInitVel * (-2.0) - mEndVel;
  mParamC = mInitVel;
  mParamD = mInitPos;

  mDuration = duration;
}

void ThirdPoly::Get(double t, Eigen::VectorXd& pos) {
  double lx = t / mDuration;

  pos = mParamA * lx * lx * lx + mParamB * lx * lx + mParamC * lx + mParamD;
}

void ThirdPoly::Get(double t, Eigen::VectorXd& pos, Eigen::VectorXd& vel) {
  double lx = t / mDuration;

  pos = mParamA * lx * lx * lx + mParamB * lx * lx + mParamC * lx + mParamD;
  vel = mParamA * lx * lx * 3.0 + mParamB * lx * 2.0 + mParamC;
}
}
