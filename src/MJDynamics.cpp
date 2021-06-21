/*
 * MJDynamics.cpp
 *
 *  Created on: May 15, 2013
 *      Author: seungsu
 */

#include "fast_gmm/MJDynamics.h"

#include <iostream>

namespace fast_gmm {

MJDynamics::MJDynamics(int dim, double dt) {
  // set environment
  mDim = dim;
  mDT = dt;

  // resize variables
  mTarget.resize(mDim);
  mState.resize(mDim);
  mStateVelocity.resize(mDim);
  mStateAccel.resize(mDim);
  mStateZerk.resize(mDim);

  mPositionLimitsUp.resize(mDim);
  mPositionLimitsDn.resize(mDim);
  mVelocityLimits.resize(mDim);
  mAccelLimits.resize(mDim);
  mZerkLimits.resize(mDim);

  a0.resize(mDim);
  a1.resize(mDim);
  a2.resize(mDim);
  a3.resize(mDim);
  a4.resize(mDim);
  a5.resize(mDim);

  // set initial values
  mState.setZero();
  mStateVelocity.setZero();
  mStateAccel.setZero();
  mStateZerk.setZero();

  mTarget.setZero();

  mPositionLimitsUp.setZero();
  mPositionLimitsDn.setZero();
  mAccelLimits.setZero();
  mZerkLimits.setZero();

  a0.setZero();
  a1.setZero();
  a2.setZero();
  a3.setZero();
  a4.setZero();
  a5.setZero();

  mReachingTime = 0.0;
}

void MJDynamics::SetState(const Eigen::VectorXd& Position) {
  if (mDim == Position.size()) {
    mState = Position - mTarget;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetState() " << std::endl;
  }
}

void MJDynamics::SetState(const Eigen::VectorXd& Position, const Eigen::VectorXd& Velocity) {
  if ((mDim == Position.size()) && (mDim == Velocity.size())) {
    mState = Position;
    mStateVelocity = Velocity;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetState() " << std::endl;
  }
}

void MJDynamics::SetTarget(const Eigen::VectorXd& target) {
  if (mDim == target.size()) {
    mTarget = target;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetTarget() " << std::endl;
  }
}

void MJDynamics::SetTarget(const Eigen::VectorXd& target, double ReachingTime) {
  if (mDim == target.size()) {
    mTarget = target;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetTarget() " << std::endl;
  }

  mReachingTime = ReachingTime;
  mCurrentTime = 0.0;
}

void MJDynamics::SetStateTarget(const Eigen::VectorXd& Position, const Eigen::VectorXd& Target) {
  SetTarget(Target);
  SetState(Position);
}

void MJDynamics::SetDt(double dt) {
  mDT = dt;
}

void MJDynamics::SetAccelLimits(const Eigen::VectorXd& accelLimits) {
  if (mDim == accelLimits.size()) {
    mAccelLimits = accelLimits;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetAccelLimits() " << std::endl;
  }
}

void MJDynamics::RemoveAccelLimits() {
  mAccelLimits.setZero();
}

double MJDynamics::GetAccelLimits(unsigned int index) {
  if (index < mDim) {
    return mAccelLimits(index);
  } else {
    return 0.0;
  }
}

void MJDynamics::SetZerkLimits(const Eigen::VectorXd& zerkLimits) {
  if (mDim == zerkLimits.size()) {
    mZerkLimits = zerkLimits;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetZerkLimits() " << std::endl;
  }
}

void MJDynamics::RemoveZerkLimits() {
  mZerkLimits.setZero();
}

double MJDynamics::GetZerkLimits(unsigned int index) {
  if (index < mDim) {
    return mZerkLimits(index);
  } else {
    return 0.0;
  }
}

void MJDynamics::SetVelocityLimits(const Eigen::VectorXd& velLimits) {
  if (mDim == velLimits.size()) {
    mVelocityLimits = velLimits;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetVelocityLimits() " << std::endl;
  }
}

void MJDynamics::RemoveVelocityLimits() {
  mVelocityLimits.setZero();
}

double MJDynamics::GetVelocityLimits(unsigned int index) {
  if (index < mDim) {
    return mVelocityLimits(index);
  } else {
    return 0.0;
  }
}

void MJDynamics::SetPositionLimits(const Eigen::VectorXd& posLimitsUp, const Eigen::VectorXd& posLimitsDn) {
  if ((mDim == posLimitsUp.size()) && (mDim == posLimitsDn.size())) {
    mPositionLimitsUp = posLimitsUp;
    mPositionLimitsDn = posLimitsDn;
  } else {
    std::cout << "Dimension error! @ MJDynamics::SetPositionLimits() " << std::endl;
  }
}

void MJDynamics::RemovePositionLimits() {
  mPositionLimitsUp.setZero();
  mPositionLimitsDn.setZero();
}

void MJDynamics::GetTarget(Eigen::VectorXd& target) {
  target = mTarget;
}

void MJDynamics::GetState(Eigen::VectorXd& Position) {
  Position = mState;
}

void MJDynamics::GetState(Eigen::VectorXd& Position, Eigen::VectorXd& Velocity) {
  Position = mState;
  Velocity = mStateVelocity;
}

void MJDynamics::GetStateAccel(Eigen::VectorXd& Accel) {
  Accel = mStateAccel;
}

void MJDynamics::Update() {
  Update(mDT);

}

void MJDynamics::Update(double dt) {
  double lR2, lR3;
  double lD2, lD3, lD4, lD5;

  if (mReachingTime < dt) {
    if (((mTarget - mState).norm() > 0.0001) || (mStateVelocity.norm() > 0.001)) {
      mReachingTime = dt * 2.0;
    }
  }

  if (mReachingTime >= dt) {
    lR2 = mReachingTime * mReachingTime;
    lR3 = lR2 * mReachingTime;
    lD2 = dt * dt;
    lD3 = lD2 * dt;
    lD4 = lD3 * dt;
    lD5 = lD4 * dt;

    a0 = mState;
    a1 = mStateVelocity;
    a2 = mStateAccel / 2.0;

    a3 = (mTarget * 10.0 - a0 * 10.0 - a1 * 6.0 * mReachingTime - a2 * 3.0 * lR2) / (lR3);
    a4 = (a1 * (-2.0) - a2 * 3.0 * mReachingTime - a3 * 3.0 * lR2) / (2.0 * lR3);
    a5 = (a2 * (-1.0) - a3 * 3.0 * mReachingTime - a4 * 6.0 * lR2) / (10.0 * lR3);

    mState = a0 + a1 * dt + a2 * lD2 + a3 * lD3 + a4 * lD4 + a5 * lD5;
    mStateVelocity = a1 + a2 * 2.0 * dt + a3 * 3.0 * lD2 + a4 * 4.0 * lD3 + a5 * 5.0 * lD4;
    mStateAccel = a2 * 2.0 + a3 * 6.0 * dt + a4 * 12.0 * lD2 + a5 * 20.0 * lD3;

    mReachingTime -= dt;
  }
}

double MJDynamics::GetReachingTime() {
  return mReachingTime;
}
}// namespace fast_gmm
