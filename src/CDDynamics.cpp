/*
 * CDDynamics.cpp
 *
 *  Created on: May 29, 2012
 *      Author: Seungsu KIM
 */

#include <iostream>

/*
 * how to use

	CDDynamics *testDyn;
	testDyn = new CDDynamics(dim, dt, wn);

	testDyn->SetVelocityLimits(velLimits);
	testDyn->SetState(initial);
	testDyn->SetTarget(target);


	start loop
		// if new target is set
		testDyn->SetTarget(target);

		// update dynamics
		testDyn->Update();

		// get state
		testDyn->GetState(state);
	end loop
 */

#include "fast_gmm/CDDynamics.h"

namespace fast_gmm {

CDDynamics::CDDynamics(int dim, double dt, double Wn) {
  // set environment
  mDim = dim;
  mDT = dt;
  mWn = Wn;

  // resize variables
  mTarget.resize(mDim);
  mTargetVelocity.resize(mDim);
  mState.resize(mDim);
  mStateVelocity.resize(mDim);
  mStateAccel.resize(mDim);

  mPositionLimits.resize(mDim);
  mVelocityLimits.resize(mDim);
  mAccelLimits.resize(mDim);

  // set initial values
  mState.setZero();
  mStateVelocity.setZero();
  mStateAccel.setZero();

  mTarget.setZero();
  mTargetVelocity.setZero();

  mPositionLimits.setZero();
  mVelocityLimits.setZero();
  mAccelLimits.setZero();

  mReachingTime = 0.0;
}

void CDDynamics::SetState(const Eigen::VectorXd& Position) {
  if (mDim == Position.size()) {
    mState = Position - mTarget;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetState() " << std::endl;
  }
}

void CDDynamics::SetState(const Eigen::VectorXd& Position, const Eigen::VectorXd& Velocity) {
  if ((mDim == Position.size()) && (mDim == Velocity.size())) {
    mState = Position - mTarget;
    mStateVelocity = Velocity - mTargetVelocity;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetState() " << std::endl;
  }
}

void CDDynamics::SetTarget(const Eigen::VectorXd& target) {
  if (mDim == target.size()) {
    mState += (mTarget - target);
    mTarget = target;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetTarget() " << std::endl;
  }
}

void CDDynamics::SetTarget(const double target[]) {
  for (int i = 0; i < mDim; i++) {
    mState[i] += (mTarget[i] - target[i]);
    mTarget[i] = target[i];
  }
}

void CDDynamics::SetTarget(const Eigen::VectorXd& target, double ReachingTime) {
  SetTarget(target);
  mReachingTime = ReachingTime;
}

void CDDynamics::SetStateTarget(const Eigen::VectorXd& Position, const Eigen::VectorXd& Target) {
  SetTarget(Target);
  SetState(Position);
}
/*
void CDDynamics::SetTarget(const Eigen::VectorXd& target, const Eigen::VectorXd& targetVel) {
  if ((mDim == target.size()) && (mDim == targetVel.size())) {
    mState += (mTarget - target);
    mStateVelocity += (mTargetVelocity - targetVel);

    mTarget = target;
    mTargetVelocity = targetVel;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetTarget() " << std::endl;
  }
}
*/
void CDDynamics::SetDt(double dt) {
  mDT = dt;
}

void CDDynamics::SetWn(double Wn) {
  mWn = Wn;
}
void CDDynamics::SetAccelLimits(const Eigen::VectorXd& accelLimits) {
  if (mDim == accelLimits.size()) {
    mAccelLimits = accelLimits;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetAccelLimits() " << std::endl;
  }
}

void CDDynamics::RemoveAccelLimits() {
  mAccelLimits.setZero();
}

double CDDynamics::GetAccelLimits(int index) {
  if (index < mDim) {
    return mAccelLimits(index);
  } else {
    return 0.0;
  }
}

void CDDynamics::SetVelocityLimits(const Eigen::VectorXd& velLimits) {
  if (mDim == velLimits.size()) {
    mVelocityLimits = velLimits;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetVelocityLimits() " << std::endl;
  }
}

void CDDynamics::RemoveVelocityLimits() {
  mVelocityLimits.setZero();
}

double CDDynamics::GetVelocityLimits(int index) {
  if (index < mDim) {
    return mVelocityLimits(index);
  } else {
    return 0.0;
  }
}

void CDDynamics::SetPositionLimits(const Eigen::VectorXd& posLimits) {
  if (mDim == posLimits.size()) {
    mPositionLimits = posLimits;
  } else {
    std::cout << "Dimension error! @ CDDynamics::SetPositionLimits() " << std::endl;
  }
}

void CDDynamics::RemovePositionLimits() {
  mPositionLimits.setZero();
}

void CDDynamics::GetTarget(Eigen::VectorXd& target) {
  target = mTarget;
}
/*
void CDDynamics::GetTarget(Eigen::VectorXd& target, Eigen::VectorXd& targetVel)
{
	target = mTarget;
	targetVel = mTargetVelocity;
}
*/
void CDDynamics::GetState(Eigen::VectorXd& Position) {
  Position = mState + mTarget;
}

void CDDynamics::GetState(double* Position) {
  for (int i = 0; i < mDim; i++) {
    Position[i] = mState[i] + mTarget[i];
  }
}

void CDDynamics::GetState(Eigen::VectorXd& Position, Eigen::VectorXd& Velocity) {
  Position = mState + mTarget;
  Velocity = mStateVelocity + mTargetVelocity;
}

void CDDynamics::GetStateAccel(Eigen::VectorXd& Accel) {
  Accel = mStateAccel;
}

void CDDynamics::Update() {
  Update(mDT, 1.0);
}

void CDDynamics::Update(double dt) {
  Update(dt, 1.0);
}

void CDDynamics::Update(double dt, double muxVel) {
  // A      = x(0);
  // B      = x_d(0) + w*x(0)
  // x(t)   = (A+Bt)e^(-w*t)
  // x_d(t) = (-w*A+(1-w*t)B ) e^(-w*t)

  double lWn = mWn * muxVel;

  mStateAccel = mState * (-lWn * lWn) + mStateVelocity * (-2.0 * lWn);

  for (int i = 0; i < mDim; i++) {
    if (mAccelLimits(i) > 0.0) {
      if (mStateAccel(i) > mAccelLimits(i)) { mStateAccel(i) = mAccelLimits(i); }
      else if (mStateAccel(i) < -mAccelLimits(i)) { mStateAccel(i) = -mAccelLimits(i); }
    }
  }

  mState += mStateVelocity * dt + mStateAccel * dt * dt * 0.5;
  mStateVelocity += mStateAccel * dt;

  for (int i = 0; i < mDim; i++) {
    if (mVelocityLimits(i) > 0) {
      if (mStateVelocity(i) > mVelocityLimits(i)) { mStateVelocity(i) = mVelocityLimits(i); }
      else if (mStateVelocity(i) < -mVelocityLimits(i)) { mStateVelocity(i) = -mVelocityLimits(i); }
    }
  }

  mReachingTime -= dt;
}

double CDDynamics::GetReachingTime(double dt, double muxVel) {
  Eigen::VectorXd lState(mDim);
  Eigen::VectorXd lStateVelocity(mDim);
  Eigen::VectorXd lStateAccel(mDim);
  Eigen::VectorXd lX(mDim);
  Eigen::VectorXd lB(mDim);
  unsigned int frame;
  double lWn;

  lState = mState;
  lStateVelocity = mStateVelocity;

  lWn = mWn * muxVel;
  for (frame = 0; frame < 500; frame++) {
    lStateAccel = lState * (-lWn * lWn) + lStateVelocity * (-2.0 * lWn);
    for (int i = 0; i < mDim; i++) {
      if (mAccelLimits(i) > 0.0) {
        if (lStateAccel(i) > mAccelLimits(i)) { lStateAccel(i) = mAccelLimits(i); }
        else if (lStateAccel(i) < -mAccelLimits(i)) { lStateAccel(i) = -mAccelLimits(i); }
      }
    }

    lState += lStateVelocity * dt;
    lStateVelocity += lStateAccel * dt;

    for (int i = 0; i < mDim; i++) {
      if (mVelocityLimits(i) > 0) {
        if (lStateVelocity(i) > mVelocityLimits(i)) { lStateVelocity(i) = mVelocityLimits(i); }
        else if (lStateVelocity(i) < -mVelocityLimits(i)) { lStateVelocity(i) = -mVelocityLimits(i); }
      }
    }

    if (lState.norm() < 0.001) {
      return (double) frame * dt;
    }
  }
  return (double) frame * dt;
}

double CDDynamics::GetTargetTime() const {
  return mReachingTime;
}
}// namespace fast_gmm