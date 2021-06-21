/*
 * CDDynamics.h
 * Critical Damped 2nd order Dynamics
 *
 *  Created on: May 29, 2012
 *      Author: Seungsu KIM
 */

#pragma once

#include <eigen3/Eigen/Core>

namespace fast_gmm {

class CDDynamics {
private :
  Eigen::VectorXd mTarget;
  Eigen::VectorXd mTargetVelocity;

  Eigen::VectorXd mState;
  Eigen::VectorXd mStateVelocity;
  Eigen::VectorXd mStateAccel;

  Eigen::VectorXd mPositionLimits;
  Eigen::VectorXd mVelocityLimits;
  Eigen::VectorXd mAccelLimits;

  int mDim;
  double mWn;
  double mDT;

  double mReachingTime;

public :

  CDDynamics(int dim, double dt, double Wn);

  void SetState(const Eigen::VectorXd& Position);
  void SetState(const Eigen::VectorXd& Position, const Eigen::VectorXd& Velocity);
  void SetTarget(const Eigen::VectorXd& target);
  void SetTarget(const double target[]);
  void SetTarget(const Eigen::VectorXd& target, double ReachingTime);
  void SetStateTarget(const Eigen::VectorXd& Position, const Eigen::VectorXd& Target);
//  void SetTarget(const Eigen::VectorXd& target, const Eigen::VectorXd& targetVel);

  void SetDt(double dt);
  void SetWn(double Wn);

  void SetAccelLimits(const Eigen::VectorXd& velLimits);
  void RemoveAccelLimits(void);
  double GetAccelLimits(int index);

  void SetVelocityLimits(const Eigen::VectorXd& velLimits);
  void RemoveVelocityLimits();
  double GetVelocityLimits(int index);

  void SetPositionLimits(const Eigen::VectorXd& posLimits);
  void RemovePositionLimits();

  void GetTarget(Eigen::VectorXd& target);
  //void GetTarget(Eigen::VectorXd& target, Eigen::VectorXd& targetVel);

  void GetState(Eigen::VectorXd& Position);
  void GetState(double* Position);
  void GetState(Eigen::VectorXd& Position, Eigen::VectorXd& Velocity);
  void GetStateAccel(Eigen::VectorXd& Accel);

  void Update();
  void Update(double dt);
  void Update(double dt, double muxVel);

  double GetReachingTime(double dt, double muxVel);
  double GetTargetTime() const;
};
}// namespace fast_gmm