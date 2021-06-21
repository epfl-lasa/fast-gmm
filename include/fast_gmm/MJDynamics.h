/*
 * MJDynamics.h
 *
 *  Created on: May 15, 2013
 *      Author: seungsu
 */

#pragma once

#include <eigen3/Eigen/Core>

namespace fast_gmm {

class MJDynamics {
private :
  Eigen::VectorXd mTarget;

  Eigen::VectorXd mState;
  Eigen::VectorXd mStateVelocity;
  Eigen::VectorXd mStateAccel;
  Eigen::VectorXd mStateZerk;

  Eigen::VectorXd mPositionLimitsUp;
  Eigen::VectorXd mPositionLimitsDn;
  Eigen::VectorXd mVelocityLimits;
  Eigen::VectorXd mAccelLimits;
  Eigen::VectorXd mZerkLimits;

  Eigen::VectorXd a0, a1, a2, a3, a4, a5;

  unsigned int mDim;
  double mDT;

  double mReachingTime;
  double mCurrentTime;
public :

  MJDynamics(int dim, double dt);

  void SetState(const Eigen::VectorXd& Position);
  void SetState(const Eigen::VectorXd& Position, const Eigen::VectorXd& Velocity);
  void SetTarget(const Eigen::VectorXd& target);
  void SetTarget(const Eigen::VectorXd& target, double ReachingTime);
  void SetStateTarget(const Eigen::VectorXd& Position, const Eigen::VectorXd& Target);

  void SetDt(double dt);

  void SetZerkLimits(const Eigen::VectorXd& velLimits);
  void RemoveZerkLimits();
  double GetZerkLimits(unsigned int index);

  void SetAccelLimits(const Eigen::VectorXd& velLimits);
  void RemoveAccelLimits();
  double GetAccelLimits(unsigned int index);

  void SetVelocityLimits(const Eigen::VectorXd& velLimits);
  void RemoveVelocityLimits();
  double GetVelocityLimits(unsigned int index);

  void SetPositionLimits(const Eigen::VectorXd& posLimitsUp, const Eigen::VectorXd& posLimitsDn);
  void RemovePositionLimits();

  void GetTarget(Eigen::VectorXd& target);

  void GetState(Eigen::VectorXd& Position);
  void GetState(Eigen::VectorXd& Position, Eigen::VectorXd& Velocity);
  void GetStateAccel(Eigen::VectorXd& Accel);

  void Update();
  void Update(double dt);

  double GetReachingTime();
};
}// namespace fast_gmm
