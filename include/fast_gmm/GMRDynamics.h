/*
 * GMRDynamics.h
 *
 *  Created on: Nov 20, 2011
 *      Author: Seungsu KIM
 */

#pragma once

#include <vector>
#include "fast_gmm/Gaussians.h"

namespace fast_gmm {

#define GMR_ERROR_TOLERANCE 0.001
#define INTEGRATION_L 5
#define REACHING_ITERATION_MAX 500
#define REAL_MIN (1e-30)

// GMR Dynamics
class GMRDynamics {
private:
  Gaussians* GMM;

  double delta_t;
  double target_t;
  double current_t;

  Eigen::VectorXd gXi;
  Eigen::VectorXd target;
  unsigned int gDim;

public:
  GMRDynamics(int nStates,
              int nVar,
              double delta_t,
              const std::vector<double> pri_vec,
              const std::vector<double> mu_vec,
              const std::vector<double> sig_vec);

  void initGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex);

  void setStateTarget(const Eigen::VectorXd& state, const Eigen::VectorXd& input_target);
  void setTarget(const Eigen::VectorXd& input_target, double input_target_t = -1.0);
  Eigen::VectorXd getTarget();
  double getTargetT();
  void setState(const Eigen::VectorXd& state);
  Eigen::VectorXd getState();
  void setCurrentTime(double current_t);
  double getCurrentTime();

  Eigen::VectorXd getVelocity(const Eigen::VectorXd& x);

  Eigen::VectorXd getNextState();
  Eigen::VectorXd getNextState(double lamda);
  double getReachingTime(double lamda);
};
}// namespace fast_gmm
