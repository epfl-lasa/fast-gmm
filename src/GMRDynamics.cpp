/*
 * GMRDynamics.cpp
 *
 *  Created on: Nov 20, 2011
 *      Author: Seungsu KIM
 */

#include "fast_gmm/GMRDynamics.h"

#include <iostream>

namespace fast_gmm {

GMRDynamics::GMRDynamics(int nStates,
                         int nVar,
                         double delta_t,
                         const std::vector<double> pri_vec,
                         const std::vector<double> mu_vec,
                         const std::vector<double> sig_vec) {
  this->delta_t = delta_t;
  GMM = new Gaussians(nStates, nVar, pri_vec, mu_vec, sig_vec);
}

void GMRDynamics::initGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex) {
  GMM->InitFastGMR(first_inindex, last_inindex, first_outindex, last_outindex);

  gDim = last_inindex - first_inindex + 1;
  if (gDim != (last_outindex - first_outindex + 1)) {
    std::cout << "dynamics dimension is not matching" << std::endl;
  }

  gXi.resize(gDim);
  target.resize(gDim);

  gXi.setZero();
  target.setZero();
}

void GMRDynamics::setStateTarget(const Eigen::VectorXd& state, const Eigen::VectorXd& input_target) {
  setTarget(input_target);
  setState(state);
}

void GMRDynamics::setTarget(const Eigen::VectorXd& input_target, double input_target_t) {
  this->target_t = input_target_t;

  //gXi += (this->target - target);
  this->target = input_target;
}

Eigen::VectorXd GMRDynamics::getTarget() {
  return target;
}

double GMRDynamics::getTargetT() {
  return target_t;
}

void GMRDynamics::setState(const Eigen::VectorXd& state) {
  gXi = state;
}

Eigen::VectorXd GMRDynamics::getState() {
  return gXi;
}

void GMRDynamics::setCurrentTime(double input_current_t) {
  this->current_t = input_current_t;
}

double GMRDynamics::getCurrentTime() {
  return current_t;
}

Eigen::VectorXd GMRDynamics::getVelocity(const Eigen::VectorXd& x) {
  return GMM->Regression(x);
}

Eigen::VectorXd GMRDynamics::getNextState() {
  return getNextState(1.0);
}

Eigen::VectorXd GMRDynamics::getNextState(double lamda) {
  // target time
  target_t -= (delta_t * lamda);

  gXi += (getVelocity(gXi - target) * (delta_t * lamda));

  return gXi;
}

double GMRDynamics::getReachingTime(double lamda) {
  unsigned int frame = 0;
  unsigned int li = 0;
  Eigen::VectorXd xi(3);
  xi = gXi;

  for (frame = 0; frame < REACHING_ITERATION_MAX; frame++) {
    for (li = 0; li < INTEGRATION_L; li++) {
      xi += getVelocity(xi - target) * delta_t / (double) INTEGRATION_L * lamda;

      if ((xi - target).norm() < GMR_ERROR_TOLERANCE) {
        return (double) (frame * INTEGRATION_L + li) * delta_t / (double) INTEGRATION_L;
      }
    }
  }
  return (double) (frame * INTEGRATION_L + li) * delta_t / (double) INTEGRATION_L;
}
}// namespace fast_gmm
