/*
 * Gaussians.h
 *
 *  Created on: Nov 19, 2011
 *      Author: Seungsu KIM
 */

#pragma once

#include <eigen3/Eigen/Core>

namespace fast_gmm {

#define GAUSSIAN_MAXIMUM_NUMBER 50

struct GMMState {
  Eigen::VectorXd Mu;
  Eigen::MatrixXd Sigma;
  double Prior;
};

struct GMMStateP {
  Eigen::VectorXd MuI;
  Eigen::MatrixXd SigmaII;
  Eigen::MatrixXd SigmaIIInv;
  double detSigmaII;

  // for GMR
  Eigen::VectorXd muO;
  Eigen::MatrixXd SigmaIO;
  Eigen::MatrixXd SigmaIOInv;
};

struct GMMs {
  unsigned int nbStates;
  unsigned int nbDim;

  GMMState States[GAUSSIAN_MAXIMUM_NUMBER];
};

class Gaussians {
private:
  GMMStateP gmmpinv[GAUSSIAN_MAXIMUM_NUMBER];

public:
  GMMs model;

  Gaussians(const int nbStates,
            const int nbDim,
            const std::vector<double> pri_vec,
            const std::vector<double> mu_vec,
            const std::vector<double> sig_vec);
  Gaussians(GMMs* model);

  void setGMMs(GMMs* model);

  // For fast computation of GaussianPDF
  Eigen::VectorXd gfDiff, gfDiffp;
  Eigen::VectorXd gDer;
  Eigen::VectorXd gPdf;
  int nbDimI;

  void InitFastGaussians(int first_inindex, int last_inindex);
  double GaussianPDFFast(int state, Eigen::VectorXd x);
  double GaussianProbFast(const Eigen::VectorXd& x);
  Eigen::VectorXd GaussianDerProbFast(const Eigen::VectorXd& x);

  void InitFastGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex);
  void Regression(const Eigen::VectorXd& indata, Eigen::VectorXd& outdata, Eigen::MatrixXd& derGMR);
  void Regression(const Eigen::VectorXd& indata, Eigen::VectorXd& outdata);
  Eigen::VectorXd Regression(const Eigen::VectorXd& indata);

};
/*
void GaussianMux(GMMs *modelK, GMMs *modelL, GMMs *modelOut);
void GaussianRotate(GMMs *model, Vector P, Matrix R, GMMs *modelOut);
*/
}// namespace fast_gmm
