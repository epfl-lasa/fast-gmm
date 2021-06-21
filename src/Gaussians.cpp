/*
 * Gaussians.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: Seungsu KIM
 */

#include "fast_gmm/Gaussians.h"

#include <cmath>
#include <iostream>

namespace fast_gmm {
/*
Gaussians::Gaussians(GMMs *model)
{
	this->model.nbStates = model->nbStates;
	this->model.nbDim    = model->nbDim;

	this->model.States = (GMMState  *)malloc(model->nbStates*sizeof(GMMState) );

	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}
*/
Gaussians::Gaussians(const int nbStates,
                     const int nbDim,
                     const std::vector<double> pri_vec,
                     const std::vector<double> mu_vec,
                     const std::vector<double> sig_vec) {

  model.nbStates = nbStates;
  model.nbDim = nbDim;

  for (int s = 0; s < nbStates; s++) {
    model.States[s].Mu.resize(nbDim);
    model.States[s].Sigma.resize(nbDim, nbDim);
  }

  for (int s = 0; s < nbStates; s++) {
    model.States[s].Prior = pri_vec[s];
  }
  // cout << endl << "Printing the constructed Priors" << endl;
  // for ( int s = 0; s < nbStates; s++ ) {
  // 	cout << model.States[s].Prio  << "\t";
  // }
  // cout << endl;

  for (int s = 0; s < nbStates; s++) {
    for (int d = 0; d < nbDim; d++) {
      model.States[s].Mu[d] = mu_vec[s * nbDim + d];
    }
  }

  // cout << endl << "Printing the constructed Mu" << endl;
  // for ( int s = 0; s < nbStates; s++ ) {
  // 	for (int d = 0; d < nbDim; d++) {
  // 		cout << model.States[s].Mu[d]  << "\t";
  // 	}
  // 	cout << endl;
  // }

  for (int s = 0; s < nbStates; s++) {
    for (int row = 0; row < nbDim; row++) {
      for (int col = 0; col < nbDim; col++) {
        int ind = s * nbDim * nbDim + row * nbDim + col;
        model.States[s].Sigma(row, col) = sig_vec[ind];
      }
    }
  }

  // cout << endl << "Printing the constructed Sigma" << endl;
  // for ( int s = 0; s < nbStates; s++ ) {
  // 	for (int row = 0; row < nbDim; row++) {
  // 		for (int col = 0; col < nbDim; col++) {
  // 			cout << model.States[s].Sigma(row, col) << "\t";
  // 		}
  // 		cout <<endl;
  // 	}
  // 	cout << endl;
  // }
}

void Gaussians::setGMMs(GMMs* model) {
  for (unsigned int s = 0; s < model->nbStates; s++) {
    this->model.States[s].Mu = model->States[s].Mu;
    this->model.States[s].Sigma = model->States[s].Sigma;
    this->model.States[s].Prior = model->States[s].Prior;
  }
}

void Gaussians::InitFastGaussians(int first_inindex, int last_inindex) {
  double det;
  int nbIN = last_inindex - first_inindex + 1;

  for (unsigned int s = 0; s < model.nbStates; s++) {
    gmmpinv[s].MuI.resize(nbIN);
    gmmpinv[s].SigmaII.resize(nbIN, nbIN);
    gmmpinv[s].SigmaIIInv.resize(nbIN, nbIN);
  }

  for (unsigned int s = 0; s < model.nbStates; s++) {
    for (int i = first_inindex; i <= last_inindex; i++) { gmmpinv[s].MuI(i - first_inindex) = model.States[s].Mu(i); }
    for (int i = first_inindex; i <= last_inindex; i++) {
      for (int j = first_inindex; j <= last_inindex; j++) {
        gmmpinv[s].SigmaII(i - first_inindex, j - first_inindex) = model.States[s].Sigma(i, j);
      }
    }

    gmmpinv[s].SigmaIIInv = gmmpinv[s].SigmaII.inverse();
    double det = gmmpinv[s].SigmaIIInv.determinant();
    if (det < 0) { det = 1e-30; }
    gmmpinv[s].detSigmaII = det;
  }

  nbDimI = last_inindex - first_inindex + 1;
  gfDiff.resize(nbDimI);
  gfDiffp.resize(nbDimI);
  gDer.resize(nbDimI);

}

double Gaussians::GaussianPDFFast(int state, Eigen::VectorXd x) {
  double p;
  gfDiff = x - gmmpinv[state].MuI;
  gfDiffp = gmmpinv[state].SigmaIIInv * gfDiff;

  p = exp(-0.5 * gfDiff.dot(gfDiffp)) / sqrt(pow(2.0 * M_PI, nbDimI) * (gmmpinv[state].detSigmaII + 1e-30));

  return p;
}

double Gaussians::GaussianProbFast(const Eigen::VectorXd& x) {
  double totalP = 0;
  for (unsigned int s = 0; s < model.nbStates; s++) {
    totalP += model.States[s].Prior * GaussianPDFFast(s, x);
  }
  return totalP;
}

Eigen::VectorXd Gaussians::GaussianDerProbFast(const Eigen::VectorXd& x) {
  gDer.setZero();
  for (unsigned int s = 0; s < model.nbStates; s++) {
    gDer += (gmmpinv[s].SigmaIIInv * (x - gmmpinv[s].MuI)) * model.States[s].Prior * GaussianPDFFast(s, x);
  }
  return -gDer;
}

void Gaussians::InitFastGMR(int first_inindex, int last_inindex, int first_outindex, int last_outindex) {
  double det;
  int nbIN = last_inindex - first_inindex + 1;
  int nbOUT = last_outindex - first_outindex + 1;

  gPdf.resize(model.nbStates);

  for (unsigned int s = 0; s < model.nbStates; s++) {
    gmmpinv[s].MuI.resize(nbIN);
    gmmpinv[s].SigmaII.resize(nbIN, nbIN);
    gmmpinv[s].SigmaIIInv.resize(nbIN, nbIN);

    gmmpinv[s].muO.resize(nbOUT);
    gmmpinv[s].SigmaIO.resize(nbIN, nbOUT);
    gmmpinv[s].SigmaIOInv.resize(nbOUT, nbOUT);
  }

  for (unsigned int s = 0; s < model.nbStates; s++) {
    for (int i = first_inindex; i <= last_inindex; i++) {
      gmmpinv[s].MuI(i - first_inindex) = model.States[s].Mu(i);

      for (int j = first_inindex; j <= last_inindex; j++) {
        gmmpinv[s].SigmaII(i - first_inindex, j - first_inindex) = model.States[s].Sigma(i, j);
      }
      for (int j = first_outindex; j <= last_outindex; j++) {
        gmmpinv[s].SigmaIO(i - first_inindex, j - first_outindex) = model.States[s].Sigma(i, j);
      }
    }

    for (int i = first_outindex; i <= last_outindex; i++) {
      gmmpinv[s].muO(i - first_outindex) = model.States[s].Mu(i);
    }

    gmmpinv[s].SigmaIIInv = gmmpinv[s].SigmaII.inverse();
    det = gmmpinv[s].SigmaIIInv.determinant();
    if (det < 0) { det = 1e-30; }
    gmmpinv[s].detSigmaII = det;
    gmmpinv[s].SigmaIOInv = gmmpinv[s].SigmaIO.transpose().inverse();
  }

  nbDimI = last_inindex - first_inindex + 1;
  gfDiff.resize(nbDimI);
  gfDiffp.resize(nbDimI);
  gDer.resize(nbDimI);

}

void Gaussians::Regression(const Eigen::VectorXd& indata, Eigen::VectorXd& outdata, Eigen::MatrixXd& derGMR) {
  Regression(indata, outdata);
  std::cout << "derivative is not implemented yet " << std::endl;
}

void Gaussians::Regression(const Eigen::VectorXd& indata, Eigen::VectorXd& outdata) {
  double pdfall;
  Eigen::VectorXd h(model.nbStates);
  Eigen::VectorXd r_diff(outdata.size());

  for (unsigned int s = 0; s < model.nbStates; s++) {
    gPdf(s) = model.States[s].Prior * GaussianPDFFast(s, indata);
  }
  pdfall = gPdf.sum();

  outdata.setZero();
  for (unsigned int s = 0; s < model.nbStates; s++) {
    //h(s) = gPdf(s)/(pdfall + 1e-30 );
    h(s) = gPdf(s) / (pdfall);
    r_diff = gmmpinv[s].SigmaIO.transpose() * gmmpinv[s].SigmaIIInv * (indata - gmmpinv[s].MuI);

    for (unsigned int i = 0; i < r_diff.size(); i++) {
      outdata(i) += h(s) * (r_diff(i) + gmmpinv[s].muO(i));
    }
  }
}

Eigen::VectorXd Gaussians::Regression(const Eigen::VectorXd& indata) {
  Eigen::VectorXd outdata(indata.size());
  Regression(indata, outdata);
  return outdata;
}


/*
#include <math.h>
#include "Gaussians.h"

#include "armadillo"

using namespace arma;
using namespace std;

Gaussians::Gaussians(GMMs *model)
{
	this->model.nbStates = model->nbStates;
	this->model.nbDim    = model->nbDim;

	this->model.States = (GMMState  *)malloc(model->nbStates*sizeof(GMMState) );

	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}

Gaussians::Gaussians(int nbStates, int nbDim, char *f_mu, char *f_sigma, char *f_prio)
{

	int s, i, j;

	model.nbStates = nbStates;
	model.nbDim    = nbDim;
	model.States = (GMMState  *)malloc(nbStates*sizeof(GMMState) );

	for( s=0; s<nbStates; s++ ){
		model.States[s].Mu       =  zeros<vec>(nbDim);
		model.States[s].Sigma    =  zeros<mat>(nbDim, nbDim );
	}

	// f_mu
	ifstream fin;
	fin.open(f_mu);
	for( i=0; i<nbDim; i++ ){
		for( s=0; s<nbStates; s++ ){
			fin >> model.States[s].Mu(i);
		}
	}
	fin.close();

	// f_sigma
	fin.open(f_sigma);
	for( s=0; s<nbStates; s++ ){
		for( i=0; i<nbDim; i++ ){
			for( j=0; j<nbDim; j++ ){
				fin >> model.States[s].Sigma(i,j);
			}
		}
	}
	fin.close();

	// f_prio
	fin.open(f_prio);
	for( s=0; s<nbStates; s++ ){
		fin >>model.States[s].Prio;
	}
	fin.close();
}

void Gaussians::setGMMs(GMMs *model)
{
	for(int s=0; s<model->nbStates; s++ ){
		this->model.States[s].Mu    = model->GMMState[s].Mu;
		this->model.States[s].Sigma = model->GMMState[s].Sigma;
		this->model.States[s].Prio  = model->GMMState[s].Prio;
	}
}


void Gaussians::InitFastGaussians(int first_inindex, int last_inindex)
{
	gmmpinv = (GMMStateP *)malloc(model.nbStates*sizeof(GMMStateP) );

	for(int s=0; s<model.nbStates; s++ ){
		gmmpinv[s].MuI = model.States[s].Mu.subvec(first_inindex, last_inindex);
		gmmpinv[s].SigmaII = model.States[s].Sigma.submat(first_inindex, first_inindex, last_inindex, last_inindex);
		gmmpinv[s].SigmaIIInv = pinv(gmmpinv[s].SigmaII);
		gmmpinv[s].detSigmaII = det(gmmpinv[s].SigmaII);
	}

	nbDimI = last_inindex - first_inindex +1;
	gfDiff  = zeros<vec>(nbDimI);
	gfDiffp = zeros<vec>(nbDimI);
	gDer    = zeros<vec>(nbDimI);
}

double Gaussians::GaussianPDFFast(int state, vec x)
{
	double p;
	gfDiff  = x - gmmpinv[state].MuI;
	gfDiffp = gmmpinv[state].SigmaIIInv * gfDiff;

	p = exp(-0.5*dot(gfDiff, gfDiffp)) / sqrt(pow(2.0*math::pi(), nbDimI)*( gmmpinv[state].detSigmaII +1e-30));

	return p;
}

double Gaussians::GaussianProbFast(vec x)
{
	double totalP=0;
	for(int s=0; s<model.nbStates; s++ ){
		totalP += model.States[s].Prio*GaussianPDFFast(s,x);
	}
	return totalP;
}

vec Gaussians::GaussianDerProbFast(vec x)
{
	gDer.zeros();
	for(int s=0; s<model.nbStates; s++ ){
		gDer += model.States[s].Prio * gmmpinv[s].SigmaIIInv *(x-gmmpinv[s].MuI)*GaussianPDFFast(s,x);
	}
	return -gDer;
}

//-------------------------------------------------------------------------------------------------------
void AllocateGMMs(GMMs *model, int nbDim, int nbStates)
{
	model->nbDim = nbDim;
	model->nbStates = nbStates;
	model->GMMState = (GMMState  *)malloc(nbStates*sizeof(GMMState) );

	for(int s=0; s<nbStates; s++ ){
		model->GMMState[s].Mu       =  zeros<vec>(nbDim);
		model->GMMState[s].Sigma    =  zeros<mat>(nbDim, nbDim );
	}
}


double GaussianPDF(vec x, vec Mu, mat Sigma)
{
	double p;
	vec diff  = x - Mu;
	vec diffp = pinv( Sigma ) * diff;
	int nbDim = x.size();

	p = exp(-0.5*dot(diff, diffp)) / sqrt(pow(2.0*math::pi(), nbDim)*( abs(det(Sigma)) +1e-30));

    if(p < 1e-30){
		return 1e-30;
    }
	else{
		return p;
	}
}

void GaussianMux(GMMs *modelK, GMMs *modelL, GMMs *modelOut)
{
	int k,l,j;
	int K = modelK->nbStates;
	int L = modelL->nbStates;
	int J = K*L;

	//modelOut->nbDim = modelK->nbDim;
	//modelOut->nbStates = J;
	//modelOut->GMMState = (GMMState *)malloc(J*sizeof(GMMState) );

	j=0;
	for(k=0; k<K; k++){
		for(l=0; l<L; l++){
			modelOut->GMMState[j].Sigma = pinv( pinv(modelK->GMMState[k].Sigma) + pinv( modelL->GMMState[l].Sigma) );
			modelOut->GMMState[j].Mu    = modelOut->GMMState[j].Sigma *( pinv(modelK->GMMState[k].Sigma) * modelK->GMMState[k].Mu + pinv(modelL->GMMState[l].Sigma) * modelL->GMMState[l].Mu );
			modelOut->GMMState[j].Prio  = modelK->GMMState[k].Prio* modelL->GMMState[l].Prio * GaussianPDF( modelK->GMMState[k].Mu, modelL->GMMState[l].Mu, modelK->GMMState[k].Sigma + modelL->GMMState[l].Sigma);
			j++;
		}
	}
}

void GaussianRotate(GMMs *model, arma::vec P, arma::mat R, GMMs *modelOut)
{
	for(int s=0; s<model->nbStates; s++){
		modelOut->GMMState[s].Mu    = R*model->GMMState[s].Mu + P;
		modelOut->GMMState[s].Sigma = R*model->GMMState[s].Sigma*trans(R);
	}
	//modelOut->nbDim = model->nbDim;
	//modelOut->nbStates = model->nbStates;
}
*/
}// namespace fast_gmm