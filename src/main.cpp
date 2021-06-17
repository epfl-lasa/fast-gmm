/*
 * main.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: Kim Seungsu
 */
#include <stdio.h>
#include "fast_gmm/GMRDynamics.h"
#include "fast_gmm/MJDynamics.h"

GMRDynamics *GMR_ready;
Gaussians *GMR_dir;
Gaussians *GMR_mag;

void initGMR(void)
{
	GMR_ready = new GMRDynamics(4, 12, 20./1000.,
			"/home/seungsu/data/model/modelLowerPosZ_mu.txt",
			"/home/seungsu/data/model/modelLowerPosZ_sigma.txt",
			"/home/seungsu/data/model/modelLowerPosZ_prio.txt");

	GMR_ready->initGMR(0, 5, 6, 11);
/*
	GMR_dir = new Gaussians(3, 3*2,
			"/home/seungsu/data/models/throwing/GMR_dir_mu.txt",
			"/home/seungsu/data/models/throwing/GMR_dir_sigma.txt",
			"/home/seungsu/data/models/throwing/GMR_dir_prio.txt" );
	GMR_dir->InitFastGMR(0, 2, 3, 5);


	GMR_mag = new Gaussians(3, 4,
			"/home/seungsu/data/models/throwing/GMR_mag_mu.txt",
			"/home/seungsu/data/models/throwing/GMR_mag_sigma.txt",
			"/home/seungsu/data/models/throwing/GMR_mag_prio.txt" );
	GMR_mag->InitFastGMR(0,2,3,3);
	*/
}

void regressionGMR_ready(void)
{
	Vector lTarget(6);
	Vector lState(6), lNextState(6);
	double lTargetDouble[] ={0.127663, -0.214431, -0.043558 , 0.059188, 0.988783, -0.137133};

	lTarget.Set(lTargetDouble, 6);

	GMR_ready->setStateTarget(lTarget, lTarget);

	for(int i=0; i<50; i++){
		lTarget = GMR_ready->getTarget();
		lState = GMR_ready->getState();
		printf("%5.3f %5.3f %5.3f , %5.3f %5.3f %5.3f\n",
				lTarget(0), lTarget(1), lTarget(2),
				lState(0), lState(1), lState(2) );

		lNextState = GMR_ready->getNextState();
		(lState-lTarget).Print("state  - target ");
		lNextState = GMR_ready->getVelocity(lState-lTarget);

		printf("vel : %5.3f %5.3f %5.3f \n", lNextState(0), lNextState(1), lNextState(2) );
	}

}

void testFDyn(void)
{
	MJDynamics *lFDyn;
	MJDynamics *lFDynL;
	Vector lPos(3), lVel(3), lAccel(3), lTarget(3);
	Vector lAccelLimits(3);
	Vector lZerkLimits(3);
	Vector lVelLimits(3);

	lPos.Zero();
	lTarget.One();
	lTarget *= 2.00;

	lFDynL = new MJDynamics(3, 1./500.);
	lFDynL->SetState(lPos);
	lFDynL->SetTarget(lTarget, 0.2);

	lVelLimits.One();
	lVelLimits *= DEG2RAD(190.0);
	lFDynL->SetVelocityLimits(lVelLimits);

	FILE *fid;
	fid =fopen("./log.txt", "w+");
	for(int i=0; i<1000; i++)
	{
		lFDynL->Update();
		lFDynL->GetState(lPos, lVel);
		lFDynL->GetStateAccel(lAccel);

		fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", lPos(0), lPos(1), lPos(2), lVel(0), lVel(1), lVel(2), lAccel(0), lAccel(1), lAccel(2));

//		lFDyn->Update();
//		lFDyn->GetState(lPos, lVel);
//		lFDyn->GetStateAccel(lAccel);
//		fprintf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n", lPos(0), lPos(1), lPos(2), lVel(0), lVel(1), lVel(2), lAccel(0), lAccel(1), lAccel(2));

		if( i==200 )
		{
			lTarget /= 2.00;
			//lFDynL->SetTarget(lTarget, 0.5);
		}
		if( i==400 )
		{
			lTarget /= -10.00;
			//lFDynL->SetTarget(lTarget, 0.5);
		}

	}
	fclose(fid);

}

int main(int argc, char** argv)
{
//	Gaussians *TestModel;
//
//	char gname_mu[1024];
//	char gname_sigma[1024];
//	char gname_prio[1024];
//
//	char dir[] = "/home/seungsu/data/models/";
//	char GraspibilityName[] = "GMM_rest";
//
//	sprintf(gname_mu   , "%s%s_mu.txt"   , dir, GraspibilityName );
//	sprintf(gname_sigma, "%s%s_sigma.txt", dir, GraspibilityName );
//	sprintf(gname_prio , "%s%s_prio.txt" , dir, GraspibilityName );
//
//	Gaussians *mGraspabilityModel;
//	mGraspabilityModel = new Gaussians(gname_mu, gname_sigma, gname_prio);
//	mGraspabilityModel->model.States[1].Sigma.Print();

	initGMR();
	//regressionGMR_ready();

	//testFDyn();

	return 0;
}

