/*
 * main.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: Kim Seungsu
 */
#include <stdio.h>
#include "GMRDynamics.h"

GMRDynamics *GMR_ready;
Gaussians *GMR_dir;
Gaussians *GMR_mag;

void initGMR(void)
{
	GMR_ready = new GMRDynamics(3, 6, 1./240.,
			"/home/seungsu/data/models/throwing/GMM_ready_mu.txt",
			"/home/seungsu/data/models/throwing/GMM_ready_sigma.txt",
			"/home/seungsu/data/models/throwing/GMM_ready_prio.txt" );
	GMR_ready->initGMR(0, 2, 3, 5);

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
	Vector lTarget(3);
	Vector lState(3), lNextState(3);

	lTarget(0) = 0.01;
	lTarget(1) = 0.01;
	lTarget(2) = 0.01;


	lState(0)= -0.1445;
	lState(1)= -0.0303;
	lState(2)= 0.3139;

	GMR_ready->setTarget(lTarget);
	GMR_ready->setState(lState);

	for(int i=0; i<50; i++){
		lNextState = GMR_ready->getNextState();
		printf("%5.3f %5.3f %5.3f \n", lNextState(0), lNextState(1), lNextState(2) );
	}

}


int main(int argc, char** argv)
{

	initGMR();
	regressionGMR_ready();

	return 0;


}
