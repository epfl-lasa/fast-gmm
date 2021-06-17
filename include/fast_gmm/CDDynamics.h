/*
 * CDDynamics.h
 * Critical Damped 2nd order Dynamics
 *
 *  Created on: May 29, 2012
 *      Author: Seungsu KIM
 */

#ifndef CDDYNAMICS_H_
#define CDDYNAMICS_H_


#include "MathLib.h"

#ifdef USE_MATHLIB_NAMESPACE
using namespace MathLib;
#endif

class CDDynamics
{
private :
	Vector mTarget;
	Vector mTargetVelocity;

	Vector mState;
	Vector mStateVelocity;
	Vector mStateAccel;

	Vector mPositionLimits;
	Vector mVelocityLimits;
	Vector mAccelLimits;

	unsigned int mDim;
	double mWn;
	double mDT;

	double mReachingTime;

public :

	CDDynamics(int dim, double dt, double Wn);

	void SetState(const Vector & Position);
	void SetState(const Vector & Position, const Vector & Velocity);
	void SetTarget(const Vector & target);
	void SetTarget(double target[]);
	void SetTarget(const Vector & target, double ReachingTime);
	void SetStateTarget(const Vector & Position, const Vector & Target);
	//void SetTarget(const Vector & target, const Vector & targetVel);

	void SetDt(double dt);
	void SetWn(double Wn);

	void SetAccelLimits(const Vector & velLimits);
	void RemoveAccelLimits(void);
	double GetAccelLimits(unsigned int index);

	void SetVelocityLimits(const Vector & velLimits);
	void RemoveVelocityLimits(void);
	double GetVelocityLimits(unsigned int index);

	void SetPositionLimits(const Vector & posLimits);
	void RemovePositionLimits(void);


	void GetTarget(Vector & target);
	//void GetTarget(Vector & target, Vector & targetVel);

	void GetState(Vector & Position);
	void GetState(double *Position);
	void GetState(Vector & Position, Vector & Velocity);
	void GetStateAccel(Vector & Accel);

	void Update();
	void Update(double dt);
	void Update(double dt, double muxVel);

	double GetReachingTime(double dt, double muxVel);
	double GetTargetTime(void);
};


#endif /* CDDYNAMICS_H_ */
