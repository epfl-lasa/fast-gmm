/*
 * MJDynamics.h
 *
 *  Created on: May 15, 2013
 *      Author: seungsu
 */

#ifndef MZDYNAMICS_H_
#define MZDYNAMICS_H_

#include "MathLib.h"

using namespace MathLib;

class MJDynamics
{
private :
	Vector mTarget;

	Vector mState;
	Vector mStateVelocity;
	Vector mStateAccel;
	Vector mStateZerk;

	Vector mPositionLimitsUp;
	Vector mPositionLimitsDn;
	Vector mVelocityLimits;
	Vector mAccelLimits;
	Vector mZerkLimits;

	Vector a0, a1, a2, a3, a4, a5;

	unsigned int mDim;
	double mDT;

	double mReachingTime;
	double mCurrentTime;
public :

	MJDynamics(int dim, double dt);

	void SetState(const Vector & Position);
	void SetState(const Vector & Position, const Vector & Velocity);
	void SetTarget(const Vector & target);
	void SetTarget(const Vector & target, double ReachingTime);
	void SetStateTarget(const Vector & Position, const Vector & Target);

	void SetDt(double dt);

	void SetZerkLimits(const Vector & velLimits);
	void RemoveZerkLimits(void);
	double GetZerkLimits(unsigned int index);

	void SetAccelLimits(const Vector & velLimits);
	void RemoveAccelLimits(void);
	double GetAccelLimits(unsigned int index);

	void SetVelocityLimits(const Vector & velLimits);
	void RemoveVelocityLimits(void);
	double GetVelocityLimits(unsigned int index);

	void SetPositionLimits(const Vector & posLimitsUp, const Vector & posLimitsDn);
	void RemovePositionLimits(void);

	void GetTarget(Vector & target);

	void GetState(Vector & Position);
	void GetState(Vector & Position, Vector & Velocity);
	void GetStateAccel(Vector & Accel);

	void Update();
	void Update(double dt);

	double GetReachingTime(void);
};

#endif /* MZDYNAMICS_H_ */
