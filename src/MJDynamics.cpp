/*
 * MJDynamics.cpp
 *
 *  Created on: May 15, 2013
 *      Author: seungsu
 */

#include "fast_gmm/MJDynamics.h"

using namespace MathLib;

MJDynamics::MJDynamics(int dim, double dt)
{
	// set environment
	mDim = dim;
	mDT = dt;

	// resize variables
	mTarget.Resize(mDim);
	mState.Resize(mDim);
	mStateVelocity.Resize(mDim);
	mStateAccel.Resize(mDim);
	mStateZerk.Resize(mDim);

	mPositionLimitsUp.Resize(mDim);
	mPositionLimitsDn.Resize(mDim);
	mVelocityLimits.Resize(mDim);
	mAccelLimits.Resize(mDim);
	mZerkLimits.Resize(mDim);

	a0.Resize(mDim);
	a1.Resize(mDim);
	a2.Resize(mDim);
	a3.Resize(mDim);
	a4.Resize(mDim);
	a5.Resize(mDim);

	// set initial values
	mState.Zero();
	mStateVelocity.Zero();
	mStateAccel.Zero();
	mStateZerk.Zero();

	mTarget.Zero();

	mPositionLimitsUp.Zero();
	mPositionLimitsDn.Zero();
	mAccelLimits.Zero();
	mZerkLimits.Zero();

	a0.Zero();
	a1.Zero();
	a2.Zero();
	a3.Zero();
	a4.Zero();
	a5.Zero();

	mReachingTime = 0.0;
}

void MJDynamics::SetState(const Vector & Position)
{
	if(mDim == Position.Size()){
		mState = Position - mTarget;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetState() "<<endl;
	}
}

void MJDynamics::SetState(const Vector & Position, const Vector & Velocity)
{
	if( (mDim == Position.Size()) && (mDim == Velocity.Size())){
		mState = Position;
		mStateVelocity = Velocity;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetState() "<<endl;
	}
}

void MJDynamics::SetTarget(const Vector & target)
{
	if(mDim == target.Size()){
		mTarget = target;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetTarget() "<<endl;
	}
}

void MJDynamics::SetTarget(const Vector & target, double ReachingTime)
{
	if(mDim == target.Size()){
		mTarget = target;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetTarget() "<<endl;
	}

	mReachingTime = ReachingTime;
	mCurrentTime = 0.0;
}


void MJDynamics::SetStateTarget(const Vector & Position, const Vector & Target)
{
	SetTarget(Target);
	SetState(Position);
}

void MJDynamics::SetDt(double dt)
{
	mDT = dt;
}


void MJDynamics::SetAccelLimits(const Vector & accelLimits)
{
	if(mDim == accelLimits.Size()){
		mAccelLimits = accelLimits;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetAccelLimits() "<<endl;
	}
}

void MJDynamics::RemoveAccelLimits(void)
{
	mAccelLimits.Zero();
}

double MJDynamics::GetAccelLimits(unsigned int index)
{
	if( index < mDim ){
		return mAccelLimits(index);
	}
	else{
		return 0.0;
	}
}


void MJDynamics::SetZerkLimits(const Vector & zerkLimits)
{
	if(mDim == zerkLimits.Size()){
		mZerkLimits.Set(zerkLimits);
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetZerkLimits() "<<endl;
	}
}

void MJDynamics::RemoveZerkLimits(void)
{
	mZerkLimits.Zero();
}

double MJDynamics::GetZerkLimits(unsigned int index)
{
	if( index < mDim ){
		return mZerkLimits(index);
	}
	else{
		return 0.0;
	}
}

void MJDynamics::SetVelocityLimits(const Vector & velLimits)
{
	if(mDim == velLimits.Size()){
		mVelocityLimits = velLimits;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetVelocityLimits() "<<endl;
	}
}

void MJDynamics::RemoveVelocityLimits(void)
{
	mVelocityLimits.Zero();
}

double MJDynamics::GetVelocityLimits(unsigned int index)
{
	if( index < mDim ){
		return mVelocityLimits(index);
	}
	else{
		return 0.0;
	}
}

void MJDynamics::SetPositionLimits(const Vector & posLimitsUp, const Vector & posLimitsDn)
{
	if( (mDim == posLimitsUp.Size()) && (mDim == posLimitsDn.Size()) ){
		mPositionLimitsUp = posLimitsUp;
		mPositionLimitsDn = posLimitsDn;
	}
	else{
		cout<<"Dimension error! @ MJDynamics::SetPositionLimits() "<<endl;
	}
}

void MJDynamics::RemovePositionLimits(void)
{
	mPositionLimitsUp.Zero();
	mPositionLimitsDn.Zero();
}

void MJDynamics::GetTarget(Vector & target)
{
	target = mTarget;
}

void MJDynamics::GetState(Vector & Position)
{
	Position = mState;
}

void MJDynamics::GetState(Vector & Position, Vector & Velocity)
{
	Position = mState;
	Velocity = mStateVelocity;
}

void MJDynamics::GetStateAccel(Vector & Accel)
{
	Accel = mStateAccel;
}


void MJDynamics::Update()
{
	Update(mDT);

}

void MJDynamics::Update(double dt)
{
	double lR2, lR3;
	double lD2, lD3, lD4, lD5;

	if( mReachingTime < dt)
	{
		if( ((mTarget-mState).Norm() > 0.0001) || (mStateVelocity.Norm() > 0.001) ){
			mReachingTime = dt*2.0;
		}
	}

	if(mReachingTime >= dt){
		lR2 = mReachingTime*mReachingTime;
		lR3 = lR2 * mReachingTime;
		lD2 = dt*dt;
		lD3 = lD2*dt;
		lD4 = lD3*dt;
		lD5 = lD4*dt;

		a0 = mState;
		a1 = mStateVelocity;
		a2 = mStateAccel/2.0;

		a3 = (mTarget*10.0 - a0*10.0   - a1*6.0*mReachingTime - a2*3.0*lR2)/(     lR3);
		a4 = (               a1*(-2.0) - a2*3.0*mReachingTime - a3*3.0*lR2)/( 2.0*lR3);
		a5 = (               a2*(-1.0) - a3*3.0*mReachingTime - a4*6.0*lR2)/(10.0*lR3);

		mState         = a0 +a1*dt +a2*lD2    +a3*lD3     +a4*lD4      +a5*lD5;
		mStateVelocity =     a1    +a2*2.0*dt +a3*3.0*lD2 +a4*4.0*lD3  +a5*5.0*lD4;
		mStateAccel    =            a2*2.0    +a3*6.0*dt  +a4*12.0*lD2 +a5*20.0*lD3;

		mReachingTime -= dt;
	}
}


double MJDynamics::GetReachingTime(void)
{
	return mReachingTime;
}

