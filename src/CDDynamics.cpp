/*
 * CDDynamics.cpp
 *
 *  Created on: May 29, 2012
 *      Author: Seungsu KIM
 */

/*
 * how to use

	CDDynamics *testDyn;
	testDyn = new CDDynamics(dim, dt, wn);

	testDyn->SetVelocityLimits(velLimits);
	testDyn->SetState(initial);
	testDyn->SetTarget(target);


	start loop
		// if new target is set
		testDyn->SetTarget(target);

		// update dynamics
		testDyn->Update();

		// get state
		testDyn->GetState(state);
	end loop
 */

#include "CDDynamics.h"

#ifdef USE_MATHLIB_NAMESPACE
using namespace MathLib;
#endif

CDDynamics::CDDynamics(int dim, double dt, double Wn)
{
	// set environment
	mDim = dim;
	mDT = dt;
	mWn = Wn;

	// resize variables
	mTarget.Resize(mDim);
	mTargetVelocity.Resize(mDim);
	mState.Resize(mDim);
	mStateVelocity.Resize(mDim);
	mStateAccel.Resize(mDim);

	mPositionLimits.Resize(mDim);
	mVelocityLimits.Resize(mDim);
	mAccelLimits.Resize(mDim);

	// set initial values
	mState.Zero();
	mStateVelocity.Zero();
	mStateAccel.Zero();

	mTarget.Zero();
	mTargetVelocity.Zero();

	mPositionLimits.Zero();
	mVelocityLimits.Zero();
	mAccelLimits.Zero();

	mReachingTime = 0.0;
}

void CDDynamics::SetState(const Vector & Position)
{
	if(mDim == Position.Size()){
		mState = Position - mTarget;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetState() "<<endl;
	}
}

void CDDynamics::SetState(const Vector & Position, const Vector & Velocity)
{
	if( (mDim == Position.Size()) && (mDim == Velocity.Size())){
		mState = Position - mTarget;
		mStateVelocity = Velocity - mTargetVelocity;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetState() "<<endl;
	}
}

void CDDynamics::SetTarget(const Vector & target)
{
	if(mDim == target.Size()){
		mState += (mTarget-target);
		mTarget = target;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetTarget() "<<endl;
	}
}

void CDDynamics::SetTarget(double target[])
{
	for(int i=0; i<mDim; i++)
	{
		mState[i] += (mTarget[i]-target[i]);
		mTarget[i] = target[i];
	}
}

void CDDynamics::SetTarget(const Vector & target, double ReachingTime)
{
	SetTarget(target);
	mReachingTime = ReachingTime;
}


void CDDynamics::SetStateTarget(const Vector & Position, const Vector & Target)
{
	SetTarget(Target);
	SetState(Position);
}

/*
void CDDynamics::SetTarget(const Vector & target, const Vector & targetVel)
{
	if( (mDim == target.Size()) && (mDim == targetVel.Size())){
		mState += (mTarget-target);
		mStateVelocity += (mTargetVelocity-targetVel);

		mTarget = target;
		mTargetVelocity = targetVel;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetTarget() "<<endl;
	}
}
*/

void CDDynamics::SetDt(double dt)
{
	mDT = dt;
}

void CDDynamics::SetWn(double Wn)
{
	mWn = Wn;
}
void CDDynamics::SetAccelLimits(const Vector & accelLimits)
{
	if(mDim == accelLimits.Size()){
		mAccelLimits = accelLimits;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetAccelLimits() "<<endl;
	}
}

void CDDynamics::RemoveAccelLimits(void)
{
	mAccelLimits.Zero();
}

double CDDynamics::GetAccelLimits(unsigned int index)
{
	if( index < mDim ){
		return mAccelLimits(index);
	}
	else{
		return 0.0;
	}
}


void CDDynamics::SetVelocityLimits(const Vector & velLimits)
{
	if(mDim == velLimits.Size()){
		mVelocityLimits = velLimits;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetVelocityLimits() "<<endl;
	}
}

void CDDynamics::RemoveVelocityLimits(void)
{
	mVelocityLimits.Zero();
}

double CDDynamics::GetVelocityLimits(unsigned int index)
{
	if( index < mDim ){
		return mVelocityLimits(index);
	}
	else{
		return 0.0;
	}
}

void CDDynamics::SetPositionLimits(const Vector & posLimits)
{
	if(mDim == posLimits.Size()){
		mPositionLimits = posLimits;
	}
	else{
		cout<<"Dimension error! @ CDDynamics::SetPositionLimits() "<<endl;
	}
}

void CDDynamics::RemovePositionLimits(void)
{
	mPositionLimits.Zero();
}

void CDDynamics::GetTarget(Vector & target)
{
	target = mTarget;
}
/*
void CDDynamics::GetTarget(Vector & target, Vector & targetVel)
{
	target = mTarget;
	targetVel = mTargetVelocity;
}
*/
void CDDynamics::GetState(Vector & Position)
{
	Position = mState + mTarget;
}

void CDDynamics::GetState(double *Position)
{
	for(int i=0; i<mDim; i++)
	{
		Position[i] = mState[i] + mTarget[i];
	}
}

void CDDynamics::GetState(Vector & Position, Vector & Velocity)
{
	Position = mState + mTarget;
	Velocity = mStateVelocity + mTargetVelocity;
}

void CDDynamics::GetStateAccel(Vector & Accel)
{
	Accel = mStateAccel;
}


void CDDynamics::Update()
{
	Update(mDT, 1.0);
}

void CDDynamics::Update(double dt)
{
	Update(dt, 1.0);
}

void CDDynamics::Update(double dt, double muxVel)
{
	// A      = x(0);
	// B      = x_d(0) + w*x(0)
	// x(t)   = (A+Bt)e^(-w*t)
	// x_d(t) = (-w*A+(1-w*t)B ) e^(-w*t)

	double lWn = mWn*muxVel;

	mStateAccel    =  mState*(-lWn*lWn) + mStateVelocity *(-2.0*lWn);

	for(unsigned int i=0; i<mDim; i++ ){
		if( mAccelLimits(i)>0.0){
			if     ( mStateAccel(i) >  mAccelLimits(i) ) mStateAccel(i) =  mAccelLimits(i);
			else if( mStateAccel(i) < -mAccelLimits(i) ) mStateAccel(i) = -mAccelLimits(i);
		}
	}

	mState         += mStateVelocity*dt + mStateAccel*dt*dt*0.5;
	mStateVelocity += mStateAccel*dt;

	for(unsigned int i=0; i<mDim; i++ ){
		if( mVelocityLimits(i)>0){
			if     ( mStateVelocity(i) >  mVelocityLimits(i) ) mStateVelocity(i) =  mVelocityLimits(i);
			else if( mStateVelocity(i) < -mVelocityLimits(i) ) mStateVelocity(i) = -mVelocityLimits(i);
		}
	}


	mReachingTime -= dt;
}


double CDDynamics::GetReachingTime(double dt, double muxVel)
{
	MathLib::Vector lState(mDim);
	MathLib::Vector lStateVelocity(mDim);
	MathLib::Vector lStateAccel(mDim);
	MathLib::Vector lX(mDim);
	MathLib::Vector lB(mDim);
	unsigned int frame;
	double lWn;

	lState = mState;
	lStateVelocity = mStateVelocity;

	lWn = mWn*muxVel;
	for(frame=0; frame<500; frame++){
		lStateAccel    = lState*(-lWn*lWn) + lStateVelocity *(-2.0*lWn);
		for(unsigned int i=0; i<mDim; i++ ){
			if( mAccelLimits(i)>0.0){
				if     ( lStateAccel(i) >  mAccelLimits(i) ) lStateAccel(i) =  mAccelLimits(i);
				else if( lStateAccel(i) < -mAccelLimits(i) ) lStateAccel(i) = -mAccelLimits(i);
			}
		}

		lState         += lStateVelocity*dt;
		lStateVelocity += lStateAccel*dt;

		for(unsigned int i=0; i<mDim; i++ ){
			if( mVelocityLimits(i)>0){
				if     ( lStateVelocity(i) >  mVelocityLimits(i) ) lStateVelocity(i) =  mVelocityLimits(i);
				else if( lStateVelocity(i) < -mVelocityLimits(i) ) lStateVelocity(i) = -mVelocityLimits(i);
			}
		}

		if( lState.Norm() < 0.001 ){
			return (double)frame*dt;
		}
	}
	return (double)frame*dt;
}


double CDDynamics::GetTargetTime(void)
{
	return mReachingTime;
}
