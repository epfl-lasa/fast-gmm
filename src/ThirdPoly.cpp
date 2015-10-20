/*
 * ThridPoly.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: seungsu
 */

#include "ThirdPoly.h"



ThirdPoly::ThirdPoly(int dim)
{
	mDim = dim;

	mInitPos.Resize(mDim);
	mInitVel.Resize(mDim);
	mEndPos.Resize(mDim);
	mEndVel.Resize(mDim);

	mParamA.Resize(mDim);
	mParamB.Resize(mDim);
	mParamC.Resize(mDim);
	mParamD.Resize(mDim);
}


void ThirdPoly::SetConstraints(MathLib::Vector &initPos, MathLib::Vector &initVel, MathLib::Vector &endPos, MathLib::Vector &endVel, double duration)
{
	mInitPos.Set(initPos);
	mInitVel.Set(initVel);
	mEndPos.Set(endPos);
	mEndVel.Set(endVel);

	// calculate polynomial
	mParamA = mInitPos*( 2.0) +mEndPos*(-2.0) +mInitVel       +mEndVel;
	mParamB = mInitPos*(-3.0) +mEndPos*( 3.0) +mInitVel*(-2.0)-mEndVel;
	mParamC = mInitVel;
	mParamD = mInitPos;

	mDuration = duration;
}

void ThirdPoly::Get(double t, MathLib::Vector &pos)
{
	double lx = t/mDuration;

	pos = mParamA*lx*lx*lx + mParamB*lx*lx + mParamC*lx + mParamD;
}

void ThirdPoly::Get(double t, MathLib::Vector &pos, MathLib::Vector &vel)
{
	double lx = t/mDuration;

	pos = mParamA*lx*lx*lx  + mParamB*lx*lx  + mParamC*lx + mParamD;
	vel = mParamA*lx*lx*3.0 + mParamB*lx*2.0 + mParamC;
}

