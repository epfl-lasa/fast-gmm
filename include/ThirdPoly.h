/*
 * ThirdPoly.h
 *
 *  Created on: Oct 30, 2012
 *      Author: seungsu
 */

#ifndef THIRDPOLY_H_
#define THIRDPOLY_H_

#include "MathLib/MathLib.h"

class ThirdPoly
{
private :
	int mDim;

	MathLib::Vector mParamA, mParamB, mParamC, mParamD;

	MathLib::Vector mInitPos, mInitVel;
	MathLib::Vector mEndPos, mEndVel;
	double mDuration;
public :
	ThirdPoly(int dim);

	void SetConstraints(MathLib::Vector &initPos, MathLib::Vector &initVel, MathLib::Vector &endPos, MathLib::Vector &endVel, double duration);

	void Get(double t, MathLib::Vector &pos);
	void Get(double t, MathLib::Vector &pos, MathLib::Vector &vel);
};


#endif /* THIRDPOLY_H_ */
