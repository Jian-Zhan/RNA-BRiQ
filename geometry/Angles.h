/*
 * Angles.h
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#ifndef GEOMETRY_ANGLES_H_
#define GEOMETRY_ANGLES_H_

#include "geometry/xyz.h"

namespace NSPgeometry {

inline float angleX(const XYZ& A, const XYZ& B){
	float x = A*B;
	float lenA = len(A);
	float lenB = len(B);
	if(lenA == 0 || lenB == 0) return 0.0;
	float cosAng = x/(lenA*lenB);
	if(cosAng > 1)
		cosAng = 1;
	else if(cosAng < -1)
		cosAng = -1;
	return acos(cosAng)*57.2957795;
}

inline float angleX(const XYZ& A, const XYZ& B, const XYZ& C)
{
	return angleX(A-B,C-B);
}


inline float dihedral(const XYZ& A, const XYZ& B, const XYZ& C, const XYZ& D)
{

	XYZ n1 = (A-B)^(C-B);
	XYZ n2 = (B-C)^(D-C);

	float ang = angleX(n1,n2);
	float symbol = (n1^n2)*(C-B);
	if(symbol >= 0)
		return ang;
	else
		return 0-ang;
}

}
#endif /* GEOMETRY_ANGLES_H_ */
