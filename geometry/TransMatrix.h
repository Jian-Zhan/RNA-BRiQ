/*
 * TransMatrix.h
 *
 *  Created on: Nov 1, 2018
 *      Author: s2982206
 */

#ifndef GEOMETRY_TRANSMATRIX_H_
#define GEOMETRY_TRANSMATRIX_H_

#include "geometry/xyz.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

namespace NSPgeometry {
using namespace std;
class TransMatrix {
public:
	double mtx[3][3];
	TransMatrix() {
		mtx[0][0] = 1;
		mtx[0][1] = 0;
		mtx[0][2] = 0;
		mtx[1][0] = 0;
		mtx[1][1] = 1;
		mtx[1][2] = 0;
		mtx[2][0] = 0;
		mtx[2][1] = 0;
		mtx[2][2] = 1;
	}

	TransMatrix(double values[][3]){
		mtx[0][0] = values[0][0];
		mtx[0][1] = values[0][1];
		mtx[0][2] = values[0][2];
		mtx[1][0] = values[1][0];
		mtx[1][1] = values[1][1];
		mtx[1][2] = values[1][2];
		mtx[2][0] = values[2][0];
		mtx[2][1] = values[2][1];
		mtx[2][2] = values[2][2];
	}

	TransMatrix(const XYZ& x, const XYZ& y, const XYZ& z){
		mtx[0][0] = x.x_;
		mtx[0][1] = y.x_;
		mtx[0][2] = z.x_;
		mtx[1][0] = x.y_;
		mtx[1][1] = y.y_;
		mtx[1][2] = z.y_;
		mtx[2][0] = x.z_;
		mtx[2][1] = y.z_;
		mtx[2][2] = z.z_;
	}

	TransMatrix& operator=(const TransMatrix& other);

	void normalize(){
		XYZ x(mtx[0][0], mtx[1][0], mtx[2][0]);
		XYZ y(mtx[0][1], mtx[1][1], mtx[2][1]);
		XYZ z(mtx[0][2], mtx[1][2], mtx[2][2]);
		x = ~x;
		y = z^x;
		y = ~y;
		z = x^y;

		mtx[0][0] = x.x_;
		mtx[0][1] = y.x_;
		mtx[0][2] = z.x_;
		mtx[1][0] = x.y_;
		mtx[1][1] = y.y_;
		mtx[1][2] = z.y_;
		mtx[2][0] = x.z_;
		mtx[2][1] = y.z_;
		mtx[2][2] = z.z_;
	}

	void copyValue(const TransMatrix& other){
		mtx[0][0] = other.mtx[0][0];
		mtx[0][1] = other.mtx[0][1];
		mtx[0][2] = other.mtx[0][2];
		mtx[1][0] = other.mtx[1][0];
		mtx[1][1] = other.mtx[1][1];
		mtx[1][2] = other.mtx[1][2];
		mtx[2][0] = other.mtx[2][0];
		mtx[2][1] = other.mtx[2][1];
		mtx[2][2] = other.mtx[2][2];
	}

	XYZ applyTransition(const XYZ& coord);

	TransMatrix transpose(){
		TransMatrix out;
		out.mtx[0][0] = mtx[0][0];
		out.mtx[0][1] = mtx[1][0];
		out.mtx[0][2] = mtx[2][0];
		out.mtx[1][0] = mtx[0][1];
		out.mtx[1][1] = mtx[1][1];
		out.mtx[1][2] = mtx[2][1];
		out.mtx[2][0] = mtx[0][2];
		out.mtx[2][1] = mtx[1][2];
		out.mtx[2][2] = mtx[2][2];
		return out;
	}

	TransMatrix multiply(TransMatrix& b);

	TransMatrix rotByX(double theta);
	TransMatrix rotByY(double theta);
	TransMatrix rotByZ(double theta);
	TransMatrix randomRot(double theta);

	void checkTM() {
		XYZ x(mtx[0][0], mtx[1][0], mtx[2][0]);
		XYZ y(mtx[0][1], mtx[1][1], mtx[2][1]);
		XYZ z(mtx[0][2], mtx[1][2], mtx[2][2]);
		printf("%11.8f\n", x.length());
		printf("%11.8f\n", y.length());
		printf("%11.8f\n", z.length());


		XYZ xx = z^x;

		printf("%10.7f %10.7f %10.7f\n", xx.x_, xx.y_, xx.z_);
		printf("%10.7f %10.7f %10.7f\n", y.x_, y.y_, y.z_);

	}

	void print() const{
		printf("%6.3f %6.3f %6.3f\n",mtx[0][0], mtx[0][1], mtx[0][2]);
		printf("%6.3f %6.3f %6.3f\n",mtx[1][0], mtx[1][1], mtx[1][2]);
		printf("%6.3f %6.3f %6.3f\n",mtx[2][0], mtx[2][1], mtx[2][2]);

	}

	virtual ~TransMatrix();
};

inline XYZ operator *(const XYZ& coord, const TransMatrix& tm) {
	return XYZ(tm.mtx[0][0]*coord.x_ + tm.mtx[0][1]*coord.y_ + tm.mtx[0][2]*coord.z_,
			tm.mtx[1][0]*coord.x_ + tm.mtx[1][1]*coord.y_ + tm.mtx[1][2]*coord.z_,
			tm.mtx[2][0]*coord.x_ + tm.mtx[2][1]*coord.y_ + tm.mtx[2][2]*coord.z_);
}

inline TransMatrix operator *(const TransMatrix& tmA, const TransMatrix& tmB) {
	double mtx[3][3];
	mtx[0][0] = tmA.mtx[0][0]*tmB.mtx[0][0] + tmA.mtx[0][1]*tmB.mtx[1][0] + tmA.mtx[0][2]*tmB.mtx[2][0];
	mtx[0][1] = tmA.mtx[0][0]*tmB.mtx[0][1] + tmA.mtx[0][1]*tmB.mtx[1][1] + tmA.mtx[0][2]*tmB.mtx[2][1];
	mtx[0][2] = tmA.mtx[0][0]*tmB.mtx[0][2] + tmA.mtx[0][1]*tmB.mtx[1][2] + tmA.mtx[0][2]*tmB.mtx[2][2];
	mtx[1][0] = tmA.mtx[1][0]*tmB.mtx[0][0] + tmA.mtx[1][1]*tmB.mtx[1][0] + tmA.mtx[1][2]*tmB.mtx[2][0];
	mtx[1][1] = tmA.mtx[1][0]*tmB.mtx[0][1] + tmA.mtx[1][1]*tmB.mtx[1][1] + tmA.mtx[1][2]*tmB.mtx[2][1];
	mtx[1][2] = tmA.mtx[1][0]*tmB.mtx[0][2] + tmA.mtx[1][1]*tmB.mtx[1][2] + tmA.mtx[1][2]*tmB.mtx[2][2];
	mtx[2][0] = tmA.mtx[2][0]*tmB.mtx[0][0] + tmA.mtx[2][1]*tmB.mtx[1][0] + tmA.mtx[2][2]*tmB.mtx[2][0];
	mtx[2][1] = tmA.mtx[2][0]*tmB.mtx[0][1] + tmA.mtx[2][1]*tmB.mtx[1][1] + tmA.mtx[2][2]*tmB.mtx[2][1];
	mtx[2][2] = tmA.mtx[2][0]*tmB.mtx[0][2] + tmA.mtx[2][1]*tmB.mtx[1][2] + tmA.mtx[2][2]*tmB.mtx[2][2];
	return TransMatrix(mtx);
}

inline TransMatrix operator !(const TransMatrix& tm) {
	/*
	 * transpose
	 */
	double mtx[3][3];
	mtx[0][0] = tm.mtx[0][0];
	mtx[0][1] = tm.mtx[1][0];
	mtx[0][2] = tm.mtx[2][0];
	mtx[1][0] = tm.mtx[0][1];
	mtx[1][1] = tm.mtx[1][1];
	mtx[1][2] = tm.mtx[2][1];
	mtx[2][0] = tm.mtx[0][2];
	mtx[2][1] = tm.mtx[1][2];
	mtx[2][2] = tm.mtx[2][2];
	return TransMatrix(mtx);
}

inline XYZ TransMatrix::applyTransition(const XYZ& coord) {
	return XYZ(mtx[0][0]*coord.x_ + mtx[0][1]*coord.y_ + mtx[0][2]*coord.z_,
			mtx[1][0]*coord.x_ + mtx[1][1]*coord.y_ + mtx[1][2]*coord.z_,
			mtx[2][0]*coord.x_ + mtx[2][1]*coord.y_ + mtx[2][2]*coord.z_);
}

inline TransMatrix& TransMatrix::operator =(const TransMatrix& other){
	mtx[0][0] = other.mtx[0][0];
	mtx[0][1] = other.mtx[0][1];
	mtx[0][2] = other.mtx[0][2];
	mtx[1][0] = other.mtx[1][0];
	mtx[1][1] = other.mtx[1][1];
	mtx[1][2] = other.mtx[1][2];
	mtx[2][0] = other.mtx[2][0];
	mtx[2][1] = other.mtx[2][1];
	mtx[2][2] = other.mtx[2][2];
	return *this;
}

inline TransMatrix TransMatrix::multiply(TransMatrix& b){
	TransMatrix out;
	out.mtx[0][0] = mtx[0][0]*b.mtx[0][0] + mtx[0][1]*b.mtx[1][0] + mtx[0][2]*b.mtx[2][0];
	out.mtx[0][1] = mtx[0][0]*b.mtx[0][1] + mtx[0][1]*b.mtx[1][1] + mtx[0][2]*b.mtx[2][1];
	out.mtx[0][2] = mtx[0][0]*b.mtx[0][2] + mtx[0][1]*b.mtx[1][2] + mtx[0][2]*b.mtx[2][2];
	out.mtx[1][0] = mtx[1][0]*b.mtx[0][0] + mtx[1][1]*b.mtx[1][0] + mtx[1][2]*b.mtx[2][0];
	out.mtx[1][1] = mtx[1][0]*b.mtx[0][1] + mtx[1][1]*b.mtx[1][1] + mtx[1][2]*b.mtx[2][1];
	out.mtx[1][2] = mtx[1][0]*b.mtx[0][2] + mtx[1][1]*b.mtx[1][2] + mtx[1][2]*b.mtx[2][2];
	out.mtx[2][0] = mtx[2][0]*b.mtx[0][0] + mtx[2][1]*b.mtx[1][0] + mtx[2][2]*b.mtx[2][0];
	out.mtx[2][1] = mtx[2][0]*b.mtx[0][1] + mtx[2][1]*b.mtx[1][1] + mtx[2][2]*b.mtx[2][1];
	out.mtx[2][2] = mtx[2][0]*b.mtx[0][2] + mtx[2][1]*b.mtx[1][2] + mtx[2][2]*b.mtx[2][2];

	return out;
}

inline TransMatrix TransMatrix::rotByX(double theta){
	theta = theta*0.017453292;
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);
	double x[3][3] = {{1,0,0}, {0, cosTheta, sinTheta}, {0, -sinTheta, cosTheta}};
	TransMatrix tx = TransMatrix(x);
	return multiply(tx);
}

inline TransMatrix TransMatrix::rotByY(double theta){
	theta = theta*0.017453292;
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);
	double x[3][3] = { {cosTheta,0,-sinTheta}, {0, 1, 0}, {sinTheta, 0, cosTheta}};
	TransMatrix tx = TransMatrix(x);
	return multiply(tx);
}

inline TransMatrix TransMatrix::rotByZ(double theta){
	theta = theta*0.017453292;
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);
	double x[3][3] = { {cosTheta,sinTheta,0}, {-sinTheta, cosTheta, 0}, {0, 0, 1}};
	TransMatrix tx = TransMatrix(x);
	return multiply(tx);
}

inline TransMatrix TransMatrix::randomRot(double theta){
	double rx = 2*theta*(1.0*rand()/RAND_MAX - 0.5);
	double ry = 2*theta*(1.0*rand()/RAND_MAX - 0.5);
	double rz = 2*theta*(1.0*rand()/RAND_MAX - 0.5);
	//cout << rx << " " << ry << " "	<< rz << endl;
	return rotByX(rx).rotByY(ry).rotByZ(rz);
}

} /* namespace NSPmodel */

#endif /* GEOMETRY_TRANSMATRIX_H_ */
