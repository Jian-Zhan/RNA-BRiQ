/*
 * RMSD.h
 *
 *  Created on: Feb 6, 2019
 *      Author: s2982206
 */

#ifndef GEOMETRY_RMSD_H_
#define GEOMETRY_RMSD_H_

#include "geometry/xyz.h"
#include <stdlib.h>
#include <math.h>

namespace NSPgeometry {

using namespace std;

class TransForm{
public:
	double r00, r01, r02, r10, r11, r12, r20, r21, r22, tx, ty, tz;
	TransForm() {
		r00 = 1.0;
		r01 = 0.0;
		r02 = 0.0;
		r10 = 0.0;
		r11 = 1.0;
		r12 = 0.0;
		r20 = 0.0;
		r21 = 0.0;
		r22 = 1.0;
		tx = 0.0;
		ty = 0.0;
		tz = 0.0;
	}

	TransForm(XYZ& q, double w) {
        double X2 = 2.0 * q.x_ * q.x_;
        double Y2 = 2.0 * q.y_ * q.y_;
        double Z2 = 2.0 * q.z_ * q.z_;
        double WX = 2.0 * w * q.x_;
        double WY = 2.0 * w * q.y_;
        double WZ = 2.0 * w * q.z_;
        double XY = 2.0 * q.x_ * q.y_;
        double XZ = 2.0 * q.x_ * q.z_;
        double YZ = 2.0 * q.y_ * q.z_;

        r00 = 1.0 - Y2 - Z2; // Diagonal members
        r11 = 1.0 - X2 - Z2; // "
        r22 = 1.0 - X2 - Y2;  // "
        r01 = XY + WZ;
        r02 = XZ - WY; // Off-diagonal members
        r10 = XY - WZ;
        r12 = YZ + WX; // "
        r20 = XZ + WY;
        r21 = YZ - WX; // "
        tx = 0.0;
        ty = 0.0;
        tz = 0.0;

	}


	TransForm& operator=(const TransForm& other){
		this->r00 = other.r00;
		this->r01 = other.r01;
		this->r02 = other.r02;
		this->r10 = other.r10;
		this->r11 = other.r11;
		this->r12 = other.r12;
		this->r20 = other.r20;
		this->r21 = other.r21;
		this->r22 = other.r22;
		this->tx = other.tx;
		this->ty = other.ty;
		this->tz = other.tz;
		return *this;
	}

	XYZ transform(XYZ& p) {
		double x =  p.x_*r00 + p.y_*r01 + p.z_*r02 + tx;
		double y =  p.x_*r10 + p.y_*r11 + p.z_*r12 + ty;
		double z = 	p.x_*r20 + p.y_*r21 + p.z_*r22 + tz;
		return XYZ(x,y,z);
	}

	void print() {
		printf("%6.3f %6.3f %6.3f\n", r00, r01, r02);
		printf("%6.3f %6.3f %6.3f\n", r10, r11, r12);
		printf("%6.3f %6.3f %6.3f\n", r20, r21, r22);
		printf("%6.3f %6.3f %6.3f\n", tx, ty, tz);

	}

	virtual ~TransForm();

};

void jRotate(double** vecA, int i, int j, int k, int l, double* gh, double s, double tau);
void computeJacobi(double** vecA, double* vecD, double** vecV, int n);

TransForm buildRotation(vector<XYZ>& points1, vector<XYZ>& points2);
double simpleRMSD(const vector<XYZ>& points1, const vector<XYZ>& points2);
XYZ getCOG(const vector<XYZ>& points);
double rmsd(const vector<XYZ>& points1, const vector<XYZ>& points2);
double subRMSD(const vector<XYZ>& points1, const vector<XYZ>& points2, const vector<XYZ>& subList1, const vector<XYZ>& subList2);
double rmsd2(const vector<XYZ>& points1, const vector<XYZ>& points2);




} /* namespace NSPgeometry */

#endif /* GEOMETRY_RMSD_H_ */
