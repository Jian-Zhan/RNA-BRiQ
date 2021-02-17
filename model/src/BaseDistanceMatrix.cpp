/*
 * BaseDistanceMatrix.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: s2982206
 */

#include "model/BaseDistanceMatrix.h"

namespace NSPmodel {

BaseDistanceMatrix::BaseDistanceMatrix() {
	for(int i=0;i<16;i++){
		this->dm[i] = 0.0;
	}
	this->hashKey = "AAAAAAAAAAAAAAAA";
}

BaseDistanceMatrix::BaseDistanceMatrix(const string& line){
	vector<string> list;
	NSPtools::splitString(line, " ", &list);
	double d;
	for(int i=0;i<16;i++){
		d = atof(list[i+2].c_str());
		this->dm[i] = d;
	}
}

BaseDistanceMatrix::BaseDistanceMatrix(LocalFrame& csA, LocalFrame& csB){
	XYZ a(2.158 ,  3.826 ,  1.427);
	XYZ b(-0.789 , -0.329 , -1.273);
	XYZ c(4.520 , -3.006 ,  1.586);
	XYZ d(6.018 ,  1.903 , -1.638);

	XYZ a1 = local2global(csA, a);
	XYZ b1 = local2global(csA, b);
	XYZ c1 = local2global(csA, c);
	XYZ d1 = local2global(csA, d);

	XYZ a2 = local2global(csB, a);
	XYZ b2 = local2global(csB, b);
	XYZ c2 = local2global(csB, c);
	XYZ d2 = local2global(csB, d);

	dm[0] = a1.distance(a2);
	dm[1] = a1.distance(b2);
	dm[2] = a1.distance(c2);
	dm[3] = a1.distance(d2);
	dm[4] = b1.distance(a2);
	dm[5] = b1.distance(b2);
	dm[6] = b1.distance(c2);
	dm[7] = b1.distance(d2);
	dm[8] = c1.distance(a2);
	dm[9] = c1.distance(b2);
	dm[10] = c1.distance(c2);
	dm[11] = c1.distance(d2);
	dm[12] = d1.distance(a2);
	dm[13] = d1.distance(b2);
	dm[14] = d1.distance(c2);
	dm[15] = d1.distance(d2);
}

BaseDistanceMatrix::BaseDistanceMatrix(CsMove& move){
	XYZ a1(2.158 ,  3.826 ,  1.427);
	XYZ b1(-0.789 , -0.329 , -1.273);
	XYZ c1(4.520 , -3.006 ,  1.586);
	XYZ d1(6.018 ,  1.903 , -1.638);

	XYZ a2 = a1 * move.tm + move.oriMove;
	XYZ b2 = b1 * move.tm + move.oriMove;
	XYZ c2 = c1 * move.tm + move.oriMove;
	XYZ d2 = d1 * move.tm + move.oriMove;

	dm[0] = a1.distance(a2);
	dm[1] = a1.distance(b2);
	dm[2] = a1.distance(c2);
	dm[3] = a1.distance(d2);
	dm[4] = b1.distance(a2);
	dm[5] = b1.distance(b2);
	dm[6] = b1.distance(c2);
	dm[7] = b1.distance(d2);
	dm[8] = c1.distance(a2);
	dm[9] = c1.distance(b2);
	dm[10] = c1.distance(c2);
	dm[11] = c1.distance(d2);
	dm[12] = d1.distance(a2);
	dm[13] = d1.distance(b2);
	dm[14] = d1.distance(c2);
	dm[15] = d1.distance(d2);
}

BaseDistanceMatrix::BaseDistanceMatrix(RNABase& baseA, RNABase& baseB){
	baseA.updateCoordSystem();
	baseB.updateCoordSystem();
	if(!baseA.hasLocalFrame || !baseB.hasLocalFrame){
		cerr << "fail to generate distance matrix: " << baseA.baseID << " " << baseB.baseID << endl;
		exit(0);
	}
	LocalFrame csA = baseA.getCoordSystem();
	LocalFrame csB = baseB.getCoordSystem();

	XYZ a(2.158 ,  3.826 ,  1.427);
	XYZ b(-0.789 , -0.329 , -1.273);
	XYZ c(4.520 , -3.006 ,  1.586);
	XYZ d(6.018 ,  1.903 , -1.638);

	XYZ a1 = local2global(csA, a);
	XYZ b1 = local2global(csA, b);
	XYZ c1 = local2global(csA, c);
	XYZ d1 = local2global(csA, d);

	XYZ a2 = local2global(csB, a);
	XYZ b2 = local2global(csB, b);
	XYZ c2 = local2global(csB, c);
	XYZ d2 = local2global(csB, d);

	dm[0] = a1.distance(a2);
	dm[1] = a1.distance(b2);
	dm[2] = a1.distance(c2);
	dm[3] = a1.distance(d2);
	dm[4] = b1.distance(a2);
	dm[5] = b1.distance(b2);
	dm[6] = b1.distance(c2);
	dm[7] = b1.distance(d2);
	dm[8] = c1.distance(a2);
	dm[9] = c1.distance(b2);
	dm[10] = c1.distance(c2);
	dm[11] = c1.distance(d2);
	dm[12] = d1.distance(a2);
	dm[13] = d1.distance(b2);
	dm[14] = d1.distance(c2);
	dm[15] = d1.distance(d2);
}

BaseDistanceMatrix::BaseDistanceMatrix(vector<double>& dList){
	for(int i=0;i<16;i++){
		dm[i] = dList[i];
	}
}

BaseDistanceMatrix::BaseDistanceMatrix(double dm[]){
	for(int i=0;i<16;i++) {
		this->dm[i] = dm[i];
	}
}

double BaseDistanceMatrix::distanceTo(BaseDistanceMatrix& other){
	double tot = 0;
	double d;
	for(int i=0;i<16;i++){
		d = this->dm[i] - other.dm[i];
		tot += d*d;
	}
	return sqrt(tot*0.0625);
}

double BaseDistanceMatrix::squareDistance(BaseDistanceMatrix& other){
	double tot = 0;
	double d;
	for(int i=0;i<9;i++){
		d = this->dm[i] - other.dm[i];
		tot += d*d;
	}
	return tot*0.0625;
}

BaseDistanceMatrix::~BaseDistanceMatrix() {
}

} /* namespace NSPtest */
