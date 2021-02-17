/*
 * CsMove.cpp
 *
 *  Created on: Nov 1, 2018
 *      Author: s2982206
 */

#include "geometry/CsMove.h"

namespace NSPgeometry {

CsMove::CsMove() {
}

CsMove::CsMove(const XYZ& oriMove, const TransMatrix& tm){
	this->oriMove = oriMove;
	this->tm = tm;
}

CsMove::CsMove(const string& s){
	vector<string> spt;
	string dim = " ";
	NSPtools::splitString(s, dim, &spt);
	this->oriMove = XYZ(atof(spt[0].c_str()), atof(spt[1].c_str()), atof(spt[2].c_str()));

	XYZ tx(atof(spt[3].c_str()), atof(spt[6].c_str()), atof(spt[9].c_str()));
	XYZ ty(atof(spt[4].c_str()), atof(spt[7].c_str()), atof(spt[10].c_str()));
	XYZ tz(atof(spt[5].c_str()), atof(spt[8].c_str()), atof(spt[11].c_str()));
	tx = ~tx;
	ty = ~(tz^tx);
	tz = tx^ty;


	this->tm.mtx[0][0] = tx.x_;
	this->tm.mtx[0][1] = ty.x_;
	this->tm.mtx[0][2] = tz.x_;
	this->tm.mtx[1][0] = tx.y_;
	this->tm.mtx[1][1] = ty.y_;
	this->tm.mtx[1][2] = tz.y_;
	this->tm.mtx[2][0] = tx.z_;
	this->tm.mtx[2][1] = ty.z_;
	this->tm.mtx[2][2] = tz.z_;
}

CsMove CsMove::reverse(){
	XYZ ori = XYZ(-this->oriMove.x_, -this->oriMove.y_, -this->oriMove.z_);

	TransMatrix revTM = tm.transpose();
	ori = revTM.applyTransition(ori);
	return CsMove(ori, revTM);
}

CsMove CsMove::add(CsMove& other){
	XYZ ori = tm.applyTransition(other.oriMove) + oriMove;
	TransMatrix tm2 =  tm.multiply(other.tm);
	return CsMove(ori, tm2);
}

CsMove CsMove::substract(CsMove& other) {
	CsMove rev = other.reverse();
	return add(rev);
}

string CsMove::toString(){
	char s[200];
	sprintf(s, "%-10.6f %-10.6f %-10.6f %-10.7f %-10.7f %-10.7f %-10.7f %-10.7f %-10.7f %-10.7f %-10.7f %-10.7f",
			oriMove.x_, oriMove.y_, oriMove.z_,
			tm.mtx[0][0], tm.mtx[0][1], tm.mtx[0][2],
			tm.mtx[1][0], tm.mtx[1][1], tm.mtx[1][2],
			tm.mtx[2][0], tm.mtx[2][1], tm.mtx[2][2]);
	return string(s);
}

void CsMove::print(){
	printf("%7.3f %7.3f %7.3f\n", oriMove.x_, oriMove.y_, oriMove.z_);
	tm.print();
}


CsMove::~CsMove() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
