/*
 * localframe.cpp
 *
 *  Created on: 2016年3月31日
 *      Author: hyliu
 */

#include "geometry/localframe.h"


namespace NSPgeometry{
	LocalFrame::LocalFrame(){
		x = XYZ(1.0, 0, 0);
	}

	LocalFrame::LocalFrame(const XYZ& origin, const XYZ& x, const XYZ& y, const XYZ& z){
		this->origin_ = origin;
		this->localToGlobalTM = TransMatrix(x,y,z);
		this->x = x;
		this->globalToLocalTM = !localToGlobalTM;
	}

	LocalFrame::LocalFrame(const XYZ& origin, const TransMatrix& m){
		this->origin_ = origin;
		this->localToGlobalTM = m;
		this->x = XYZ(m.mtx[0][0], m.mtx[1][0], m.mtx[2][0]);
		this->globalToLocalTM = !m;
	}

	LocalFrame::LocalFrame(const XYZ& A, const XYZ& B, const XYZ& C){
		XYZ x = ~(C-B);
		XYZ z = ~((A-C)^x);
		XYZ y = z^x;
		this->origin_ = C;
		this->x = x;
		this->localToGlobalTM = TransMatrix(x,y,z);
		this->globalToLocalTM = !this->localToGlobalTM;
	}

	XYZ LocalFrame::coordNext(double bondLength, double angle, double torsion){

		TransMatrix tm2 = localToGlobalTM.rotByX(-torsion).rotByZ(180-angle);
		XYZ inner = XYZ(bondLength, 0, 0);
		return tm2.applyTransition(inner) + this->origin_;
	}

	LocalFrame LocalFrame::csNext(double bondLength, double angle, double torsion){
		LocalFrame lf;
		lf.localToGlobalTM = localToGlobalTM.rotByX(-torsion).rotByZ(180-angle);
		lf.x = XYZ(lf.localToGlobalTM.mtx[0][0], lf.localToGlobalTM.mtx[1][0], lf.localToGlobalTM.mtx[2][0]);
		lf.origin_ = lf.localToGlobalTM.applyTransition(XYZ(bondLength, 0, 0)) + origin_;
		lf.globalToLocalTM = lf.localToGlobalTM.transpose();
		return lf;
	}

	LocalFrame LocalFrame::add(CsMove& move) {
		LocalFrame lf;
		lf.origin_ = this->local2globalcrd(move.oriMove);
		lf.localToGlobalTM = this->localToGlobalTM.multiply(move.tm);
		lf.x = XYZ(lf.localToGlobalTM.mtx[0][0], lf.localToGlobalTM.mtx[1][0], lf.localToGlobalTM.mtx[2][0]);
		lf.globalToLocalTM = lf.localToGlobalTM.transpose();
		return lf;
	}

	XYZ LocalFrame::local2globalcrd(XYZ& p){
		return this->localToGlobalTM.applyTransition(p) + this->origin_;
	}

	XYZ LocalFrame::global2localcrd(XYZ& p) {
		XYZ t = p - this->origin_;
		return this->globalToLocalTM.applyTransition(t);
	}

}

