/*
 * localframe.h
 *
 */

#ifndef LOCALFRAME_H_
#define LOCALFRAME_H_

#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include <vector>
namespace NSPgeometry {

class LocalFrame {

public:
	XYZ origin_;
	XYZ x;
	TransMatrix localToGlobalTM;
	TransMatrix globalToLocalTM;

	LocalFrame();
	LocalFrame(const XYZ& origin, const XYZ& x, const XYZ& y, const XYZ& z);
	LocalFrame(const XYZ& origin, const TransMatrix& m);
	LocalFrame(const XYZ& A, const XYZ& B, const XYZ& C);
	bool isStandardCS() {
		if(origin_.x_ != 0 || origin_.y_ != 0 || origin_.z_ != 0)
			return false;
		if(localToGlobalTM.mtx[0][0] != 1 || localToGlobalTM.mtx[0][1] != 0 || localToGlobalTM.mtx[0][2] != 0)
			return false;
		if(localToGlobalTM.mtx[1][0] != 0 || localToGlobalTM.mtx[1][1] != 1 || localToGlobalTM.mtx[1][2] != 0)
			return false;
		if(localToGlobalTM.mtx[2][0] != 0 || localToGlobalTM.mtx[2][1] != 0 || localToGlobalTM.mtx[2][2] != 1)
			return false;
		return true;
	}

	LocalFrame& operator=(const LocalFrame& other){
		this->origin_ = other.origin_;
		this->x = other.x;
		this->localToGlobalTM = other.localToGlobalTM;
		this->globalToLocalTM = other.globalToLocalTM;
		return *this;
	}

	CsMove getMove(LocalFrame other){
		XYZ ori = global2localcrd(other.origin_);
		TransMatrix tm = this->globalToLocalTM.multiply(other.localToGlobalTM);
		return CsMove(ori, tm);
	}

	bool equalTo(const LocalFrame& other){
		if(squareDistance(origin_, other.origin_) > 0.0000001)
			return false;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(abs(localToGlobalTM.mtx[i][j] - other.localToGlobalTM.mtx[i][j]) > 0.0001) return false;
				if(abs(globalToLocalTM.mtx[i][j] - other.globalToLocalTM.mtx[i][j]) > 0.0001) return false;
			}
		}
		return true;
	}


	XYZ coordNext(double bondLength, double angle, double torsion);
	LocalFrame csNext(double bondLength, double angle, double torsion);
	XYZ local2globalcrd(XYZ& p);
	XYZ global2localcrd(XYZ&  p);
	LocalFrame add(CsMove& move);

	void print(){
		printf("%7.3f %7.3f %7.3f\n", origin_.x_, origin_.y_, origin_.z_);
		localToGlobalTM.print();
	}
};

inline LocalFrame operator +(const LocalFrame& cs, const CsMove& move){
	return LocalFrame(move.oriMove * cs.localToGlobalTM + cs.origin_, cs.localToGlobalTM * move.tm );
}

inline LocalFrame operator -(const LocalFrame& cs, const CsMove& move){
	XYZ oriMove = -move.oriMove;
	TransMatrix revTM = !move.tm;
	oriMove = oriMove * revTM;
	return LocalFrame(oriMove * cs.localToGlobalTM + cs.origin_, cs.localToGlobalTM * revTM);
}

inline CsMove operator -(const LocalFrame& csA, const LocalFrame& csB) {
	XYZ ori = (csA.origin_ - csB.origin_) * csB.globalToLocalTM;
	TransMatrix tm = csB.globalToLocalTM * csA.localToGlobalTM;
	tm.normalize();
	return CsMove(ori, tm);
}

inline XYZ local2global(const LocalFrame& cs, const XYZ& p) {
	return p * cs.localToGlobalTM + cs.origin_;
}

inline XYZ global2local(const LocalFrame& cs, const XYZ& p) {
	return (p-cs.origin_) * cs.globalToLocalTM;
}

}


#endif /* LOCALFRAME_H_ */
