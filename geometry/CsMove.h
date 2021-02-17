/*
 * CsMove.h
 *
 *  Created on: Nov 1, 2018
 *      Author: s2982206
 */

#ifndef GEOMETRY_CSMOVE_H_
#define GEOMETRY_CSMOVE_H_

#include "geometry/TransMatrix.h"
#include "geometry/xyz.h"
#include "tools/StringTool.h"
#include <string>
#include <vector>

namespace NSPgeometry {
using namespace std;

class CsMove {
public:
	TransMatrix tm;
	XYZ oriMove;

	CsMove();

	CsMove(const XYZ& oriMove, const TransMatrix& tm);
	CsMove(const string& s);

	CsMove& operator=(const CsMove& other){
		this->tm = other.tm;
		this->oriMove = other.oriMove;
		return *this;
	}

	CsMove reverse();
	CsMove add(CsMove& other);
	CsMove substract(CsMove& other);

	string toString();
	void print();
	void checkCs(){
		tm.checkTM();
	}
	virtual ~CsMove();
};

inline CsMove operator+(const CsMove& cmA, const CsMove& cmB) {
	return CsMove(cmB.oriMove * cmA.tm + cmA.oriMove, cmA.tm * cmB.tm);
}

inline CsMove operator-(const CsMove& cmA, const CsMove& cmB){
	TransMatrix revTM = !cmB.tm;
	return CsMove((-cmB.oriMove)*revTM*cmA.tm + cmA.oriMove, cmA.tm * revTM);
}

inline CsMove operator-(const CsMove& cm) {
	TransMatrix revTM = !cm.tm;
	return CsMove((-cm.oriMove)*revTM, revTM);
}

} /* namespace NSPmodel */

#endif /* GEOMETRY_CSMOVE_H_ */
