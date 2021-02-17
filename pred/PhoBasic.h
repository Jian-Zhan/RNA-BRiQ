/*
 * PhoBasic.h
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#ifndef PRED_PHOBASIC_H_
#define PRED_PHOBASIC_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"

namespace NSPpred {
using namespace NSPgeometry;

class PhoBasicLocal {
public:
	XYZ t;
	double e;
	PhoBasicLocal(){
		t = XYZ(  -1.184 , -0.805 , -4.546);
		e = 0.0;
	}

	PhoBasicLocal(const XYZ& t, double e){
		this->t = t;
		this->e = e;
	}

	PhoBasicLocal& operator=(const PhoBasicLocal& other){
		this->t = other.t;
		this->e = other.e;
		return *this;
	}

	virtual ~PhoBasicLocal();
};


class PhoBasic {
public:
	XYZ t;
	double e;
	PhoBasic(){
		this->e = 0.0;
	}

	PhoBasic(const PhoBasicLocal& pl, const LocalFrame& cs){
		this->t = local2global(cs, pl.t);
		this->e = pl.e;
	}

	PhoBasic& operator=(const PhoBasic& other){
		this->t = other.t;
		this->e = other.e;
		return *this;
	}

	virtual ~PhoBasic();
};

} /* namespace NSPpred */

#endif /* PRED_PHOBASIC_H_ */
