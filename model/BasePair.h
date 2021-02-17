/*
 * BasePair.h
 *
 *  Created on: 2020Äê11ÔÂ17ÈÕ
 *      Author: pengx
 */

#ifndef MODEL_BASEPAIR_H_
#define MODEL_BASEPAIR_H_

#include "model/ProteinRep.h"
#include "model/BaseDistanceMatrix.h"
#include "model/RnaAtomLib.h"

namespace NSPmodel {

class BasePair {
public:
	RNABase* baseA;
	RNABase* baseB;
	BaseDistanceMatrix dm;
	string type;
	bool isHbondPair;
	int hbNum;

	BasePair(RNABase* baseA, RNABase* baseB, RnaAtomLib* atLib);
	bool isWCPair();
	bool isHbondedPair();
	double distanceToWCPair();

	virtual ~BasePair();
};

} /* namespace NSPmodel */

#endif /* MODEL_BASEPAIR_H_ */
