/*
 * RNABaseLib.h
 *
 *  Created on: Nov 18, 2019
 *      Author: s2982206
 */

#ifndef MODEL_RNABASELIB_H_
#define MODEL_RNABASELIB_H_

#include "model/RnaAtomLib.h"
#include "model/ProteinRep.h"
#include "geometry/localframe.h"
#include "geometry/RMSD.h"
#include "model/BaseRotamer.h"
#include "model/PhophateGroup.h"

namespace NSPmodel {

using namespace NSPgeometry;

class RNABaseLib {
public:
	RnaAtomLib atLib;
	RNABaseLib();

	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs);
	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, BaseRotamer* rot);
	RNABase* getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, BaseRotamer* rot, PhophateGroupLocal* pl);

	RNABase* toStandardBase(RNABase* base);

	virtual ~RNABaseLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_RNABASELIB_H_ */
