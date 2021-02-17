/*
 * RNARibose.h
 *
 *  Created on: Apr 3, 2019
 *      Author: s2982206
 */

#ifndef MODEL_RNARIBOSE_H_
#define MODEL_RNARIBOSE_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/ProteinRep.h"

namespace NSPmodel {

using namespace NSPgeometry;

class RNARibose {
public:
	LocalFrame cs;
	LocalFrame tmpCs;

	CsMove move12;
	CsMove move13;
	CsMove move31;
	CsMove move32;

	XYZ* localTernList;
	double improperAng;

	RNARibose();
	RNARibose(RNABase* base);
	RNARibose(const string& line);
	RNARibose& operator=(const RNARibose& other);

	void updateRibose(RNARibose& preRibo, CsMove& riboRibomove);
	void updateRibose(LocalFrame& preBaseCs, CsMove& baseRiboMove);

	virtual ~RNARibose();
};

} /* namespace NSPtools */

#endif /* MODEL_RNARIBOSE_H_ */
