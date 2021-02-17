/*
 * CsMoveTo6DKey.h
 *
 *  Created on: Sep 6, 2019
 *      Author: s2982206
 */

#ifndef PRED_CSMOVETO6DKEY_H_
#define PRED_CSMOVETO6DKEY_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "model/PhophateGroup.h"
#include "geometry/Angles.h"
#include <fstream>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"


namespace NSPpred {
using namespace std;
using namespace NSPgeometry;

class CsMoveTo6DKey {

private:
	map<int, int> sphereKeyMap;

public:

	CsMoveTo6DKey();
	string toKey(const LocalFrame& csA, const LocalFrame& csB);

	virtual ~CsMoveTo6DKey();
};

} /* namespace NSPforcefield */

#endif /* PRED_CSMOVETO6DKEY_H_ */
