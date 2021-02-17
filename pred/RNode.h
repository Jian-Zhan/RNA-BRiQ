/*
 * RNode.h
 *
 *  Created on: Apr 3, 2019
 *      Author: s2982206
 */

#ifndef PRED_RNODE_H_
#define PRED_RNODE_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/Angles.h"
#include "model/ProteinRep.h"
#include "model/RnaAtomLib.h"
#include "pred/BNode.h"
#include "pred/BRNode.h"
#include <vector>

namespace NSPpred {

using namespace NSPgeometry;

class RNode : public BRNode {
public:

	RNode();
	RNode(RNABase* base, int seqID);
	RNode(XYZ* localTernList, int seqID, int baseType);


};

}

#endif /* PRED_RNODE_H_ */
