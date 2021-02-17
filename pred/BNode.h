/*
 * BNode.h
 *
 *  Created on: Apr 3, 2019
 *      Author: s2982206
 */

#ifndef PRED_BNODE_H_
#define PRED_BNODE_H_

#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "model/ProteinRep.h"
#include "model/RnaAtomLib.h"
#include "model/RNARibose.h"
#include <vector>
#include "pred/BRNode.h"

namespace NSPpred {
using namespace NSPgeometry;
using namespace std;
using namespace NSPmodel;


class BNode : public BRNode{
public:

	BNode();
	BNode(LocalFrame& cs, int type, int seqID);
};



} /* namespace NSPpred */

#endif /* PRED_BNODE_H_ */
