/*
 * PNode.h
 *
 *  Created on: Apr 23, 2019
 *      Author: s2982206
 */

#ifndef PRED_PNODE_H_
#define PRED_PNODE_H_

#include "pred/BRNode.h"

namespace NSPpred {

class PNode :public BRNode {
public:
	PNode();
	PNode(LocalFrame& cs, int type, int seqID);
	virtual ~PNode();
};

} /* namespace NSPpred */

#endif /* PRED_PNODE_H_ */
