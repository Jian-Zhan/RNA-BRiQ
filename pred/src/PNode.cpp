/*
 * PNode.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: s2982206
 */

#include <pred/PNode.h>

namespace NSPpred {

PNode::PNode() {
	// TODO Auto-generated constructor stub
	this->seqID = 0;
	this->baseType = 0;
	this->fixed = false;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->midChild = NULL;
	this->pChild = NULL;
	this->father = NULL;
	this->upConnection = NULL;
	this->nodeType = 2;
}

PNode::PNode(LocalFrame & cs, int type, int seqID) {
	this->seqID = seqID;
	this->cs = cs;
	this->tmpCs = cs;
	this->nodeType = 2;
	this->baseType = type;
	nodeCenterLocal = XYZ(0.0,0.0,0.0);
	nodeCenter = cs.origin_;
	tmpNodeCenter = nodeCenter;
	this->father = NULL;
	this->leftChild = NULL;
	this->midChild = NULL;
	this->rightChild = NULL;
	this->pChild = NULL;
	this->upConnection = NULL;
	this->fixed = false;
}

PNode::~PNode() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpred */
