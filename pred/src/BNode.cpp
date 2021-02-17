/*
 * BNode.cpp
 *
 *  Created on: Apr 3, 2019
 *      Author: s2982206
 */

#include <pred/BNode.h>

namespace NSPpred {

BNode::BNode() {
	this->seqID = 0;
	this->baseType = 0;
	this->fixed = false;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->midChild = NULL;
	this->pChild = NULL;
	this->father = NULL;
	this->upConnection = NULL;
	this->nodeType = 0;
	// TODO Auto-generated constructor stub
}

BNode::BNode(LocalFrame& cs,int type, int seqID) {
	this->seqID = seqID;
	this->cs = cs;
	this->tmpCs = cs;
	this->nodeType = 0;

	this->baseType = type;
	if(type == 0){
		nodeCenterLocal = XYZ(3.264, 1.976, 0.0);
	}
	else if(type == 1) {
		nodeCenterLocal = XYZ(2.874, 0.008, 0.0);
	}
	else if(type == 2) {
		nodeCenterLocal = XYZ(3.238, 1.918, 0.0);
	}
	else if(type == 3) {
		nodeCenterLocal = XYZ(2.831, 0.044, 0.0);
	}

	this->nodeCenter = cs.local2globalcrd(nodeCenterLocal);
	this->tmpNodeCenter = nodeCenter;

	this->father = NULL;
	this->leftChild = NULL;
	this->midChild = NULL;
	this->rightChild = NULL;
	this->pChild = NULL;
	this->upConnection = NULL;
	this->fixed = false;
}



} /* namespace NSPpred */
