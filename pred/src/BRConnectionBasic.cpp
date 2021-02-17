/*
 * BRConnectionBasic.cpp
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#include <pred/BRConnectionBasic.h>

namespace NSPpred {

BRConnectionBasic::BRConnectionBasic() {
	// TODO Auto-generated constructor stub
	this->connectionType = 0;
	this->fatherNode = NULL;
	this->childNode = NULL;
	this->treeSize = 0;
	this->childOrNotChild = new bool[1];
	this->fixed = false;
	this->mvLib = NULL;
	this->mvRotType = 0;
}

BRConnectionBasic::BRConnectionBasic(BaseMoveLibrary* mvLib, BRNodeBasic* fatherNode, BRNodeBasic* childNode, int seqLen){
	this->connectionType = connectionType;
	this->fatherNode = fatherNode;
	this->childNode = childNode;
	this->treeSize = seqLen;
	this->childOrNotChild = new bool[treeSize];
	this->fixed = false;
	this->cm = childNode->cs1 - fatherNode->cs1;
	this->mvRotType = 0;
	this->mvLib = mvLib;
	this->connectionType = 0;
	//this->mutator = new MoveMutator(mutatorType, fatherNode->baseType, childNode->baseType);
}

BRConnectionBasic::~BRConnectionBasic() {
	delete [] childOrNotChild;
}

} /* namespace NSPforcefield */
