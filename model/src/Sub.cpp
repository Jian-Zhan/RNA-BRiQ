/*
 * Sub.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: s2982206
 */

#include "model/Sub.h"

namespace NSPmodel {

Sub::Sub() {
	this->seqID = 0;
	this->type = 'A';
	this->fatherIndex = -1;
	this->childIndex[0] = -1;
	this->childIndex[1] = -1;
	this->childIndex[2] = -1;
	this->isRoot = false;
	// TODO Auto-generated constructor stub
}

Sub::Sub(LocalFrame& cs, char type, int seqID) {
	this->cs = cs;
	this->type = type;
	this->seqID = seqID;
	this->fatherIndex = -1;
	this->childIndex[0] = -1;
	this->childIndex[1] = -1;
	this->childIndex[2] = -1;
	this->isRoot = false;
}

Sub Sub::applyMove(CsMove& move){
	LocalFrame cs = this->cs.add(move);
	Sub newSub = Sub(cs, type, seqID);
	newSub.tPho = tPho;
	newSub.fatherIndex = fatherIndex;
	newSub.isRoot = isRoot;
	newSub.childIndex[0] = childIndex[0];
	newSub.childIndex[1] = childIndex[1];
	newSub.childIndex[2] = childIndex[2];
	return newSub;
}

Sub Sub::copy(){
	Sub newSub = Sub(cs, type, seqID);
	newSub.tPho = tPho;
	newSub.fatherIndex = fatherIndex;
	newSub.isRoot = isRoot;
	newSub.childIndex[0] = childIndex[0];
	newSub.childIndex[1] = childIndex[1];
	newSub.childIndex[2] = childIndex[2];
	return newSub;
}

Sub::~Sub() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPtest */
