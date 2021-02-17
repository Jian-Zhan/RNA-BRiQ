/*
 * BaseNode.h
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#ifndef PRED_BASENODE0_H_
#define PRED_BASENODE0_H_

#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "model/ProteinRep.h"
#include "model/RnaAtomLib.h"
#include <vector>

namespace NSPpred {

using namespace NSPgeometry;
using namespace std;
using namespace NSPmodel;

class BaseNode {
public:
	int seqID;
	LocalFrame cs;
	LocalFrame tmpCs;

	XYZ tPho;
	XYZ tmpPho;
	XYZ baseCenter;
	XYZ tmpBaseCenter;
	XYZ baseCenterLocal;
	int baseType;
	bool fixed;

	BaseNode* father;
	BaseNode* leftChild; //3' sequential neighbor
	BaseNode* rightChild; //Watson-Crick pairing partner

	BaseNode();

	BaseNode(LocalFrame& cs, int type, int seqID){
		this->cs = cs;
		this->tmpCs = cs;

		this->baseType = type;
		if(type == 0){
			baseCenterLocal = XYZ(3.287, 1.127, 0.0);
		}
		else if(type == 1) {
			baseCenterLocal = XYZ(3.022, 0.318, 0.0);
		}
		else if(type == 2) {
			baseCenterLocal = XYZ(3.210, 1.444, 0.0);
		}
		else if(type == 3) {
			baseCenterLocal = XYZ(3.012, 0.334, 0.0);
		}
		XYZ localP = XYZ(0.960, -5.091, -0.948);
		this->baseCenter = cs.local2globalcrd(baseCenterLocal);
		this->tmpBaseCenter = baseCenter;
		this->seqID = seqID;
		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->fixed = false;
		this->tPho = cs.local2globalcrd(localP);
		this->tmpPho = this->tPho;
	}

	BaseNode(LocalFrame& cs, XYZ& tp, int type, int seqID){
		this->cs = cs;
		this->tmpCs = cs;

		this->baseType = type;
		if(type == 0){
			baseCenterLocal = XYZ(3.287, 1.127, 0.0);
		}
		else if(type == 1) {
			baseCenterLocal = XYZ(3.022, 0.318, 0.0);
		}
		else if(type == 2) {
			baseCenterLocal = XYZ(3.210, 1.444, 0.0);
		}
		else if(type == 3) {
			baseCenterLocal = XYZ(3.012, 0.334, 0.0);
		}

		this->baseCenter = cs.local2globalcrd(baseCenterLocal);
		this->tmpBaseCenter = baseCenter;
		this->seqID = seqID;
		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->fixed = false;
		this->tPho = tp;
		this->tmpPho = this->tPho;
	}

	BaseNode& operator=(const BaseNode& other) {
		this->seqID = other.seqID;
		this->cs = other.cs;
		this->tmpCs = other.tmpCs;
		this->tPho = other.tPho;
		this->tmpPho = other.tmpPho;
		this->baseCenter = other.baseCenter;
		this->tmpBaseCenter = other.tmpBaseCenter;
		this->baseType = other.baseType;
		this->father = other.father;
		this->leftChild = other.leftChild;
		this->rightChild = other.rightChild;
		this->fixed = other.fixed;

		return *this;
	}

	void copyValueFrom(BaseNode* other) {
		this->seqID = other->seqID;
		this->cs = other->cs;
		this->tmpCs = other->tmpCs;
		this->tPho = other->tPho;
		this->tmpPho = other->tmpPho;
		this->baseCenter = other->baseCenter;
		this->tmpBaseCenter = other->tmpBaseCenter;
		this->baseType = other->baseType;
		this->father = other->father;
		this->leftChild = other->leftChild;
		this->rightChild = other->rightChild;
		this->fixed = other->fixed;
	}

	BaseNode* copy() {
		BaseNode* node = new BaseNode(cs, this->baseType, this->seqID);
		node->father = this->father;
		node->tmpCs = this->tmpCs;
		node->tmpBaseCenter = this->tmpBaseCenter;
		node->tPho = this->tPho;
		node->tmpPho = this->tmpPho;
		node->leftChild = this->leftChild;
		node->rightChild = this->rightChild;
		node->baseCenter = this->baseCenter;
		node->fixed = this->fixed;
		return node;
	}

	void tmpMove(CsMove& cm) {
		this->tmpCs = this->cs.add(cm);
		this->tmpBaseCenter = this->tmpCs.local2globalcrd(this->baseCenterLocal);
	}

	void move(CsMove& cm){
		this->cs = this->cs.add(cm);
		this->baseCenter = this->cs.local2globalcrd(this->baseCenterLocal);
	}

	void updateCs(LocalFrame& cs);

	void updateTmpCs(LocalFrame& cs);

	void updatePho(XYZ& tP){
		//cout << "update pho: " << this->seqID << endl;
		this->tPho = tP;
	}

	void updateTmpPho(XYZ& tP){
		//cout << "update tmpPho: " << this->seqID << endl;
		this->tmpPho = tP;
	}

	vector<Atom*> toAtomList(RnaAtomLib& atLib);
	vector<Atom*> toTmpAtomList(RnaAtomLib& atLib);

	virtual ~BaseNode();
};

inline void BaseNode::updateCs(LocalFrame& cs) {
	this->cs = cs;
	this->baseCenter = this->cs.local2globalcrd(this->baseCenterLocal);
}

inline void BaseNode::updateTmpCs(LocalFrame& cs) {
	this->tmpCs = cs;
	this->tmpBaseCenter = cs.local2globalcrd(this->baseCenterLocal);
}

} /* namespace NSPpred */

#endif /* PRED_BASENODE0_H_ */
