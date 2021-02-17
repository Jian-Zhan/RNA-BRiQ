/*
 * RNode.cpp
 *
 *  Created on: Apr 3, 2019
 *      Author: s2982206
 */

#include <pred/RNode.h>

namespace NSPpred {

RNode::RNode() {
	this->seqID = 0;
	this->nodeType = 1;
	this->baseType = 0;
	this->fixed = false;
	this->midChild = NULL;
	this->atomNum = 8;
	this->atomCoordLocal = new XYZ[8];
	this->atomCoords = new XYZ[8];
	this->atomCoordsTmp = new XYZ[8];
	this->pseudoAtomCoords = new XYZ[3];
	this->pseudoAtomCoordsTmp = new XYZ[3];
	this->pseudoAtomLocal = new XYZ[3];
	this->improperAng = 0;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->father = NULL;
	this->upConnection = NULL;
}

RNode::RNode(RNABase* base, int seqID) {
	this->nodeType = 1;
	this->baseType = base->baseTypeInt;
	this->fixed = false;

	Atom* a = base->getAtom("C1'");
	Atom* b = base->getAtom("C2'");
	Atom* c = base->getAtom("C3'");
	Atom* d = base->getAtom("C4'");
	Atom* e = base->getAtom("O4'");
	Atom* f = base->getAtom("O2'");
	Atom* g = base->getAtom("O3'");
	Atom* h = base->getAtom("C5'");
	this->atomNum = 8;
	this->atomCoordLocal = new XYZ[8];
	this->atomCoords = new XYZ[8];
	this->atomCoordsTmp = new XYZ[8];
	this->pseudoAtomCoords = new XYZ[3];
	this->pseudoAtomCoordsTmp = new XYZ[3];
	this->pseudoAtomLocal = new XYZ[3];


	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->midChild = NULL;
	this->upConnection = NULL;

	if(a!=NULL && b!=NULL && c!=NULL && d!=NULL && e!=NULL && f!=NULL && g!=NULL && h!=NULL) {
		this->cs = LocalFrame(d->coord, b->coord, a->coord);
		this->tmpCs = cs;
		LocalFrame cs2(b->coord, c->coord, g->coord);
		LocalFrame cs3(c->coord, d->coord, h->coord);
		move12 = cs2 - cs;
		move31 = cs3 - cs;
		atomCoordLocal[0] = cs.global2localcrd(a->coord);
		atomCoordLocal[1] = cs.global2localcrd(b->coord);
		atomCoordLocal[2] = cs.global2localcrd(c->coord);
		atomCoordLocal[3] = cs.global2localcrd(d->coord);
		atomCoordLocal[4] = cs.global2localcrd(e->coord);
		atomCoordLocal[5] = cs.global2localcrd(f->coord);
		atomCoordLocal[6] = cs.global2localcrd(g->coord);
		atomCoordLocal[7] = cs.global2localcrd(h->coord);

		atomCoords[0] = a->coord;
		atomCoords[1] = b->coord;
		atomCoords[2] = c->coord;
		atomCoords[3] = d->coord;
		atomCoords[4] = e->coord;
		atomCoords[5] = f->coord;
		atomCoords[6] = g->coord;
		atomCoords[7] = h->coord;

		atomCoordsTmp[0] = a->coord;
		atomCoordsTmp[1] = b->coord;
		atomCoordsTmp[2] = c->coord;
		atomCoordsTmp[3] = d->coord;
		atomCoordsTmp[4] = e->coord;
		atomCoordsTmp[5] = f->coord;
		atomCoordsTmp[6] = g->coord;
		atomCoordsTmp[7] = h->coord;

		this->nodeCenterLocal = XYZ(-0.747, -0.744, 0);
		this->nodeCenter = local2global(cs, nodeCenterLocal);
		this->tmpNodeCenter = nodeCenter;

		this->pseudoAtomLocal[0] = XYZ(5.722 ,  4.303 , -1.758);
		this->pseudoAtomLocal[1] = XYZ(-2.635 ,  2.489  , 0.351);
		this->pseudoAtomLocal[2] = XYZ(6.605 , -1.044  , 3.374);

		this->pseudoAtomCoords[0] = local2global(cs, pseudoAtomLocal[0]);
		this->pseudoAtomCoords[1] = local2global(cs, pseudoAtomLocal[1]);
		this->pseudoAtomCoords[2] = local2global(cs, pseudoAtomLocal[2]);

		this->pseudoAtomCoordsTmp[0] = this->pseudoAtomCoords[0];
		this->pseudoAtomCoordsTmp[1] = this->pseudoAtomCoords[1];
		this->pseudoAtomCoordsTmp[2] = this->pseudoAtomCoords[2];

		this->seqID = seqID;
		this->improperAng = dihedral(a->coord, b->coord, d->coord, c->coord);
	}
	else {
		cerr << "backbone atom incomplete: " << base->chainID << " " << base->baseID << endl;
		exit(1);
	}
}

RNode::RNode(XYZ* localTernList, int seqID, int baseType) {
	this->nodeType = 1;
	this->atomNum = 8;
	this->atomCoordLocal = new XYZ[8];
	this->atomCoords = new XYZ[8];
	this->atomCoordsTmp = new XYZ[8];
	this->pseudoAtomCoords = new XYZ[3];
	this->pseudoAtomCoordsTmp = new XYZ[3];
	this->pseudoAtomLocal = new XYZ[3];

	for(int i=0;i<8;i++) {
		this->atomCoordLocal[i] = localTernList[i];
		this->atomCoords[i] = localTernList[i];
		this->atomCoordsTmp[i] = localTernList[i];
	}
	this->father = NULL;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->cs = LocalFrame(localTernList[3], localTernList[1], localTernList[0]);
	this->tmpCs = cs;
	LocalFrame cs2(localTernList[1], localTernList[2], localTernList[6]);
	LocalFrame cs3(localTernList[2], localTernList[3], localTernList[7]);
	move12 = cs.getMove(cs2);
	move31 = cs3.getMove(cs);

	this->nodeCenterLocal = XYZ(-0.747, -0.744, 0);
	this->nodeCenter = local2global(cs, nodeCenterLocal);
	this->tmpNodeCenter = nodeCenter;

	this->pseudoAtomLocal[0] = XYZ(5.722 ,  4.303 , -1.758);
	this->pseudoAtomLocal[1] = XYZ(-2.635 ,  2.489  , 0.351);
	this->pseudoAtomLocal[2] = XYZ(6.605 , -1.044  , 3.374);

	this->pseudoAtomCoords[0] = local2global(cs, pseudoAtomLocal[0]);
	this->pseudoAtomCoords[1] = local2global(cs, pseudoAtomLocal[1]);
	this->pseudoAtomCoords[2] = local2global(cs, pseudoAtomLocal[2]);

	this->pseudoAtomCoordsTmp[0] = this->pseudoAtomCoords[0];
	this->pseudoAtomCoordsTmp[1] = this->pseudoAtomCoords[1];
	this->pseudoAtomCoordsTmp[2] = this->pseudoAtomCoords[2];

	this->seqID = seqID;
	this->baseType = baseType;
	this->fixed = false;
	this->midChild = NULL;
	this->upConnection = NULL;
	this->improperAng = dihedral(localTernList[0], localTernList[1], localTernList[3], localTernList[2]);
}

} /* namespace NSPpred */
