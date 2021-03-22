/*
 * GenerateInitModel.cpp
 *
 *  Created on: Feb 12, 2020
 *      Author: s2982206
 */

#include "pred/BRFoldingTree.h"
#include "pred/BRFoldingTreeBasic.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "forcefield/NeighborConnectionGenerator.h"
#include "pred/MCRun.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;


int main(int argc, char** argv){

	if(argc !=  3 || argv[1] == "-h"){
		cout << "Usage: BRiQ_initPDB $SEQUENCE $OUTPDBFILE" << endl;
		exit(0);
	}
	BaseRotamerLib rotLib;

	string seq = string(argv[1]);
	string outpdb = string(argv[2]);

	CsMove cm("1.179875   3.813016   -3.749098  0.8472664  0.5174403  -0.1199800 -0.5178250 0.8549461  0.0304035  0.1183084  0.0363688  0.9923107");

	vector<int> typeList;
	vector<LocalFrame> csList;

	map<char, int> c2i;
	c2i['A'] = 0;
	c2i['U'] = 1;
	c2i['G'] = 2;
	c2i['C'] = 3;

	LocalFrame cs;

	csList.push_back(cs);
	typeList.push_back(c2i[seq[0]]);

	for(int i=1;i<seq.length();i++){
		cs = cs + cm;

		char c = seq[i];
		if(c2i.find(c) == c2i.end())
			continue;
		else {
			cout << c << endl;
			typeList.push_back(c2i[c]);
			csList.push_back(cs);
		}
	}

	int len = typeList.size();
	bool connectToNeighbor[len];
	int baseTypeSeq[len];
	int id = 0;
	for(int i=0;i<seq.length();i++){
		if(c2i.find(seq[i]) == c2i.end()){
			continue;
		}
		else{
			baseTypeSeq[id] = c2i[seq[i]];
			id++;
		}
	}

	id = 0;
	for(int i=1;i<seq.length();i++){
		if(c2i.find(seq[i]) == c2i.end()){
			connectToNeighbor[id] = false;
			id++;
			i++;
		}
		else {
			connectToNeighbor[id] = true;
			id++;
		}
	}
	connectToNeighbor[len-1] = false;


	BRNode** nodeList = new BRNode*[len];
	for(int i=0;i<len;i++){
		BRNode* node = new BRNode(typeList[i], i);
		node->cs1 = csList[i];
		node->rot = rotLib.rotLib[2][163];
		node->cs2 = node->cs1 + node->rot->mv12;
		for(int j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = local2global(node->cs1, node->atomCoordLocal[j]);
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
		for(int j=0;j<11;j++){
			node->riboAtomCoords[j] = local2global(node->cs1, node->rot->tList1[j]);
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
		nodeList[i] = node;
	}

	Parameter para;
	AtomicEnergyTable* at = new AtomicEnergyTable(&para);
	RiboseOxygenEnergyTable* rET = new RiboseOxygenEnergyTable();
	PO3Builder* pb = new PO3Builder(&para, at, rET);
	for(int i=0;i<len-1;i++){

		BRNode* nodeA = nodeList[i];
		BRNode* nodeB = nodeList[i+1];
		nodeA->connectToNeighbor = connectToNeighbor[i];
		nodeA->phoLocal = pb->getPhoLocal(nodeA->cs2, nodeB->riboAtomCoords, nodeA->rot->improper, nodeB->rot->improper);
		nodeA->pho = PhophateGroup(nodeA->phoLocal, nodeA->cs2);
	}

	BRTreeInfo* info = new BRTreeInfo(len, baseTypeSeq, connectToNeighbor, nodeList, 0.0);
	info->printPDB(outpdb);

}


