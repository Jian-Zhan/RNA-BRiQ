/*
 * FoldingTree.h
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#ifndef PRED_FOLDINGTREE_H_
#define PRED_FOLDINGTREE_H_

#include <pred/BaseNode0.h>
#include "pred/BaseConnection.h"
#include "pred/MoveMutator.h"
#include "model/RNABaseName.h"
#include "forcefield/EnergyTable.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include <iostream>
#include <set>

namespace NSPpred {

using namespace std;
using namespace NSPmodel;
using namespace NSPforcefield;

class treeInfo {
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BaseNode** nodes;
	double ene;

	treeInfo(int seqLen, int* seq, bool* con, BaseNode** nodeList, double ene){
		this->seqLen = seqLen;
		this->seq = new int[seqLen];
		this->connectToDownstream = new bool[seqLen];
		for(int i=0;i<seqLen;i++){
			this->seq[i] = seq[i];
			this->connectToDownstream[i] = con[i];
		}
		this->nodes = new BaseNode*[seqLen];
		for(int i=0;i<seqLen;i++) {
			this->nodes[i] = new BaseNode();
			this->nodes[i]->copyValueFrom(nodeList[i]);
		}
		this->ene = ene;
	}

	double rmsd(treeInfo* other);

	void printPDB(const string& outputFile);

	virtual ~treeInfo();

};

class FoldingTree {
public:
	int seqLen;
	int* seq;
	int* wcPairPosID;
	int* nwcPairPosID;

	bool* fixed;
	bool* connectToDownstream;


	vector<set<int>> fixedGroups;

	BaseNode** nodes;
	BaseConnection** connectionToFatherNode;

	vector<BaseConnection*> pseudoConnectList;
	double pseudoConnectEne;


	//vector<BaseNode*> fixedNodes;

	vector<BaseConnection*> moveFixedConnectList;
	vector<BaseConnection*> childFixedConnecList;

	vector<BaseConnection*> flexNbConnectList;
	vector<BaseConnection*> flexWcConnectList;

	int** connectMap;
	EnergyTable* et;

	treeInfo* initTreeInfo;

	FoldingTree(){
		this->seqLen = 0;
		this->seq = NULL;
		this->wcPairPosID = NULL;
		this->nwcPairPosID = NULL;
		this->fixed = NULL;
		this->connectToDownstream = NULL;
		this->nodes = NULL;
		this->connectMap = NULL;
		this->et = NULL;

		this->connectionToFatherNode = NULL;
		this->initTreeInfo = NULL;
		this->pseudoConnectEne = 0.0;
	}

	FoldingTree(const string& baseSeq, const string& ssSeq);
	FoldingTree(const string& inputFile);



	void buildFrom(BaseNode* node);
	void updatePsudoConnection();
	void buildConnectionNearWildType(BaseNode* node, double cutoff);

	void updateConnectMap();
	void updateNeighborConnections() {
		int n = flexNbConnectList.size();
		for(int i=0;i<n;i++){
			flexNbConnectList[i]->updateContactList(connectMap);
		}
	}
	void updateWCConnections() {
		int n = flexWcConnectList.size();
		for(int i=0;i<n;i++){
			flexWcConnectList[i]->updateContactList(connectMap);
		}
	}

	void initFromTreeInfo(treeInfo* ti);
	void randomInit();

	double partialEnergy(BaseConnection* selectConnect);

	void updateChildTmpCs(BaseNode* node);
	void updateChildTmpPho(BaseConnection* selectConnect);
	void updateChildCs(BaseNode* node);
	void updateChildCsAndTmpCs(BaseNode* node);
	void updateChildPho(BaseConnection* selectConnect);

	void updateChildTmpCs(BaseConnection* selectConnect, CsMove* mutMove);

	void updateChildCs(BaseConnection* selectConnect, CsMove* mutMove);

	double partialMutEnergy(BaseConnection* selectConnect, CsMove* mutMove);

	double dPartialMutPseudoConnectEnergy(BaseConnection* selectConnect, CsMove* mutMove);

	void acceptMove(BaseConnection* selectConnect, CsMove* mutMove);
	void accetpMoveWithoutPho(BaseConnection* selectConnect, CsMove* mutMove);
	double getPseudoConnectionEnergy() {
		double e = 0.0;
		for(int i=0;i<pseudoConnectList.size();i++) {
			BaseDistanceMatrix dm(pseudoConnectList[i]->fatherNode->cs, pseudoConnectList[i]->childNode->cs);
			double d = et->nearestDistanceNeighbor(dm, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
			e += d*d;
		}
		return e;
	}

	void printTree(unsigned int index);
	void printTreeInfo() {
		cout << "fixed info:" ;
		for(int i=0;i<seqLen;i++){
			if(fixed[i])
				cout << i << " fixed" << endl;
			else
				cout << i << " flexible" << endl;
		}

		printTree(0);
	}

	void printConnection(unsigned int index=0);

	treeInfo* getTreeInfo();
	void printPDB(const string& outputFile);
	void printTmpPDB(const string& outputFile);
	double totalEnergy();
	void printEnergyDetail();
	void printPho();

	virtual ~FoldingTree();
};

} /* namespace NSPpred */

#endif /* PRED_FOLDINGTREE_H_ */
