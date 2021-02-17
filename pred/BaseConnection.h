/*
 * BaseConnection.h
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#ifndef PRED_BASECONNECTION_H_
#define PRED_BASECONNECTION_H_
#include <pred/BaseNode0.h>
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "pred/IndexPair.h"
#include "pred/MoveMutator.h"
#include <vector>

namespace NSPpred {

using namespace NSPgeometry;
using namespace std;

class BaseConnection {
public:
	CsMove cm;
	BaseNode* fatherNode;
	BaseNode* childNode;
	int treeSize;
	vector<int> childNodeIndexList;
	bool* childOrNotChild;
	vector<int> nonChildNodeIndexList;

	vector<BaseConnection*> childConnectionList;

	vector<int> updatedPhoIndexList;
	vector<int> unchangedPhoIndexList;

	vector<IndexPair> contactList;
	bool moveFixed;
	bool childFixed;
	MoveMutator* mutator;
	string connectTye;



	BaseConnection(){
		this->mutator = NULL;
		this->fatherNode = NULL;
		this->childNode = NULL;
		this->treeSize = 0;
		this->moveFixed = false;
		this->childFixed = false;
		this->childOrNotChild = NULL;
	}

	BaseConnection(const string& connectType, BaseNode* fatherNode, BaseNode* childNode, int treeSize) {
		this->connectTye = connectType;
		this->fatherNode = fatherNode;
		this->childNode = childNode;
		this->treeSize = treeSize;
		this->mutator = new MoveMutator(connectType, fatherNode->baseType, childNode->baseType);
		this->cm = *(this->mutator->randomMove());
		this->moveFixed = false;
		this->childFixed = false;
		this->childOrNotChild = new bool[treeSize];
		for(int i=0;i<treeSize;i++)
			this->childOrNotChild[i] = false;

	}

	BaseConnection(const string& connectType, CsMove* wtMove, double cutoff, BaseNode* fatherNode, BaseNode* childNode, int treeSize) {
		this->connectTye = connectType;
		this->fatherNode = fatherNode;
		this->childNode = childNode;
		this->treeSize = treeSize;
		this->mutator = new MoveMutator(connectType, wtMove, fatherNode->baseType, childNode->baseType, cutoff);
		this->cm = *(this->mutator->randomMove());
		this->moveFixed = false;
		this->childFixed = false;
		this->childOrNotChild = new bool[treeSize];
		for(int i=0;i<treeSize;i++)
			this->childOrNotChild[i] = false;
	}

	bool isNeighborConnect() {
		return this->fatherNode->seqID == (this->childNode->seqID - 1);
	}

	void setNativeMove() {
		this->cm = fatherNode->cs.getMove(childNode->cs);
	}

	void setMoveFixed(){
		this->moveFixed = true;
		this->cm = fatherNode->cs.getMove(childNode->cs);
	}

	void setChildFixed() {
		this->childFixed = true;
	}

	void updateChildConnectionInfo(BaseConnection** connectionToFatherNode) {
		this->childConnectionList.clear();
		addChildConnectionFrom(this->childNode, connectionToFatherNode);
	}

	void updateChildInfo() {
		this->childNodeIndexList.clear();
		this->nonChildNodeIndexList.clear();
		this->childNodeIndexList.push_back(childNode->seqID);
		addChildNodeFrom(childNode);
		int indexList[treeSize];
		for(int i=0;i<treeSize;i++)
			indexList[i] = 0;
		for(unsigned int i=0;i<childNodeIndexList.size();i++){
			indexList[childNodeIndexList[i]] =1;
			this->childOrNotChild[childNodeIndexList[i]] = true;
		}
		for(int i=0;i<treeSize;i++) {
			if(indexList[i] == 0)
				nonChildNodeIndexList.push_back(i);
		}
		for(int i=treeSize-1;i>0;i--) {
			if(indexList[i-1] == 1)
				indexList[i] = 1;
		}
		for(int i=0;i<treeSize;i++) {
			if(indexList[i] == 0)
				unchangedPhoIndexList.push_back(i);
			else
				updatedPhoIndexList.push_back(i);
		}
	}

	void initRandom() {
		this->cm = *(this->mutator->randomMove());
	}


	void updateContactList(int** connectMap){
		int m = childNodeIndexList.size();
		int n = nonChildNodeIndexList.size();
		int i,j, ni, nj;
		this->contactList.clear();
		for(i=0;i<m;i++){
			ni = childNodeIndexList[i];
			for(j=0;j<n;j++) {
				nj = nonChildNodeIndexList[j];
				if(connectMap[ni][nj] > 0)
					this->contactList.push_back(IndexPair(ni,nj));
			}
		}
	}

	void acceptMove(CsMove* move) {
		this->cm = *move;
	}

	void addChildNodeFrom(BaseNode* node) {
		if(node->leftChild != NULL && !node->leftChild->fixed){
			childNodeIndexList.push_back(node->leftChild->seqID);
			addChildNodeFrom(node->leftChild);
		}
		if(node->rightChild != NULL && !node->rightChild->fixed){
			childNodeIndexList.push_back(node->rightChild->seqID);
			addChildNodeFrom(node->rightChild);
		}
	}

	void addChildConnectionFrom(BaseNode* node, BaseConnection** connectionToFatherNode){
		if(node->leftChild != NULL && !node->leftChild->fixed) {
			childConnectionList.push_back(connectionToFatherNode[node->leftChild->seqID]);
			addChildConnectionFrom(node->leftChild, connectionToFatherNode);
		}
		if(node->rightChild != NULL && !node->rightChild->fixed) {
			childConnectionList.push_back(connectionToFatherNode[node->rightChild->seqID]);
			addChildConnectionFrom(node->rightChild, connectionToFatherNode);
		}
	}

	void printPartition() {
		cout << "partA: ";
		for(int i=0;i<childNodeIndexList.size();i++) {
			cout << childNodeIndexList[i] << " ";
		}
		cout << endl;
		cout << "partB: ";
		for(int i=0;i<nonChildNodeIndexList.size();i++) {
			cout << nonChildNodeIndexList[i] << " ";
		}
		cout << endl;
	}

	void printChildConnection() {
		for(int i=0;i<childConnectionList.size();i++) {
			BaseConnection* bc = childConnectionList[i];
			cout << bc->fatherNode->seqID << "->" << bc->childNode->seqID << " ";
		}
		cout << endl;
	}

	void printContactList() {
		for(int i=0;i<contactList.size();i++) {
			cout << contactList[i].idA_ << " " << contactList[i].idB_ << endl;
		}
	}


	virtual ~BaseConnection();
};

} /* namespace NSPpred */

#endif /* PRED_BASECONNECTION_H_ */
