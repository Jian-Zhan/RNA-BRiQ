/*
 * BRConnectionBasic.h
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRCONNECTIONBASIC_H_
#define PRED_BRCONNECTIONBASIC_H_


#include "pred/BRNodeBasic.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "pred/IndexPair.h"
#include "pred/MoveMutator.h"
#include "pred/BaseMoveLibrary.h"
#include <vector>

namespace NSPpred {

using namespace NSPgeometry;
using namespace std;

class BRConnectionBasic {
public:


	/*
	 * type 0: base-base connect
	 * type 1: reverse connect
	 */
	int connectionType;

	CsMove cm;
	int mvRotType;

	BRNodeBasic* fatherNode;
	BRNodeBasic* childNode;
	int treeSize;

	bool* childOrNotChild;
	vector<int> childList;
	vector<int> nonChildList;

	vector<int> phoLocalList;
	vector<int> phoChildList;
	vector<int> phoNonChildList;

	vector<BRConnectionBasic*> childConnectionList;

	bool fixed;

	BaseMoveLibrary* mvLib;

	BRConnectionBasic();
	BRConnectionBasic(BaseMoveLibrary* mvLib, BRNodeBasic* fatherNode, BRNodeBasic* childNode, int seqLen);

	void updateChildInfo(BRNodeBasic** allNodes) {
		this->childList.clear();
		this->nonChildList.clear();
		this->childConnectionList.clear();
		this->childList.push_back(this->childNode->seqID);

		addNodeFrom(childNode);

		for(int i=0;i<treeSize;i++)
			childOrNotChild[i] = false;
		for(int i=0;i<childList.size();i++) {
			childOrNotChild[childList[i]] = true;
		}
		for(int i=0;i<treeSize;i++) {
			if(!childOrNotChild[i])
				nonChildList.push_back(i);
		}



		for(int i=0;i<treeSize-1;i++) {
			if(childOrNotChild[i] && childOrNotChild[i+1])
				phoChildList.push_back(i);
			else if(!childOrNotChild[i] && childOrNotChild[i+1]){
				phoChildList.push_back(i);
				phoLocalList.push_back(i);
			}
			else if(childOrNotChild[i] && !childOrNotChild[i+1]){
				phoChildList.push_back(i);
				phoLocalList.push_back(i);
			}
			else {
				phoNonChildList.push_back(i);
			}
		}


		if(childOrNotChild[treeSize-1])
			phoChildList.push_back(treeSize-1);

		phoNonChildList.push_back(treeSize-1);
	}

	void addNodeFrom(BRNodeBasic* node){

		BRNodeBasic* leftChild = node->leftChild;
		BRNodeBasic* midChild = node->midChild;
		BRNodeBasic* revChild = node->reverseChild;
		BRConnectionBasic* ct;
		if(leftChild != NULL && !leftChild->fixed) {
			childList.push_back(leftChild->seqID);
			ct = leftChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(leftChild);
		}
		if(midChild != NULL && !midChild->fixed) {
			childList.push_back(midChild->seqID);
			ct = midChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(midChild);
		}
		if(revChild != NULL && !revChild->fixed) {
			childList.push_back(revChild->seqID);
			ct = revChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(revChild);
		}
	}

	void printPartition() {
		cout << "connection: " << this->fatherNode->seqID << "->" << this->childNode->seqID << endl;

		cout << "cNodes:";
		for(int i=0;i<childList.size();i++){
			cout << " " << childList[i];
		}

		cout << " ncNodes:";
		for(int i=0;i<nonChildList.size();i++){
			cout << " " << nonChildList[i];
		}
		cout << endl;

		cout << "localPho: ";
		for(int i=0;i<phoLocalList.size();i++){
			cout << phoLocalList[i] << " ";
		}
		cout << "pho: ";
		for(int i=0;i<phoChildList.size();i++){
			cout << phoChildList[i] << " ";
		}
		cout << endl;
	}

	void printChildConnections() {
		for(int i=0;i<childConnectionList.size();i++) {
			BRConnectionBasic* ct = childConnectionList[i];
			cout << ct->fatherNode->seqID << "-";
			cout << ct->childNode->seqID << " ";
		}
		cout << endl;
	}

	virtual ~BRConnectionBasic();
};

} /* namespace NSPpred */
#endif /* PRED_BRCONNECTIONBASIC_H_ */
