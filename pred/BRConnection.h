/*
 * BRConnection.h
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRCONNECTION_H_
#define PRED_BRCONNECTION_H_

#include "pred/BRNode.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "pred/IndexPair.h"
#include "pred/MoveMutator.h"
#include "pred/FragmentLibrary.h"
#include <vector>

namespace NSPpred {

using namespace NSPgeometry;
using namespace std;

class BRConnection {
public:

	CsMove cm;
	CsMove cmTmp;

	F2Fragment* f2Frag;
	F2Fragment* f2FragTmp;

	BaseRotamer* childRot;
	BaseRotamer* childRotTmp;

	int mvRotType;
	BRNode* fatherNode;
	BRNode* childNode;
	int treeSize;


	bool* childOrNotChild;
	int* riboseInfo; //0->groupA 1->groupB 2->groupC
	int* phoInfo; //0->groupA 1->groupB 2->groupC

	int* ctPhoInfo;
	string ctType;

	vector<int> ctBaseGroupA;
	vector<int> ctBaseGroupB;
	vector<int> ctRiboseGroupA;
	vector<int> ctRiboseGroupB;
	vector<int> ctPhoGroupA;
	vector<int> ctPhoGroupB;
	vector<int> ctPhoGroupC;

	vector<int> f2BaseGroupA; //coordinate not changed
	vector<int> f2BaseGroupB; //coordinate changed

	vector<int> chainBreaks;

	vector<int> f2RiboseGroupA; //coordinate not changed
	vector<int> f2RiboseGroupB; //coordinate changed, rotamer not changed
	vector<int> f2RiboseGroupC; //rotamer changed

	vector<int> f2PhoGroupA; //coordinate not changed
	vector<int> f2PhoGroupB; //coordinate changed, pho local not changed
	vector<int> f2PhoGroupC; //pho local changed

	vector<BRConnection*> childConnectionList;

	bool hasThreeBaseFragment;
	int tbFragmentIndex;

	int* f3BaseInfo; //0->groupA 1->groupB 2->groupC
	int* f3RiboseInfo; //0->groupA 1->groupB 2->groupC
	int* f3PhoInfo; //0->groupA 1->groupB 2->groupC

	vector<int> f3BaseGroupA; //coordinate not changed
	vector<int> f3BaseGroupB; //coordinate changed, relative coordinate not changed
	vector<int> f3BaseGroupC; //middle base in the f3 fragment

	vector<int> f3ChainBreaks;
	vector<int> f3RiboseGroupA;
	vector<int> f3RiboseGroupB;
	vector<int> f3RiboseGroupC;

	vector<int> f3PhoGroupA;
	vector<int> f3PhoGroupB;
	vector<int> f3PhoGroupC;

	bool fixed;

	FragmentLibrary* fragLib;
	F2FragmentLib* f2Lib;

	F3FragmentLib* f3Lib;

	BRConnection();
	BRConnection(FragmentLibrary* fragLib, const string& connectType, BRNode* fatherNode, BRNode* childNode, int seqLen);

	void updateChildInfo(BRNode** allNodes) {
		if(f2BaseGroupA.size() >0 || f2BaseGroupB.size() > 0) return;

		cout << "update f2 child info: " << this->fatherNode->seqID << " " << this->childNode->seqID << endl;

		this->f2BaseGroupB.push_back(this->childNode->seqID);
		addNodeFrom(childNode);


		for(int i=0;i<treeSize;i++) {
			childOrNotChild[i] = false;
			riboseInfo[i] = -1;
			phoInfo[i] = -1;
		}


		for(int i=0;i<f2BaseGroupB.size();i++) {
			childOrNotChild[f2BaseGroupB[i]] = true;
		}

		f2BaseGroupB.clear();

		for(int i=0;i<treeSize;i++) {
			if(allNodes[i]->fixed)
				childOrNotChild[i] = false;
		}


		for(int i=0;i<treeSize;i++) {
			if(!childOrNotChild[i]) {
				f2BaseGroupA.push_back(i);
				ctBaseGroupA.push_back(i);
				ctRiboseGroupA.push_back(i);
			}
			else {
				f2BaseGroupB.push_back(i);
				ctBaseGroupB.push_back(i);
				ctRiboseGroupB.push_back(i);
			}
		}

		for(int i=0;i<treeSize-1;i++){
			if(allNodes[i]->connectToNeighbor){
				if(childOrNotChild[i] && !childOrNotChild[i+1] && allNodes[i]->connectToNeighbor) this->chainBreaks.push_back(i);
				else if(!childOrNotChild[i] && childOrNotChild[i+1] && allNodes[i]->connectToNeighbor) this->chainBreaks.push_back(i);
			}
		}

		if(f2Lib != NULL && f2Lib->hasRibose) {
			riboseInfo[childNode->seqID] = 2;
		}


		for(int i=0;i<treeSize;i++){
			if(childOrNotChild[i] && riboseInfo[i]<0) riboseInfo[i] = 1;
			if(!childOrNotChild[i] && riboseInfo[i]<0) riboseInfo[i] = 0;
		}



		for(int i=0;i<treeSize;i++){
			phoInfo[i] = riboseInfo[i];
			if(!childOrNotChild[i])
				ctPhoInfo[i] = 0;
			else
				ctPhoInfo[i] = 1;
		}

		for(int i=0;i<chainBreaks.size();i++){
			phoInfo[chainBreaks[i]] = 2;
			ctPhoInfo[chainBreaks[i]] = 2;
		}

		for(int i=0;i<treeSize;i++){
			if(riboseInfo[i] == 2) {
				if(i>0) phoInfo[i-1] = 2;
				phoInfo[i] = 2;
			}
		}


		for(int i=0;i<treeSize;i++){
			if(riboseInfo[i] == 0)
				f2RiboseGroupA.push_back(i);
			else if(riboseInfo[i] == 1)
				f2RiboseGroupB.push_back(i);
			else if(riboseInfo[i] == 2)
				f2RiboseGroupC.push_back(i);

			if(ctPhoInfo[i] == 0 && allNodes[i]->connectToNeighbor)
				ctPhoGroupA.push_back(i);
			else if(ctPhoInfo[i] == 1 && allNodes[i]->connectToNeighbor)
				ctPhoGroupB.push_back(i);
			else if(ctPhoInfo[i] == 2 && allNodes[i]->connectToNeighbor)
				ctPhoGroupC.push_back(i);

			if(phoInfo[i] == 0 && allNodes[i]->connectToNeighbor)
				f2PhoGroupA.push_back(i);
			else if(phoInfo[i] == 1 && allNodes[i]->connectToNeighbor)
				f2PhoGroupB.push_back(i);
			else if(phoInfo[i] == 2 && allNodes[i]->connectToNeighbor)
				f2PhoGroupC.push_back(i);
		}


		if(fatherNode->seqID != childNode->seqID-1)
			this->hasThreeBaseFragment = false;
		else if(fatherNode->groupID >= 0 && fatherNode->groupID == childNode->groupID)
			this->hasThreeBaseFragment = false;
		else if(childNode->groupID >= 0 && childNode->midChild != NULL && childNode->midChild->groupID == childNode->groupID)
			this->hasThreeBaseFragment = false;
		else if(fatherNode->groupID >= 0 && childNode->midChild != NULL && childNode->midChild->groupID == fatherNode->groupID)
			this->hasThreeBaseFragment = false;
		else if(!fatherNode->fixed && !childNode->fixed && childNode->midChild != NULL && !childNode->midChild->fixed && (fatherNode->leftChild == NULL || childNode->leftChild == NULL || childNode->midChild->leftChild == NULL)) {
			this->hasThreeBaseFragment = true;
			this->tbFragmentIndex = fatherNode->baseType*16 + childNode->baseType*4 + childNode->midChild->baseType;
			this->f3Lib = fragLib->acF3Lib[this->tbFragmentIndex];
		}
		else
			this->hasThreeBaseFragment = false;


		if(hasThreeBaseFragment){
			for(int i=0;i<treeSize;i++){
				if(childOrNotChild[i])
					f3BaseInfo[i] = 2;
				else
					f3BaseInfo[i] = 0;
			}

			f3BaseGroupB.push_back(this->childNode->midChild->seqID);
			f3AddNodeFrom(this->childNode->midChild);


			for(int i=0;i<f3BaseGroupB.size();i++){
				f3BaseInfo[f3BaseGroupB[i]] = 1;
			}

			for(int i=0;i<treeSize;i++){
				if(f3BaseInfo[i] == 0)
					f3BaseGroupA.push_back(i);
				else if(f3BaseInfo[i] == 2)
					f3BaseGroupC.push_back(i);
			}

			if(this->fatherNode->seqID > 0 && allNodes[this->fatherNode->seqID-1]->connectToNeighbor)
				this->f3ChainBreaks.push_back(this->fatherNode->seqID-1);

			if(this->childNode->midChild->connectToNeighbor)
				this->f3ChainBreaks.push_back(this->childNode->midChild->seqID);

			for(int i=0;i<treeSize-1;i++){
				if(allNodes[i]->connectToNeighbor){
					if(f3BaseInfo[i] != f3BaseInfo[i+1]) this->f3ChainBreaks.push_back(i);
				}
			}

			for(int i=0;i<treeSize;i++){
				this->f3RiboseInfo[i] = this->f3BaseInfo[i];
				this->f3PhoInfo[i] = this->f3BaseInfo[i];
			}

			for(int i=0;i<f3ChainBreaks.size();i++){
				this->f3PhoInfo[f3ChainBreaks[i]] = 2;
			}

			this->f3RiboseInfo[fatherNode->seqID] = 2;
			this->f3RiboseInfo[childNode->seqID] = 2;
			this->f3RiboseInfo[childNode->midChild->seqID] = 2;


			for(int i=0;i<treeSize;i++){
				if(f3RiboseInfo[i] == 0) this->f3RiboseGroupA.push_back(i);
				else if(f3RiboseInfo[i] == 1) this->f3RiboseGroupB.push_back(i);
				else if(f3RiboseInfo[i] == 2) this->f3RiboseGroupC.push_back(i);
			}

			for(int i=0;i<treeSize;i++){
				if(f3PhoInfo[i] == 0 && allNodes[i]->connectToNeighbor) this->f3PhoGroupA.push_back(i);
				else if(f3PhoInfo[i] == 1 && allNodes[i]->connectToNeighbor) this->f3PhoGroupB.push_back(i);
				else if(f3PhoInfo[i] == 2 && allNodes[i]->connectToNeighbor) this->f3PhoGroupC.push_back(i);
			}

		}
	}

	void f3AddNodeFrom(BRNode* node){
		BRNode* leftChild = node->leftChild;
		BRNode* midChild = node->midChild;
		BRNode* rightChild = node->rightChild;
		BRNode* revChild = node->reverseChild;
		BRNode* bulge13Child = node->bulge13Child;
		BRNode* bulge14Child = node->bulge14Child;
		BRNode* revBulge13Child = node->revBulge13Child;
		BRNode* revBulge14Child = node->revBulge14Child;

		BRConnection* ct;
		if(leftChild != NULL && !leftChild->fixed) {
			f3BaseGroupB.push_back(leftChild->seqID);
			ct = leftChild->upConnection;
			f3AddNodeFrom(leftChild);
		}
		if(midChild != NULL && !midChild->fixed) {
			f3BaseGroupB.push_back(midChild->seqID);
			ct = midChild->upConnection;
			f3AddNodeFrom(midChild);
		}
		if(rightChild != NULL && !rightChild->fixed) {
			f3BaseGroupB.push_back(rightChild->seqID);
			ct = rightChild->upConnection;
			f3AddNodeFrom(rightChild);
		}
		if(revChild != NULL && !revChild->fixed) {
			f3BaseGroupB.push_back(revChild->seqID);
			ct = revChild->upConnection;
			f3AddNodeFrom(revChild);
		}
		if(bulge13Child != NULL && !bulge13Child->fixed){
			f3BaseGroupB.push_back(bulge13Child->seqID);
			ct = bulge13Child->upConnection;
			f3AddNodeFrom(bulge13Child);
		}
		if(bulge14Child != NULL && !bulge14Child->fixed){
			f3BaseGroupB.push_back(bulge14Child->seqID);
			ct = bulge14Child->upConnection;
			f3AddNodeFrom(bulge14Child);
		}
		if(revBulge13Child != NULL && !revBulge13Child->fixed){
			f3BaseGroupB.push_back(revBulge13Child->seqID);
			ct = revBulge13Child->upConnection;
			f3AddNodeFrom(revBulge13Child);
		}
		if(revBulge14Child != NULL && !revBulge14Child->fixed){
			f3BaseGroupB.push_back(revBulge14Child->seqID);
			ct = revBulge14Child->upConnection;
			f3AddNodeFrom(revBulge14Child);
		}
	}

	void addNodeFrom(BRNode* node){
		BRNode* leftChild = node->leftChild;
		BRNode* midChild = node->midChild;
		BRNode* rightChild = node->rightChild;
		BRNode* revChild = node->reverseChild;
		BRNode* bulge13Child = node->bulge13Child;
		BRNode* bulge14Child = node->bulge14Child;
		BRNode* revBulge13Child = node->revBulge13Child;
		BRNode* revBulge14Child = node->revBulge14Child;
		BRConnection* ct;
		if(leftChild != NULL && !leftChild->fixed) {
			f2BaseGroupB.push_back(leftChild->seqID);
			ct = leftChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(leftChild);
		}
		if(midChild != NULL && !midChild->fixed) {
			f2BaseGroupB.push_back(midChild->seqID);
			ct = midChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(midChild);
		}
		if(rightChild != NULL && !rightChild->fixed) {
			f2BaseGroupB.push_back(rightChild->seqID);
			ct = rightChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(rightChild);
		}
		if(revChild != NULL && !revChild->fixed) {
			f2BaseGroupB.push_back(revChild->seqID);
			ct = revChild->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(revChild);
		}
		if(bulge13Child != NULL && !bulge13Child->fixed) {
			f2BaseGroupB.push_back(bulge13Child->seqID);
			ct = bulge13Child->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(bulge13Child);
		}
		if(bulge14Child != NULL && !bulge14Child->fixed) {
			f2BaseGroupB.push_back(bulge14Child->seqID);
			ct = bulge14Child->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(bulge14Child);
		}
		if(revBulge13Child != NULL && !revBulge13Child->fixed) {
			f2BaseGroupB.push_back(revBulge13Child->seqID);
			ct = revBulge13Child->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(revBulge13Child);
		}
		if(revBulge14Child != NULL && !revBulge14Child->fixed) {
			f2BaseGroupB.push_back(revBulge14Child->seqID);
			ct = revBulge14Child->upConnection;
			childConnectionList.push_back(ct);
			addNodeFrom(revBulge14Child);
		}
	}

	void printF3Partition() {
		if(hasThreeBaseFragment){
			cout << "fragment: " << this->fatherNode->seqID << "->" << this->childNode->seqID << "->" << this->childNode->seqID + 1 << endl;
		}
		cout << "base groupA:";
		for(int i=0;i<f3BaseGroupA.size();i++){
			cout << " " << f3BaseGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f3BaseGroupB.size();i++){
			cout << " " << f3BaseGroupB[i];
		}
		cout << " groupC:";
		for(int i=0;i<f3BaseGroupC.size();i++){
			cout << " " << f3BaseGroupC[i];
		}
		cout << endl;

		cout << "ribose groupA:";
		for(int i=0;i<f3RiboseGroupA.size();i++){
			cout << " " << f3RiboseGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f3RiboseGroupB.size();i++){
			cout << " " << f3RiboseGroupB[i];
		}
		cout << " groupC:";
		for(int i=0;i<f3RiboseGroupC.size();i++){
			cout << " " << f3RiboseGroupC[i];
		}
		cout << endl;

		cout << "pho groupA:";
		for(int i=0;i<f3PhoGroupA.size();i++){
			cout << " " << f3PhoGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f3PhoGroupB.size();i++){
			cout << " " << f3PhoGroupB[i];
		}
		cout << " groupC:";
		for(int i=0;i<f3PhoGroupC.size();i++){
			cout << " " << f3PhoGroupC[i];
		}
		cout << endl;
	}

	void printPartition() {
		cout << "connection: " << this->fatherNode->seqID << "->" << this->childNode->seqID << endl;

		cout << "base groupA:";
		for(int i=0;i<f2BaseGroupA.size();i++){
			cout << " " << f2BaseGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f2BaseGroupB.size();i++){
			cout << " " << f2BaseGroupB[i];
		}
		cout << endl;

		cout << "ribose groupA:";
		for(int i=0;i<f2RiboseGroupA.size();i++){
			cout << " " << f2RiboseGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f2RiboseGroupB.size();i++){
			cout << " " << f2RiboseGroupB[i];
		}
		cout << " groupC:";
		for(int i=0;i<f2RiboseGroupC.size();i++){
			cout << " " << f2RiboseGroupC[i];
		}
		cout << endl;

		cout << "pho groupA:";
		for(int i=0;i<f2PhoGroupA.size();i++){
			cout << " " << f2PhoGroupA[i];
		}
		cout << " groupB:";
		for(int i=0;i<f2PhoGroupB.size();i++){
			cout << " " << f2PhoGroupB[i];
		}
		cout << " groupC:";
		for(int i=0;i<f2PhoGroupC.size();i++){
			cout << " " << f2PhoGroupC[i];
		}
		cout << endl;
	}

	void printChildConnections() {
		for(int i=0;i<childConnectionList.size();i++) {
			BRConnection* ct = childConnectionList[i];
			cout << ct->fatherNode->seqID << "-";
			cout << ct->childNode->seqID << " ";
		}
		cout << endl;
	}

	virtual ~BRConnection();
};

} /* namespace NSPpred */

#endif /* PRED_BRCONNECTION_H_ */
