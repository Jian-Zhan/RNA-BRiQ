/*
 * BRNodeBasic.h
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRNODEBASIC_H_
#define PRED_BRNODEBASIC_H_
#include "pred/PhoBasic.h"
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/Angles.h"
#include "model/ProteinRep.h"
#include "model/RnaAtomLib.h"

namespace NSPpred {

class BRConnectionBasic;

using namespace NSPgeometry;
using namespace NSPmodel;

class BRNodeBasic {
public:

	int baseType;
	int seqID;
	LocalFrame cs1;
	LocalFrame tmpCs1;

	bool fixed;
	BRNodeBasic* father;

	BRNodeBasic* leftChild; //pairing base
	BRNodeBasic* midChild; //sequential base
	BRNodeBasic* reverseChild; //previous base

	BRConnectionBasic* upConnection; //pairing move, base move, ribo move
	int baseAtomNum;
	XYZ *atomCoordLocal;
	XYZ *baseAtomCoords;
	XYZ *baseAtomCoordsTmp;

	/*
	int ringNum;
	int *centerTypes;
	XYZ *ringCenterLocal;
	XYZ *ringCenterCoords;
	XYZ *ringCenterCoordsTmp;
	*/

	/*
	PhoBasicLocal phoLocal;
	PhoBasicLocal phoLocalTmp;
	PhoBasic pho;
	PhoBasic phoTmp;
	*/

	vector<bool> childOrNotChild;

	vector<int> childList;
	vector<int> nonChildList;

	vector<int> threeBaseFragListA;
	vector<int> threeBaseFragListB;
	vector<int> threeBaseFragListC;

	/*
	vector<int> phoLocalList;
	vector<int> phoChildList;
	vector<int> phoNonChildList;
	*/

	BRNodeBasic(int baseType, int seqID) {
		this->baseType = baseType;
		this->seqID = seqID;
		this->fixed = false;
		if(baseType == 0){
			baseAtomNum = 7; //center, N1, C2, N3, C8, N7, N6,
			atomCoordLocal = new XYZ[7];
			atomCoordLocal[0] = XYZ(3.291 ,  1.110 ,  0.002); //center
			atomCoordLocal[1] = XYZ(4.276 ,  2.864 ,  0.000); //N1
			atomCoordLocal[2] = XYZ(2.979 ,  3.188 , -0.003); //C2
			atomCoordLocal[3] = XYZ(1.914 ,  2.392 , -0.004); //N3
			atomCoordLocal[4] = XYZ(2.308 , -1.084 ,  0.004); //C8
			atomCoordLocal[5] = XYZ(3.579 , -0.769 ,  0.006); //N7
			atomCoordLocal[6] = XYZ(5.905 ,  1.234 ,  0.008); //N6

			baseAtomCoords = new XYZ[7];
			baseAtomCoordsTmp = new XYZ[7];



		}
		else if(baseType == 1){
			baseAtomNum = 6; //center, O2, N3, O4, C5, C6
			atomCoordLocal = new XYZ[6];
			atomCoordLocal[0] = XYZ(3.015 ,  0.301 ,  0.002); //center
			atomCoordLocal[1] = XYZ(1.534 ,  2.284 , -0.003); //O2
			atomCoordLocal[2] = XYZ(3.495 ,  1.155 ,  0.002); //N3
			atomCoordLocal[3] = XYZ(5.494 ,  0.117 ,  0.006); //O4
			atomCoordLocal[4] = XYZ(3.528 , -1.209 ,  0.005); //C5
			atomCoordLocal[5] = XYZ(2.193 , -1.174 ,  0.003); //C6

			baseAtomCoords = new XYZ[6];
			baseAtomCoordsTmp = new XYZ[6];

		}
		else if(baseType == 2){
			baseAtomNum = 7; //center, N1, N2, N3, C8, N7, O6
			atomCoordLocal = new XYZ[7];
			atomCoordLocal[0] = XYZ(3.218 ,  1.421 , -0.004); //center
			atomCoordLocal[1] = XYZ(4.219 ,  2.835 , -0.005); //N1
			atomCoordLocal[2] = XYZ(2.684 ,  4.547 ,  0.003); //N2
			atomCoordLocal[3] = XYZ(1.884 ,  2.390 ,  0.003); //N3
			atomCoordLocal[4] = XYZ(2.297 , -1.093 , -0.005); //C8
			atomCoordLocal[5] = XYZ(3.563 , -0.778 , -0.009); //N7
			atomCoordLocal[6] = XYZ(5.865 ,  1.269 , -0.013); //O6

			baseAtomCoords = new XYZ[7];
			baseAtomCoordsTmp = new XYZ[7];


		}
		else{
			baseAtomNum = 6; //center, O2, N3, N4, C5, C6
			atomCoordLocal = new XYZ[6];
			atomCoordLocal[0] = XYZ(3.004 ,  0.328 , -0.002); //center
			atomCoordLocal[1] = XYZ(1.496 ,  2.273 ,  0.001); //O2
			atomCoordLocal[2] = XYZ(3.507 ,  1.237 , -0.002); //N3
			atomCoordLocal[3] = XYZ(5.513 ,  0.145 , -0.005); //N4
			atomCoordLocal[4] = XYZ(3.521 , -1.173 , -0.003); //C5
			atomCoordLocal[5] = XYZ(2.183 , -1.171 , -0.001); //C6
			baseAtomCoords = new XYZ[6];
			baseAtomCoordsTmp = new XYZ[6];

		}

		this->father = NULL;
		this->leftChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->upConnection = NULL;


	}

	BRNodeBasic(RNABase* base){
		this->baseType = base->baseTypeInt;
		this->seqID = base->baseSeqID;
		this->fixed = false;
		LocalFrame cs = base->getCoordSystem();
		this->cs1 = cs;
		this->tmpCs1 = cs1;

		if(baseType == 0){ //base A
			baseAtomNum = 7; //center, N1, C2, N3, C8, N7, N6,
			atomCoordLocal = new XYZ[7];
			atomCoordLocal[0] = XYZ(3.291 ,  1.110 ,  0.002);
			atomCoordLocal[1] = XYZ(4.276 ,  2.864 ,  0.000);
			atomCoordLocal[2] = XYZ(2.979 ,  3.188 , -0.003);
			atomCoordLocal[3] = XYZ(1.914 ,  2.392 , -0.004);
			atomCoordLocal[4] = XYZ(2.308 , -1.084 ,  0.004);
			atomCoordLocal[5] = XYZ(3.579 , -0.769 ,  0.006);
			atomCoordLocal[6] = XYZ(5.905 ,  1.234 ,  0.008);
			baseAtomCoords = new XYZ[7];
			baseAtomCoordsTmp = new XYZ[7];
			baseAtomCoords[0] = cs.local2globalcrd(atomCoordLocal[0]);
			baseAtomCoords[1] = base->getAtom("N1")->coord;
			baseAtomCoords[2] = base->getAtom("C2")->coord;
			baseAtomCoords[3] = base->getAtom("N3")->coord;
			baseAtomCoords[4] = base->getAtom("C8")->coord;
			baseAtomCoords[5] = base->getAtom("N7")->coord;
			baseAtomCoords[6] = base->getAtom("N6")->coord;
			for(int i=0;i<7;i++)
				baseAtomCoordsTmp[i] = baseAtomCoords[i];

		}
		else if(baseType == 1){ //base U
			baseAtomNum = 6; //center, O2, N3, O4, C5, C6
			atomCoordLocal = new XYZ[6];
			atomCoordLocal[0] = XYZ(3.015 ,  0.301 ,  0.002);
			atomCoordLocal[1] = XYZ(1.534 ,  2.284 , -0.003);
			atomCoordLocal[2] = XYZ(3.495 ,  1.155 ,  0.002);
			atomCoordLocal[3] = XYZ(5.494 ,  0.117 ,  0.006);
			atomCoordLocal[4] = XYZ(3.528 , -1.209 ,  0.005);
			atomCoordLocal[5] = XYZ(2.193 , -1.174 ,  0.003);
			baseAtomCoords = new XYZ[6];
			baseAtomCoords[0] = cs.local2globalcrd(atomCoordLocal[0]);
			baseAtomCoords[1] = base->getAtom("O2")->coord;
			baseAtomCoords[2] = base->getAtom("N3")->coord;
			baseAtomCoords[3] = base->getAtom("O4")->coord;
			baseAtomCoords[4] = base->getAtom("C5")->coord;
			baseAtomCoords[5] = base->getAtom("C6")->coord;
			baseAtomCoordsTmp = new XYZ[6];
			for(int i=0;i<6;i++)
				baseAtomCoordsTmp[i] = baseAtomCoords[i];
		}
		else if(baseType == 2){ //base G
			baseAtomNum = 7; //center, N1, N2, N3, C8, N7, O6
			atomCoordLocal = new XYZ[7];
			atomCoordLocal[0] = XYZ(3.218 ,  1.421 , -0.004);
			atomCoordLocal[1] = XYZ(4.219 ,  2.835 , -0.005);
			atomCoordLocal[2] = XYZ(2.684 ,  4.547 ,  0.003);
			atomCoordLocal[3] = XYZ(1.884 ,  2.390 ,  0.003);
			atomCoordLocal[4] = XYZ(2.297 , -1.093 , -0.005);
			atomCoordLocal[5] = XYZ(3.563 , -0.778 , -0.009);
			atomCoordLocal[6] = XYZ(5.865 ,  1.269 , -0.013);
			baseAtomCoords = new XYZ[7];
			baseAtomCoords[0] = cs.local2globalcrd(atomCoordLocal[0]);
			baseAtomCoords[1] = base->getAtom("N1")->coord;
			baseAtomCoords[2] = base->getAtom("N2")->coord;
			baseAtomCoords[3] = base->getAtom("N3")->coord;
			baseAtomCoords[4] = base->getAtom("C8")->coord;
			baseAtomCoords[5] = base->getAtom("N7")->coord;
			baseAtomCoords[6] = base->getAtom("O6")->coord;
			baseAtomCoordsTmp = new XYZ[7];
			for(int i=0;i<7;i++)
				baseAtomCoordsTmp[i] = baseAtomCoords[i];


		}
		else{ //base C
			baseAtomNum = 6; //center, O2, N3, N4, C5, C6
			atomCoordLocal = new XYZ[6];
			atomCoordLocal[0] = XYZ(3.004 ,  0.328 , -0.002);
			atomCoordLocal[1] = XYZ(1.496 ,  2.273 ,  0.001);
			atomCoordLocal[2] = XYZ(3.507 ,  1.237 , -0.002);
			atomCoordLocal[3] = XYZ(5.513 ,  0.145 , -0.005);
			atomCoordLocal[4] = XYZ(3.521 , -1.173 , -0.003);
			atomCoordLocal[5] = XYZ(2.183 , -1.171 , -0.001);
			baseAtomCoords = new XYZ[6];
			baseAtomCoords[0] = cs.local2globalcrd(atomCoordLocal[0]);
			baseAtomCoords[1] = base->getAtom("O2")->coord;
			baseAtomCoords[2] = base->getAtom("N3")->coord;
			baseAtomCoords[3] = base->getAtom("N4")->coord;
			baseAtomCoords[4] = base->getAtom("C5")->coord;
			baseAtomCoords[5] = base->getAtom("C6")->coord;
			baseAtomCoordsTmp = new XYZ[6];
			for(int i=0;i<6;i++)
				baseAtomCoordsTmp[i] = baseAtomCoords[i];

		}

		this->father = NULL;
		this->leftChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->upConnection = NULL;

	}

	void updateChildInfo(BRNodeBasic** allNodes, int treeSize) {
		this->childList.clear();
		this->nonChildList.clear();
		this->childOrNotChild.clear();
		for(int i=0;i<treeSize;i++){
			this->childOrNotChild.push_back(false);
		}

		if(leftChild != NULL)
			addNodeFrom(leftChild);
		if(midChild != NULL)
			addNodeFrom(midChild);
		if(reverseChild != NULL)
			addNodeFrom(reverseChild);

		for(int i=0;i<childList.size();i++) {
			childOrNotChild[childList[i]] = true;
		}

		childOrNotChild[seqID] = true;

		for(int i=0;i<treeSize;i++) {
			if(!childOrNotChild[i])
				nonChildList.push_back(i);
		}
		childOrNotChild[seqID] = false;
	}

	void updateThreeBaseChildInfo(BRNodeBasic** allNodes, int treeSize){
		/*
		 * require child info
		 */
		this->threeBaseFragListA.clear();
		this->threeBaseFragListB.clear();
		this->threeBaseFragListC.clear();

		bool x[treeSize];
		for(int i=0;i<treeSize;i++){
			x[i] = false;
		}

		if(midChild != NULL) {
			for(int i=0;i<treeSize;i++){
				if(!midChild->childOrNotChild[i]) {
					x[i] = true;
				}
			}
			x[midChild->seqID] = false;

			for(int i=0;i<treeSize;i++){
				if(x[i])
					threeBaseFragListA.push_back(i);
			}
		}

		if(midChild != NULL && midChild->midChild != NULL){
			for(int i=0;i<treeSize;i++){
				if(midChild->midChild->childOrNotChild[i]) {
					this->threeBaseFragListB.push_back(i);
					x[i] = true;
				}
			}
			this->threeBaseFragListB.push_back(midChild->midChild->seqID);
			x[midChild->midChild->seqID] = true;
		}

		for(int i=0;i<treeSize;i++){
			if(!x[i])
				this->threeBaseFragListC.push_back(i);
		}

	}

	void addNodeFrom(BRNodeBasic* node){

		childList.push_back(node->seqID);
		BRNodeBasic* leftChild = node->leftChild;
		BRNodeBasic* midChild = node->midChild;
		BRNodeBasic* revChild = node->reverseChild;
		BRConnectionBasic* ct;
		if(leftChild != NULL && !leftChild->fixed) {
			addNodeFrom(leftChild);
		}
		if(midChild != NULL && !midChild->fixed) {
			addNodeFrom(midChild);
		}
		if(revChild != NULL && !revChild->fixed) {
			addNodeFrom(revChild);
		}
	}

	BRNodeBasic& operator=(const BRNodeBasic& other);
	void copyValueFrom(const BRNodeBasic& other);
	vector<Atom*> toAtomList(RnaAtomLib& atLib);
	vector<Atom*> toTmpAtomList(RnaAtomLib& atLib);

	void printPartition(){
		cout << "node" << seqID << endl << "listA:";
		for(int i=0;i<threeBaseFragListA.size();i++){
			cout << " " << threeBaseFragListA[i];
		}

		cout << " listB";
		for(int i=0;i<threeBaseFragListB.size();i++){
			cout << " " << threeBaseFragListB[i];
		}

		cout << " listC";
		for(int i=0;i<threeBaseFragListC.size();i++){
			cout << " " << threeBaseFragListC[i];
		}
		cout << endl;

	}

	virtual ~BRNodeBasic();
};

} /* namespace NSPforcefield */

#endif /* PRED_BRNODEBASIC_H_ */
