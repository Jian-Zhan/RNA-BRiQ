/*
 * BRNode.h
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRNODE_H_
#define PRED_BRNODE_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "geometry/TransMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/Angles.h"
#include "model/ProteinRep.h"
#include "model/RnaAtomLib.h"
#include "model/BaseRotamer.h"
#include "model/PhophateGroup.h"

namespace NSPpred {

class BRConnection;

using namespace NSPgeometry;
using namespace NSPmodel;

class BRNode {
public:

	int baseType;
	int seqID;
	LocalFrame cs1;   //coordinate system of base
	LocalFrame tmpCs1;

	LocalFrame cs2;   //coordinate system C2'-C3'-O3'
	LocalFrame tmpCs2;


	LocalFrame cs3;  //coordinate system O4'-C4'-C5'
	LocalFrame tmpCs3;


	bool fixed;
	bool connectToNeighbor;
	int groupID;

	BaseRotamer* rot;
	BaseRotamer* tmpRot;

	BRNode* father;

	BRNode* leftChild; //pairing base
	BRNode* midChild; //sequential base, connect by base move
	BRNode* rightChild; //jump
	BRNode* bulge13Child; //bulge13
	BRNode* bulge14Child; //bulge14
	BRNode* reverseChild; //previous base
	BRNode* revBulge13Child; //reverse bulge13
	BRNode* revBulge14Child; //reverse bulge14

	BRConnection* upConnection; //pairing move, base move, ribo move

	int baseAtomNum;
	XYZ *atomCoordLocal;
	XYZ *baseAtomCoords;
	XYZ *baseAtomCoordsTmp;

	XYZ riboAtomCoords[11];
	XYZ riboAtomCoordsTmp[11];

	PhophateGroupLocal phoLocal;
	PhophateGroupLocal phoLocalTmp;
	PhophateGroup pho;
	PhophateGroup phoTmp;

	vector<int> chainBreaks;

	vector<int> baseGroupA;
	vector<int> baseGroupC;

	vector<int> riboseGroupA;
	vector<int> riboseGroupC;

	vector<int> phoGroupA;
	vector<int> phoGroupC;

	BRNode(int baseType, int seqID) {
		this->baseType = baseType;
		this->seqID = seqID;
		this->fixed = false;
		this->connectToNeighbor = false;
		this->groupID = -1;
		if(baseType == 0){
			baseAtomNum = 10; //N9, C8, N7, C5, C6, N6, N1, C2, N3, C4
			atomCoordLocal = new XYZ[10];
			atomCoordLocal[0] = XYZ(1.468,   0.000,   0.000); //N9
			atomCoordLocal[1] = XYZ(2.306,  -1.084,   0.000); //C8
			atomCoordLocal[2] = XYZ(3.577,  -0.770,   0.000); //N7
			atomCoordLocal[3] = XYZ(3.576,   0.615,   0.000); //C5
			atomCoordLocal[4] = XYZ(4.614,   1.556,   0.000); //C6
			atomCoordLocal[5] = XYZ(5.904,   1.230,   0.000); //N6
			atomCoordLocal[6] = XYZ(4.276,   2.861,   0.000); //N1
			atomCoordLocal[7] = XYZ(2.980,   3.184,   0.000); //C2
			atomCoordLocal[8] = XYZ(1.914,   2.391,   0.000); //N3
			atomCoordLocal[9] = XYZ(2.285,   1.103,   0.000); //C4

			baseAtomCoords = new XYZ[10];
			baseAtomCoordsTmp = new XYZ[10];

		}
		else if(baseType == 1){
			baseAtomNum = 8; //N1, C2, O2, N3, C4, O4, C5, C6
			atomCoordLocal = new XYZ[8];
			atomCoordLocal[0] = XYZ(1.478,   0.000,   0.000); //N1
			atomCoordLocal[1] = XYZ(2.122,   1.221,   0.000); //C2
			atomCoordLocal[2] = XYZ(1.528,   2.282,   0.000); //O2
			atomCoordLocal[3] = XYZ(3.491,   1.159,   0.000); //N3
			atomCoordLocal[4] = XYZ(4.265,   0.020,   0.000); //C4
			atomCoordLocal[5] = XYZ(5.490,   0.123,   0.000); //O4
			atomCoordLocal[6] = XYZ(3.526,  -1.204,   0.000); //C5
			atomCoordLocal[7] = XYZ(2.191,  -1.173,   0.000); //C6

			baseAtomCoords = new XYZ[8];
			baseAtomCoordsTmp = new XYZ[8];

		}
		else if(baseType == 2){
			baseAtomNum = 11;
			atomCoordLocal = new XYZ[11];
			atomCoordLocal[0] = XYZ(1.468,   0.000,   0.000); //N9
			atomCoordLocal[1] = XYZ(2.295,  -1.094,   0.000); //C8
			atomCoordLocal[2] = XYZ(3.560,  -0.779,   0.000); //N7
			atomCoordLocal[3] = XYZ(3.570,   0.606,   0.000); //C5
			atomCoordLocal[4] = XYZ(4.655,   1.514,   0.000); //C6
			atomCoordLocal[5] = XYZ(5.864,   1.265,   0.000); //O6
			atomCoordLocal[6] = XYZ(4.221,   2.832,   0.000); //N1
			atomCoordLocal[7] = XYZ(2.909,   3.225,   0.000); //C2
			atomCoordLocal[8] = XYZ(2.690,   4.543,   0.000); //N2
			atomCoordLocal[9] = XYZ(1.886,   2.389,   0.000); //N3
			atomCoordLocal[10] = XYZ(2.287,   1.103,   0.000); //C4

			baseAtomCoords = new XYZ[11];
			baseAtomCoordsTmp = new XYZ[11];

		}
		else{
			baseAtomNum = 8;
			atomCoordLocal = new XYZ[8];
			atomCoordLocal[0] = XYZ(1.478,   0.000,   0.000); //N1
			atomCoordLocal[1] = XYZ(2.151,   1.224,   0.000); //C2
			atomCoordLocal[2] = XYZ(1.490,   2.271,   0.000); //O2
			atomCoordLocal[3] = XYZ(3.503,   1.239,   0.000); //N3
			atomCoordLocal[4] = XYZ(4.178,   0.091,   0.000); //C4
			atomCoordLocal[5] = XYZ(5.508,   0.150,   0.000); //N4
			atomCoordLocal[6] = XYZ(3.519,  -1.170,   0.000); //C5
			atomCoordLocal[7] = XYZ(2.181,  -1.170,   0.000); //C6
			baseAtomCoords = new XYZ[8];
			baseAtomCoordsTmp = new XYZ[8];
		}

		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->revBulge13Child = NULL;
		this->revBulge14Child = NULL;
		this->bulge13Child = NULL;
		this->bulge14Child = NULL;
		this->upConnection = NULL;
		this->rot = NULL;
		this->tmpRot = NULL;

	}

	BRNode(RNABase* base, BaseRotamer* rot){
		this->baseType = base->baseTypeInt;
		this->seqID = base->baseSeqID;
		this->fixed = false;
		this->connectToNeighbor = false;
		this->rot = rot;
		this->groupID = -1;
		this->tmpRot = rot;
		LocalFrame cs = base->getCoordSystem();
		this->cs1 = cs;
		this->cs2 = cs1 + rot->mv12;
		this->cs3 = cs1 + rot->mv13;

		this->tmpCs1 = cs1;
		this->tmpCs2 = cs2;
		this->tmpCs3 = cs3;

		if(baseType == 0){ //base A
			baseAtomNum = 10; //N9, C8, N7, C5, C6, N6, N1, C2, N3, C4
			atomCoordLocal = new XYZ[10];
			atomCoordLocal[0] = XYZ(1.468 ,  0.000 , -0.000); //N9
			atomCoordLocal[1] = XYZ(2.307 , -1.084 ,  0.004); //C8
			atomCoordLocal[2] = XYZ(3.578 , -0.769 ,  0.006); //N7
			atomCoordLocal[3] = XYZ(3.576 ,  0.616 ,  0.003); //C5
			atomCoordLocal[4] = XYZ(4.614 ,  1.559 ,  0.004); //C6
			atomCoordLocal[5] = XYZ(5.905 ,  1.233 ,  0.008); //N6
			atomCoordLocal[6] = XYZ(4.275 ,  2.864 , -0.000); //N1
			atomCoordLocal[7] = XYZ(2.978 ,  3.187 , -0.004); //C2
			atomCoordLocal[8] = XYZ(1.912  , 2.392 , -0.004); //N3
			atomCoordLocal[9] = XYZ(2.284 ,  1.103 , -0.000); //C4

			baseAtomCoords = new XYZ[10];
			baseAtomCoordsTmp = new XYZ[10];

			/*
			baseAtomCoords[0] = base->getAtom("N9")->coord;
			baseAtomCoords[1] = base->getAtom("C8")->coord;
			baseAtomCoords[2] = base->getAtom("N7")->coord;
			baseAtomCoords[3] = base->getAtom("C5")->coord;
			baseAtomCoords[4] = base->getAtom("C6")->coord;
			baseAtomCoords[5] = base->getAtom("N6")->coord;
			baseAtomCoords[6] = base->getAtom("N1")->coord;
			baseAtomCoords[7] = base->getAtom("C2")->coord;
			baseAtomCoords[8] = base->getAtom("N3")->coord;
			baseAtomCoords[9] = base->getAtom("C4")->coord;
			*/


			for(int i=0;i<10;i++) {
				baseAtomCoords[i] = local2global(cs1, atomCoordLocal[i]);
				baseAtomCoordsTmp[i] = baseAtomCoords[i];
			}
		}
		else if(baseType == 1){ //base U
			baseAtomNum = 8; //N1, C2, O2, N3, C4, O4, C5, C6
			atomCoordLocal = new XYZ[8];
			atomCoordLocal[0] = XYZ(1.478 , -0.000 ,  0.000); //N1
			atomCoordLocal[1] = XYZ(2.124 ,  1.221 ,  0.000); //C2
			atomCoordLocal[2] = XYZ(1.531 ,  2.284 , -0.003); //O2
			atomCoordLocal[3] = XYZ(3.493 ,  1.156 ,  0.002); //N3
			atomCoordLocal[4] = XYZ(4.267 ,  0.017 ,  0.005); //C4
			atomCoordLocal[5] = XYZ(5.491 ,  0.119 ,  0.007); //O4
			atomCoordLocal[6] = XYZ(3.527 , -1.207 ,  0.006); //C5
			atomCoordLocal[7] = XYZ(2.192 , -1.174 ,  0.004); //C6

			baseAtomCoords = new XYZ[8];
			baseAtomCoordsTmp = new XYZ[8];

			/*
			baseAtomCoords[0] = base->getAtom("N1")->coord;
			baseAtomCoords[1] = base->getAtom("C2")->coord;
			baseAtomCoords[2] = base->getAtom("O2")->coord;
			baseAtomCoords[3] = base->getAtom("N3")->coord;
			baseAtomCoords[4] = base->getAtom("C4")->coord;
			baseAtomCoords[5] = base->getAtom("O4")->coord;
			baseAtomCoords[6] = base->getAtom("C5")->coord;
			baseAtomCoords[7] = base->getAtom("C6")->coord;
			*/

			for(int i=0;i<8;i++) {
				baseAtomCoords[i] = local2global(cs1, atomCoordLocal[i]);
				baseAtomCoordsTmp[i] = baseAtomCoords[i];
			}

		}
		else if(baseType == 2){ //base G
			baseAtomNum = 11;
			atomCoordLocal = new XYZ[11];
			atomCoordLocal[0] = XYZ(1.468 , -0.000 , -0.000); //N9
			atomCoordLocal[1] = XYZ(2.295 , -1.094 , -0.006); //C8
			atomCoordLocal[2] = XYZ(3.561 , -0.780 , -0.010); //N7
			atomCoordLocal[3] = XYZ(3.570 ,  0.606 , -0.007); //C5
			atomCoordLocal[4] = XYZ(4.656 ,  1.514 , -0.010); //C6
			atomCoordLocal[5] = XYZ(5.866 ,  1.265 , -0.016); //O6
			atomCoordLocal[6] = XYZ(4.221 ,  2.833 , -0.006); //N1
			atomCoordLocal[7] = XYZ(2.908 ,  3.226 ,  0.001); //C2
			atomCoordLocal[8] = XYZ(2.688 ,  4.545 ,  0.004); //N2
			atomCoordLocal[9] = XYZ(1.886 ,  2.389 ,  0.004); //N3
			atomCoordLocal[10] = XYZ(2.287 ,  1.103 , -0.000); //C4

			baseAtomCoords = new XYZ[11];
			baseAtomCoordsTmp = new XYZ[11];

			/*
			baseAtomCoords[0] = base->getAtom("N9")->coord;
			baseAtomCoords[1] = base->getAtom("C8")->coord;
			baseAtomCoords[2] = base->getAtom("N7")->coord;
			baseAtomCoords[3] = base->getAtom("C5")->coord;
			baseAtomCoords[4] = base->getAtom("C6")->coord;
			baseAtomCoords[5] = base->getAtom("O6")->coord;
			baseAtomCoords[6] = base->getAtom("N1")->coord;
			baseAtomCoords[7] = base->getAtom("C2")->coord;
			baseAtomCoords[8] = base->getAtom("N2")->coord;
			baseAtomCoords[9] = base->getAtom("N3")->coord;
			baseAtomCoords[10] = base->getAtom("C4")->coord;
			*/

			for(int i=0;i<11;i++) {
				baseAtomCoords[i] = local2global(cs1, atomCoordLocal[i]);
				baseAtomCoordsTmp[i] = baseAtomCoords[i];
			}

		}
		else{ //base C
			baseAtomNum = 8;
			atomCoordLocal = new XYZ[8];
			atomCoordLocal[0] = XYZ(1.479  , 0.000 ,  0.000); //N1
			atomCoordLocal[1] = XYZ(2.152  , 1.225 ,  0.000); //C2
			atomCoordLocal[2] = XYZ(1.492  , 2.273 ,  0.000); //O2
			atomCoordLocal[3] = XYZ(3.504  , 1.238 , -0.001); //N3
			atomCoordLocal[4] = XYZ(4.180  , 0.090 , -0.002); //C4
			atomCoordLocal[5] = XYZ(5.511  , 0.148 , -0.003); //N4
			atomCoordLocal[6] = XYZ(3.520  ,-1.172 , -0.001); //C5
			atomCoordLocal[7] = XYZ(2.182  ,-1.171 , -0.000); //C6
			baseAtomCoords = new XYZ[8];
			baseAtomCoordsTmp = new XYZ[8];

			/*
			baseAtomCoords[0] = base->getAtom("N1")->coord;
			baseAtomCoords[1] = base->getAtom("C2")->coord;
			baseAtomCoords[2] = base->getAtom("O2")->coord;
			baseAtomCoords[3] = base->getAtom("N3")->coord;
			baseAtomCoords[4] = base->getAtom("C4")->coord;
			baseAtomCoords[5] = base->getAtom("N4")->coord;
			baseAtomCoords[6] = base->getAtom("C5")->coord;
			baseAtomCoords[7] = base->getAtom("C6")->coord;
			*/

			for(int i=0;i<8;i++) {
				baseAtomCoords[i] = local2global(cs1, atomCoordLocal[i]);
				baseAtomCoordsTmp[i] = baseAtomCoords[i];
			}
		}

		this->father = NULL;
		this->leftChild = NULL;
		this->rightChild = NULL;
		this->midChild = NULL;
		this->reverseChild = NULL;
		this->upConnection = NULL;
		this->revBulge13Child = NULL;
		this->revBulge14Child = NULL;
		this->bulge13Child = NULL;
		this->bulge14Child = NULL;

		for(int i=0;i<11;i++){
			riboAtomCoords[i] = local2global(cs1, rot->tList1[i]);
			riboAtomCoordsTmp[i] = riboAtomCoords[i];
		}

		this->pho = PhophateGroup(this->cs2);
		this->phoTmp = this->pho;
	}

	void updateRotamer(BaseRotamer* rot){
		this->rot = rot;
		this->tmpRot = rot;
		this->cs2 = cs1 + rot->mv12;
		this->cs3 = cs1 + rot->mv13;

		this->tmpCs2 = cs2;
		this->tmpCs3 = cs3;
		for(int i=0;i<11;i++){
			riboAtomCoords[i] = local2global(cs1, rot->tList1[i]);
			riboAtomCoordsTmp[i] = riboAtomCoords[i];
		}
		this->pho = PhophateGroup(this->cs2);
		this->phoTmp = this->pho;
	}

	bool baseConsistent();
	bool riboConsistent();
	bool rotamerConsistent();

	bool consistent();
	bool phoConsistent();
	bool phoLocalConsistent();

	void updateChildInfo(BRNode** allNodes, int treeSize){

		if(this->seqID > 0 && allNodes[seqID-1]->connectToNeighbor){
			this->chainBreaks.push_back(seqID-1);
		}


		if(connectToNeighbor)
			this->chainBreaks.push_back(seqID);

		for(int i=0;i<treeSize;i++){
			if(i != this->seqID)
				this->baseGroupA.push_back(i);
			else
				this->baseGroupC.push_back(i);
		}

		int* riboseInfo = new int[treeSize];
		int* phoInfo = new int[treeSize];
		for(int i=0;i<treeSize;i++){
			riboseInfo[i] = 0;
			phoInfo[i] = 0;
		}
		riboseInfo[seqID] = 2;

		for(int i=0;i<chainBreaks.size();i++){
			phoInfo[chainBreaks[i]] = 2;
		}
		for(int i=0;i<treeSize;i++){
			if(riboseInfo[i]==0)
				this->riboseGroupA.push_back(i);
			else
				this->riboseGroupC.push_back(i);
		}
		for(int i=0;i<treeSize;i++){
			if(phoInfo[i] == 0 && allNodes[i]->connectToNeighbor)
				this->phoGroupA.push_back(i);
			if(phoInfo[i] == 2 && allNodes[i]->connectToNeighbor)
				this->phoGroupC.push_back(i);
		}
		delete riboseInfo;
		delete phoInfo;
	}

	void printPartition(){
		cout << "single node: " << this->seqID << endl;
		cout << "base groupA:";
		for(int i=0;i<baseGroupA.size();i++){
			cout << " " << baseGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<baseGroupC.size();i++){
			cout << " " << baseGroupC[i];
		}
		cout << endl;

		cout << "ribose groupA:";
		for(int i=0;i<riboseGroupA.size();i++){
			cout << " " << riboseGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<riboseGroupC.size();i++){
			cout << " " << riboseGroupC[i];
		}
		cout << endl;

		cout << "pho groupA:";
		for(int i=0;i<phoGroupA.size();i++){
			cout << " " << phoGroupA[i];
		}
		cout << " groupC:";
		for(int i=0;i<phoGroupC.size();i++){
			cout << " " << phoGroupC[i];
		}
		cout << endl;
	}

	BRNode& operator=(const BRNode& other);
	void copyValueFrom(const BRNode& other);

	vector<Atom*> toAtomList(RnaAtomLib& atLib);
	vector<Atom*> toBaseAtomList(RnaAtomLib& atLib);
	vector<Atom*> phoAtoms();
	vector<Atom*> toTmpAtomList(RnaAtomLib& atLib);

	void checkRotamer();
	void checkTmpRotamer();

	virtual ~BRNode();
};


} /* namespace NSPpred */

#endif /* PRED_BRNODE_H_ */
