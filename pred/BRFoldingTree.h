/*
 * BRFoldingTree.h
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRFOLDINGTREE_H_
#define PRED_BRFOLDINGTREE_H_

#include "pred/MoveMutator.h"
#include "pred/BRNode.h"
#include "pred/FragmentLibrary.h"
#include "pred/BRConnection.h"
#include "forcefield/RiboConnectHashMap.h"
#include "forcefield/RiboConnectToPO3.h"
#include "model/RNABaseName.h"
#include "model/BaseRotamerLib.h"
#include "forcefield/EnergyTable.h"
#include "forcefield/PO3Builder.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include "pred/BaseMoveLibrary.h"
#include <iostream>
#include <set>

namespace NSPpred {

using namespace std;
using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPgeometry;

class BRTreeInfo {
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BRNode** nodes;
	double ene;
	vector<int> freeNodeIDs;

	BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, double ene);

	BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, vector<BRNode*>& flexibleNodes, double ene);
	//BRTreeInfo(const string& pdbFile);

	double rmsd(BRTreeInfo* other);

	void printPDB(ofstream& of, int modelID);
	void printPDB(const string& outputFile);
	void printTmpPDB(const string& outputFile);
	virtual ~BRTreeInfo();
};



class BRFoldingTree {
public:
	int seqLen;
	int* seq;
	int* wcPairPosID;
	int* nwcPairPosID;

	bool* fixed;
	bool* connectToDownstream;
	bool* loopRevPoints;
	bool* chainBreakPoints;

	bool* nodeConnectMatrix;

	vector<int> freeNodeIDs; //nodes for rmsd calculation
	vector<int> bulgeList;
	vector<int> fixedList;
	vector<int> agList;

	vector<set<int>> fixedGroups;

	BRNode* pseudoNode;


	BaseRotamer** initRotList;

	BRNode** nodes;
	BRNode* rootNode;

	vector<BRConnection*> fixedConnectionList; //move fixed or child fixed
	vector<BRConnection*> flexibleConnectionList; //could be changed during sampling
	vector<BRNode*> flexibleNodes; //base position could be changed during sampling
	vector<BRNode*> riboFlexibleNodes; //ribose rotamer could be changed during sampling


	BRTreeInfo* initTreeInfo;

	FragmentLibrary* fragLib;
	BaseRotamerLib* rotLib;

	int* sepTable;


	double* allBaseClashE;
	double* tmpBaseClashE;
	double* allBaseBaseE;
	double* tmpBaseBaseE;
	double* allBaseRiboseE;
	double* tmpBaseRiboseE;
	double* allBasePhoE;
	double* tmpBasePhoE;
	double* allRiboseRiboseE;
	double* tmpRiboseRiboseE;
	double* allRibosePhoE;
	double* tmpRibosePhoE;
	double* allPhoPhoE;
	double* tmpPhoPhoE;
	double* allRotE;
	double* tmpRotE;
	double* allRcE;
	double* tmpRcE;

	EnergyTable* et;

	BRFoldingTree(const string& inputFile, EnergyTable* et);

	BRFoldingTree(const string& inputFile);
	void randInit();

	void initFromKey(const string& keyInfo);

	void printBaseConnection(){
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<seqLen;j++){
				if(nodeConnectMatrix[i*seqLen+j])
					cout << "1 ";
				else
					cout << "0 ";
			}
			cout << endl;
		}
	}

	void buildFixedNodes();
	void buildGroupNodes();
	void buildAGPairs();
	void buildBasePairOfFixedNodes();
	void buildPairs();
	void buildNeighbor();
	void buildBulged13();


	void buildFrom2(BRNode* node);


	void findMidPoint(){
		for(int i=0;i<seqLen;i++){
			loopRevPoints[i] = false;
			chainBreakPoints[i] = false;
		}
		int pairID[this->seqLen];
		for(int i=0;i<seqLen;i++){
			pairID[i] = this->wcPairPosID[i];
		}

		for(int i=0;i<seqLen;i++){
			if(this->nwcPairPosID[i] >=0) {
				if(pairID[i] >= 0 || pairID[nwcPairPosID[i]] >=0)
					continue;
				pairID[i] = nwcPairPosID[i];
				pairID[nwcPairPosID[i]] = i;
			}
		}


		int loopStart = -1;
		int loopEnd = -1;
		int loopRegion[seqLen];
		for(int i=0;i<seqLen;i++){
			if(fixed[i])
				loopRegion[i] = 2; //fixed region
			else if(wcPairPosID[i] == -1 && nwcPairPosID[i] == -1)
				loopRegion[i] = 0; //loop region
			else if((wcPairPosID[i] > -1 && wcPairPosID[i] > i) || (nwcPairPosID[i] > -1 && nwcPairPosID[i] > i))
				loopRegion[i] = 1; //right bracket
			else
				loopRegion[i] = -1; // left bracket
		}

		for(int i=0;i<seqLen-1;i++){
			if(loopRegion[i] !=0  && loopRegion[i+1] == 0){
				loopStart = i+1;
				loopEnd = i+1;
			}

			if(loopStart > -1 && loopRegion[i] == 0 && loopRegion[i+1] != 0){
				loopEnd = i;
				if((loopEnd-loopStart) >=0 && nodes[i+1]->father != NULL){
					int breakPoint = loopStart + rand()%(loopEnd-loopStart+1);
					chainBreakPoints[breakPoint] = true;
					cout << "chain break: " << breakPoint << endl;
					for(int j=breakPoint;j<=loopEnd;j++){
						loopRevPoints[j] = true;
						cout << "loopRev" << j << endl;
					}
				}

				loopStart = -1;
				loopEnd = -1;
			}
		}
	}

	void buildFrom3(BRNode* node);
	void buildReverseNodes();

	void keyInfoToRMS(const string& keyFile, const string& outFile);

	void updateCoordinate(BRNode* node);
	void updatePhoGroups();
	void updateEnergies(double shift);

	string toCtKey();
	string toCtDetailString();

	pair<double,double> ctMutEnergy(BRConnection* selectConnect, double breakCTWT, double connectWT, double clashWT, double shift, bool verbose);

	vector<double> getF2MutEnergyDetail(BRConnection* selectConnect, double breakCTWT,double connectWT, double clashWT, double shift, bool verbose);

	pair<double,double> f2MutEnergy(BRConnection* selectConnect, double breakCTWT,double connectWT, double clashWT,double shift, bool verbose);

	pair<double,double> f3MutEnergy(BRConnection* ct1,  double breakCTWT,double connectWT, double clashWT,double shift, bool verbose);

	pair<double,double> singleBaseMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT,double shift, bool verbose);

	pair<double,double> baseRotamerMutEnergy(BRNode* node, double breakCTWT,double connectWT, double clashWT,double shift, bool verbose);


	double totalEnergy(double breakCTWT, double connectWT,double clashWT, double shift, bool verbose);

	void printDetailEnergy();

	void printDetailEnergy(ofstream& of);

	void updateCtChildTmpCs(BRConnection* ct, CsMove& cm, bool verbose);

	void clearCtChildTmpCs(BRConnection* ct, bool verbose);

	void acceptCtChildTmpCs(BRConnection* ct, bool verbose);

	void updateF2ChildCs(BRConnection* ct, F2Fragment* frag, bool verbose);

	void updateF2ChildTmpCs(BRConnection* ct, F2Fragment* frag, bool verbose);

	void clearF2ChildTmpCs(BRConnection* ct, bool verbose);

	void acceptF2ChildTmpCs(BRConnection* ct, bool verbose);

	void updateF3ChildTmpCs(BRConnection* ct1, F3Fragment* frag, bool verbose);

	void clearF3ChildTmpCs(BRConnection* ct1, bool verbose);

	void acceptF3ChildTmpCs(BRConnection* ct1,  bool verbose);

	void updateBaseRotamerTmp(BRNode* node, BaseRotamer* rot, bool verbose);

	void clearBaseRotamerTmp(BRNode* node, bool verbose);

	void acceptBaseRotamerTmp(BRNode* node, bool verbose);

	void updateReverseBaseRotamerTmp(BRNode* node, BaseRotamer* rot, bool verbose);
	void clearReverseBaseRotamerTmp(BRNode* node, bool verbose);
	void acceptReverseBaseRotamerTmp(BRNode* node, bool verbose);

	void updateSingleBaseCoordTmp(BRNode* node, CsMove& move, bool verbose);

	void clearSingleBaseCoordTmp(BRNode* node, bool verbose);

	void acceptSingleBaseCoordTmp(BRNode* node, bool verbose);

	void trackCoordinateChangeCt(BRConnection* selectConnect);
	void trackCoordinateChangeF2(BRConnection* selectConnect);
	void trackCoordinateChangeF3(BRConnection* selectConnect);
	void trackCoordinateChangeSingleBase(BRNode* node);
	void trackCoordinateChangeRotamer(BRNode* node);

	void checkConnection();
	void checkNode();
	void checkRibose();
	void checkPho();

	void energyChange();
	void checkTotalEnergy(double shift);
	void checkTmpTotalEnergy(double shift);

	BRTreeInfo* getTreeInfo();

	void printTree(int index);
	int printConnections();

	virtual ~BRFoldingTree();
};

inline double getBaseBaseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID)
		return getBaseBaseEnergy(nodeB, nodeA, abs(sep), et, verbose);
	double bpEnergy = 0.0;
	double minDD, dd;
	int i,j;
	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->baseAtomCoords[0]) < 200.0) {
		minDD = 999999.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->baseAtomCoords[j]);
				if(dd < minDD){
					minDD = dd;
				}
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->cs1, nodeB->cs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}
	if(verbose){
		printf("minDD: %8.3f\n", minDD);
		double len = nodeA->cs1.origin_.distance(nodeB->cs1.origin_);
		//int index = et->bpET->cm2Key.toIndex(nodeA->cs1, nodeB->cs1, len);
		//printf("len: %7.3f index: %6d\n",len, index);

		printf("base %d base %d bpEnergy: %7.3f\n",nodeA->seqID, nodeB->seqID, bpEnergy);
	}
	return bpEnergy;
}

inline double baseBaseClash(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, double shift, bool verbose){
	if(nodeA->seqID > nodeB->seqID)
		return baseBaseClash(nodeB, nodeA, abs(sep), et, shift, verbose);
	double clashEnergy = 0.0;
	double minDD, dd;
	int i,j;
	if(sep == 0) return 0;
	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->baseAtomCoords[0]) < 256.0) {
		clashEnergy = et->atET->getBaseBaseClashEnergy(nodeA->baseType, nodeA->cs1, nodeA->baseAtomCoords, nodeB->baseType, nodeB->cs1, nodeB->baseAtomCoords, shift);
	}

	if(verbose){
		cout << "bb clash energy: " << endl;
		double clashEnergy2 = et->atET->getBaseBaseClashEnergyVerbose(nodeA->baseType, nodeA->cs1, nodeA->baseAtomCoords, nodeB->baseType, nodeB->cs1, nodeB->baseAtomCoords, shift);
		printf("base %d base %d clashEnergy: %7.3f clashEnergy2: %7.3f shift: %5.3f\n",nodeA->seqID, nodeB->seqID, clashEnergy, clashEnergy2, shift);
	}


	return clashEnergy;
}

inline double getBaseBaseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID)
		return getBaseBaseEnergyTmp(nodeB, nodeA, abs(sep), et, verbose);
	double bpEnergy = 0.0;
	double minDD, dd;
	int i,j;


	if(squareDistance(nodeA->baseAtomCoordsTmp[0], nodeB->baseAtomCoordsTmp[0]) < 200.0) {
		minDD = 999999.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoordsTmp[i], nodeB->baseAtomCoordsTmp[j]);
				if(dd < minDD){
					minDD = dd;
				}
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->tmpCs1, nodeB->tmpCs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}
	if(verbose){
	//	printf("base %d base %d bpEnergy: %7.3f clashEnergy: %7.3f\n",nodeA->seqID, nodeB->seqID, bpEnergy, clashEnergy);
	}
	return bpEnergy;
}

inline double baseBaseClashTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, double shift, bool verbose){
	if(nodeA->seqID > nodeB->seqID)
		return baseBaseClashTmp(nodeB, nodeA, abs(sep), et, shift, verbose);
	double clashEnergy = 0.0;
	double minDD, dd;
	int i,j;
	if(sep == 0)
		return 0.0;
	if(squareDistance(nodeA->baseAtomCoordsTmp[0], nodeB->baseAtomCoordsTmp[0]) < 256.0) {
		clashEnergy = et->atET->getBaseBaseClashEnergy(nodeA->baseType, nodeA->tmpCs1, nodeA->baseAtomCoordsTmp, nodeB->baseType, nodeB->tmpCs1, nodeB->baseAtomCoordsTmp, shift);
	}
	if(verbose){
		cout << "bb clash tmp energy: " << endl;
		clashEnergy = et->atET->getBaseBaseClashEnergyVerbose(nodeA->baseType, nodeA->tmpCs1, nodeA->baseAtomCoordsTmp, nodeB->baseType, nodeB->tmpCs1, nodeB->baseAtomCoordsTmp, shift);
		printf("base %d base %d clashEnergy: %7.3f shift: %5.3f\n",nodeA->seqID, nodeB->seqID, clashEnergy, shift);
	}
	return clashEnergy;
}

inline double getBaseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	int i,j;
	double dd;
	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->riboAtomCoords[2]) < 144.0){
		if(abs(sep) > 0) {
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->cs1, nodeB->riboAtomCoords[5]), global2local(nodeA->cs1, nodeB->riboAtomCoords[8]), sep); //O2'
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->cs1, nodeB->riboAtomCoords[6]), global2local(nodeA->cs1, nodeB->riboAtomCoords[9]), sep); //O3'
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->cs1, nodeB->riboAtomCoords[4]), global2local(nodeA->cs1, nodeB->riboAtomCoords[10]), sep); //O4'
			for(i=0;i<nodeA->baseAtomNum;i++){
				for(j=0;j<8;j++){
					dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->riboAtomCoords[j]);
					if(dd < 36)
						clashEnergy += et->atET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(verbose && (baseOxyEnergy+clashEnergy)>100){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
	}
	return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
}



inline double getBaseRiboseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	int i,j;
	double dd;
	if(squareDistance(nodeA->baseAtomCoordsTmp[0], nodeB->riboAtomCoordsTmp[2]) < 144.0){
		if(abs(sep) > 0) {
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 0, global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[5]), global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[8]), sep); //O2'
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 1, global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[6]), global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[9]), sep); //O3'
			baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 2, global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[4]), global2local(nodeA->tmpCs1, nodeB->riboAtomCoordsTmp[10]), sep); //O4'
			for(i=0;i<nodeA->baseAtomNum;i++){
				for(j=0;j<8;j++){
					dd = squareDistance(nodeA->baseAtomCoordsTmp[i], nodeB->riboAtomCoordsTmp[j]);
					if(dd < 36)
						clashEnergy += et->atET->getBaseRiboseEnergy(nodeA->baseType, i, j, dd, sep);
				}
			}
		}
	}
	if(verbose && (baseOxyEnergy+clashEnergy)>100){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
	}
	return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
}

inline double getBasePhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(!nodeB->connectToNeighbor) return 0;
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	int i,j;
	double dd;

	if(sep == 0) return 0.0;

	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 3, global2local(nodeA->cs1, nodeB->pho.tList[1]), global2local(nodeA->cs1, nodeB->pho.tList[4]), sep);
	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 4, global2local(nodeA->cs1, nodeB->pho.tList[2]), global2local(nodeA->cs1, nodeB->pho.tList[5]), sep);
	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 4, global2local(nodeA->cs1, nodeB->pho.tList[3]), global2local(nodeA->cs1, nodeB->pho.tList[6]), sep);

	for(i=0;i<nodeA->baseAtomNum;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->pho.tList[j]);
			if(dd < 36)
				clashEnergy += et->atET->getBasePhoEnergy(nodeA->baseType, i, j, dd, sep);
		}
	}

	if(verbose && (baseOxyEnergy+clashEnergy)>100){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
	}
	return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;

}

inline double getBasePhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(!nodeB->connectToNeighbor) return 0;
	double baseOxyEnergy = 0.0;
	double clashEnergy = 0.0;
	int i,j;
	double dd;
	if(sep == 0) return 0.0;

	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 3, global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[1]), global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[4]), sep);
	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 4, global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[2]), global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[5]), sep);
	baseOxyEnergy += et->roET->getEnergy(nodeA->baseType, 4, global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[3]), global2local(nodeA->tmpCs1, nodeB->phoTmp.tList[6]), sep);

	for(i=0;i<nodeA->baseAtomNum;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->baseAtomCoordsTmp[i], nodeB->phoTmp.tList[j]);
			clashEnergy += et->atET->getBasePhoEnergy(nodeA->baseType, i, j, dd, sep);
		}
	}

	if(verbose && (baseOxyEnergy+clashEnergy)>100){
		cout << "base oxy: " << baseOxyEnergy << endl;
		cout << "clash: " << clashEnergy << endl;
	}
	return baseOxyEnergy*et->para.wtBaseOxygen + clashEnergy;
}

inline double getRiboseRiboseEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID){
		return getRiboseRiboseEnergy(nodeB, nodeA, -sep, et, verbose);
	}
	int i,j;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;
	if(squareDistance(nodeA->riboAtomCoords[2], nodeB->riboAtomCoords[2]) < 81.0) {

		if(abs(sep) >= 2)
			hbondEnergy = et->atET->getRiboseRiboseHbondEnergy(nodeA->riboAtomCoords, nodeB->riboAtomCoords);

		for(i=0;i<8;i++){
			for(j=0;j<8;j++){
				dd = squareDistance(nodeA->riboAtomCoords[i], nodeB->riboAtomCoords[j]);
				if(dd < 36)
					clashEnergy += et->atET->getRiboseRiboseEnergy(i,j,dd, sep);
				if(verbose){
					double d = squareDistance(nodeA->riboAtomCoords[i], nodeB->riboAtomCoords[j]);
					double e = et->atET->getRiboseRiboseEnergy(i,j,squareDistance(nodeA->riboAtomCoords[i], nodeB->riboAtomCoords[j]),sep);
					if(abs(e) > 0.1){
						printf("%d %d %5.3f %8.3f\n",i,j,d,e);
					}
				}
			}
		}
	}
	return clashEnergy + hbondEnergy;
}

inline double getRiboseRiboseEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	int i,j;
	if(nodeA->seqID > nodeB->seqID) {
		return getRiboseRiboseEnergyTmp(nodeB, nodeA, -sep, et, verbose);
	}
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;
	if(squareDistance(nodeA->riboAtomCoordsTmp[2], nodeB->riboAtomCoordsTmp[2]) < 81.0) {

		if(abs(sep) >= 2)
			hbondEnergy = et->atET->getRiboseRiboseHbondEnergy(nodeA->riboAtomCoordsTmp, nodeB->riboAtomCoordsTmp);


		for(i=0;i<8;i++){
			for(j=0;j<8;j++){
				dd = squareDistance(nodeA->riboAtomCoordsTmp[i], nodeB->riboAtomCoordsTmp[j]);
				if(dd < 36)
					clashEnergy += et->atET->getRiboseRiboseEnergy(i,j,dd, sep);
			}
		}
	}
	return clashEnergy + hbondEnergy;
}

inline double getRibosePhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(sep == 0 || sep == -1) return 0;
	if(!nodeB->connectToNeighbor) return 0;
	int i,j;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	if(sep > 2)
		hbondEnergy = et->atET->getRibosePhoHbondEnergy(nodeA->riboAtomCoords, nodeB->pho.tList);

	for(i=0;i<8;i++){
		for(j=1;j<4;j++){

			dd = squareDistance(nodeA->riboAtomCoords[i], nodeB->pho.tList[j]);
			if(dd < 36)
				clashEnergy += et->atET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy + hbondEnergy;
}

inline double getRibosePhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(sep == 0 || sep == -1) return 0;
	if(!nodeB->connectToNeighbor) return 0;
	int i,j;
	double clashEnergy = 0;
	double hbondEnergy = 0;
	double dd;

	if(sep > 2)
		hbondEnergy = et->atET->getRibosePhoHbondEnergy(nodeA->riboAtomCoordsTmp, nodeB->phoTmp.tList);

	for(i=0;i<8;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->riboAtomCoordsTmp[i], nodeB->phoTmp.tList[j]);
			if(dd < 36)
				clashEnergy += et->atET->getRibosePhoEnergy(i,j,dd, sep);
		}
	}
	return clashEnergy + hbondEnergy;
}


inline double getPhoPhoEnergy(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){

	if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

	if(nodeA->seqID > nodeB->seqID){
		return getPhoPhoEnergy(nodeB, nodeA, -sep, et, verbose);
	}
	if(squareDistance(nodeA->pho.tList[0], nodeB->pho.tList[0]) > 36.0) return 0.0;

	int i,j;
	double clashEnergy = 0;
	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->pho.tList[i], nodeB->pho.tList[j]);
			if(dd < 25)
				clashEnergy += et->atET->getPhoPhoEnergy(i,j,dd,sep);
		}
	}
	return clashEnergy;
}

inline double getPhoPhoEnergyTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable* et, bool verbose){
	if(nodeA->seqID > nodeB->seqID) {
		return getPhoPhoEnergyTmp(nodeB, nodeA, -sep, et, verbose);
	}
	if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;
	if(squareDistance(nodeA->phoTmp.tList[0], nodeB->phoTmp.tList[0]) > 36.0) return 0.0;
	int i,j;

	double clashEnergy = 0;
	double dd;
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			dd = squareDistance(nodeA->phoTmp.tList[i], nodeB->phoTmp.tList[j]);
			if(dd < 25)
				clashEnergy += et->atET->getPhoPhoEnergy(i,j,dd,sep);
		}
	}
	return clashEnergy;
}

}


#endif /* PRED_BRFOLDINGTREE_H_ */
