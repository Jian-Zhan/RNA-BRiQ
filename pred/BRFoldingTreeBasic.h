/*
 * BRFoldingTreeBasic.h
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#ifndef PRED_BRFOLDINGTREEBASIC_H_
#define PRED_BRFOLDINGTREEBASIC_H_

#include "pred/MoveMutator.h"
#include "pred/BRNodeBasic.h"
#include "pred/BRConnectionBasic.h"
#include "pred/ThreeBaseMoveLibrary.h"
#include "forcefield/RiboConnectHashMap.h"
#include "forcefield/RiboConnectToPO3.h"
#include "model/RNABaseName.h"
#include "forcefield/EnergyTable.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include "pred/BaseMoveLibrary.h"
#include <iostream>
#include <set>

namespace NSPpred {

using namespace std;
using namespace NSPmodel;
using namespace NSPforcefield;

class BRTreeInfoBasic {
public:
	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BRNodeBasic** nodes;
	double ene;

	BRTreeInfoBasic(int seqLen, int* seq, bool* con, BRNodeBasic** nodes, double ene);

	double rmsd(BRTreeInfoBasic* other);
	void printPDB(const string& outputFile);
	void printTmpPDB(const string& outputFile);
	virtual ~BRTreeInfoBasic();
};

class BRFoldingTreeBasic {
public:
	int seqLen;
	int* seq;
	int* wcPairPosID;
	int* nwcPairPosID;

	bool* fixed;
	bool* connectToDownstream;

	vector<set<int>> fixedGroups;


	BRNodeBasic** nodes;

	vector<BRConnectionBasic*> fixedConnectionList; //move fixed or child fixed

	vector<BRConnectionBasic*> flexibleConnectionList; //could be changed during sampling
	vector<BRNodeBasic*> flexibleNodes; //base rotamer could be changed during sampling

	BRTreeInfoBasic* initTreeInfo;

	BaseMoveLibrary* wcLib;
	BaseMoveLibrary* wcNbLib;
	BaseMoveLibrary* nwcNbLib;
	BaseMoveLibrary* loopNbLib;
	BaseMoveLibrary* riboConnectLib;
	BaseMoveLibrary* revConnectLib;
	BaseMoveLibrary** nwcLibs;

	vector<BRNodeBasic*> threeBaseFragList;
	vector<ThreeBaseMoveLibrary*> threeBaseFragMoveList;



	double* allBaseBaseE;
	double* tmpBaseBaseE;

	EnergyTable* et;

	BRFoldingTreeBasic(const string& inputFile, EnergyTable* et);

	void randInit();
	void buildFrom(BRNodeBasic* node);

	void buildFrom2(BRNodeBasic* node);
	void buildReverseNodes();
	void updateEnergies();


	double connectionMutEnergy(BRConnectionBasic* selectConnect, bool verbose);

	double threeBaseMutEnergy(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose);

	double totalEnergy(bool verbose);

	void printDetailEnergy();

	void updateConnectionChildCs(BRConnectionBasic* ct, CsMove& move, bool verbose);

	void updateConnectionChildTmpCs(BRConnectionBasic* ct, CsMove& move, bool verbose);

	void clearConnectionChildTmpCs(BRConnectionBasic* ct, bool verbose);

	void acceptConnectionChildTmpCs(BRConnectionBasic* ct, CsMove& move, bool verbose);

	void updateThreeBaseFragChildCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose);

	void updateThreeBaseFragChildTmpCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose);

	void clearThreeBaseFragChildTmpCs(BRNodeBasic* node, bool verbose);

	void acceptThreeBaseFragChildTmpCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose);

	void checkConnection();

	void checkCoordinate();
	void checkTotalEnergy();
	void checkTmpTotalEnergy();

	BRTreeInfoBasic* getTreeInfo();

	void printTree(int index);
	void printConnections();
	//double totalEnergy();
	void printNodesPartition(){
		for(int i=0;i<seqLen;i++){
			nodes[i]->printPartition();
		}
	}
	virtual ~BRFoldingTreeBasic();
};


inline string toBaseConnectHashKeyBasic(XYZ* baseA, XYZ* baseB){
	double dm[9];
	dm[0] = squareDistance(baseA[0], baseB[0]);
	dm[1] = squareDistance(baseA[0], baseB[1]);
	dm[2] = squareDistance(baseA[0], baseB[2]);
	dm[3] = squareDistance(baseA[1], baseB[0]);
	dm[4] = squareDistance(baseA[1], baseB[1]);
	dm[5] = squareDistance(baseA[1], baseB[2]);
	dm[6] = squareDistance(baseA[2], baseB[0]);
	dm[7] = squareDistance(baseA[2], baseB[1]);
	dm[8] = squareDistance(baseA[2], baseB[2]);
	char ss[10];
	ss[9] = '\0';
	double d;
	for(int i=0;i<9;i++) {
		d = dm[i];
		if(d<368.64) {
			if(d<92.16) {
				if(d < 23.04) {
					if(d < 5.76) {
						if(d < 1.44)
							ss[i] = 'A';
						else
							ss[i] = 'B';
					}
					else {
						if(d < 12.96)
							ss[i] = 'C';
						else
							ss[i] = 'D';
					}
				}
				else {
					if(d < 51.84) {
						if(d < 36)
							ss[i] = 'E';
						else
							ss[i] = 'F';
					}
					else {
						if(d < 70.56)
							ss[i] = 'G';
						else
							ss[i] = 'H';
					}
				}
			}
			else{
				if(d < 207.36) {
					if(d < 144) {
						if(d < 116.64)
							ss[i] = 'I';
						else
							ss[i] = 'J';
					}
					else {
						if(d < 174.24)
							ss[i] = 'K';
						else
							ss[i] = 'L';
					}
				}
				else {
					if(d < 282.24) {
						if(d < 243.36)
							ss[i] = 'M';
						else
							ss[i] = 'N';
					}
					else {
						if(d < 324)
							ss[i] = 'O';
						else
							ss[i] = 'P';
					}
				}
			}
		}
		else {
			if(d<829.44) {
				if(d < 576) {
					if(d < 466.56) {
						if(d < 416.16)
							ss[i] = 'Q';
						else
							ss[i] = 'R';
					}
					else {
						if(d < 519.84)
							ss[i] = 'S';
						else
							ss[i] = 'T';
					}
				}
				else {
					if(d < 696.96) {
						if(d < 635.04)
							ss[i] = 'U';
						else
							ss[i] = 'V';
					}
					else {
						if(d < 761.76)
							ss[i] = 'W';
						else
							ss[i] = 'X';
					}
				}
			}
			else if(d < 900)
				ss[i] = 'Y';
			else
				ss[i] = 'Z';
		}
	}
	return string(ss);
}

inline double getNonNeighborBaseBaseEnergyBasic(BRNodeBasic* nodeA, BRNodeBasic* nodeB, int sep, EnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return getNonNeighborBaseBaseEnergyBasic(nodeB, nodeA, sep, et, verbose);

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double d1, d2, meanD, minDD, dd;
	int i,j;
	double stackingEnergy = 0.0;

	double wtNearestD, wtPlaneD;
	/*
	 * base base energy
	 */

	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->baseAtomCoords[0]) < 81.0) {
		minDD = 990009.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->baseAtomCoords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
			}
		}

		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->cs1, nodeB->cs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}

	if(verbose && bpEnergy  + clashEnergy + stackingEnergy !=0){
		printf("baseA: %d baseB: %d bpEnergy: %7.3f stackingEnergy: %7.3f clashEnergy: %7.3f\n", nodeA->seqID, nodeB->seqID, bpEnergy ,  stackingEnergy, clashEnergy);
	}

	return bpEnergy  + clashEnergy + stackingEnergy;
}

inline double getNonNeighborBaseBaseEnergyTmpBasic(BRNodeBasic* nodeA, BRNodeBasic* nodeB, int sep, EnergyTable* et, bool verbose){

	if(nodeA->seqID > nodeB->seqID)
		return getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, sep, et, verbose);

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double d1, d2, meanD, minDD, dd;
	int i,j;
	double stackingEnergy = 0.0;

	double wtNearestD, wtPlaneD;
	/*
	 * base base energy
	 */


	if(squareDistance(nodeA->baseAtomCoordsTmp[0], nodeB->baseAtomCoordsTmp[0]) < 81.0) {

		minDD = 990009.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoordsTmp[i], nodeB->baseAtomCoordsTmp[j]);
				if(dd < minDD){
					minDD = dd;
				}
				clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
			}
		}

		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->tmpCs1, nodeB->tmpCs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}

	if(verbose && bpEnergy  + clashEnergy + stackingEnergy !=0){
		printf("baseA: %d baseB: %d bpEnergy: %7.3f stackingEnergy: %7.3f clashEnergy: %7.3f\n", nodeA->seqID, nodeB->seqID, bpEnergy ,  stackingEnergy, clashEnergy);
	}

	return bpEnergy  + clashEnergy + stackingEnergy;
}

inline double getNeighborBaseBaseEnergyBasic(BRNodeBasic* nodeA, BRNodeBasic* nodeB, EnergyTable* et, bool verbose){

	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double connectionEnergy = 0.0;

	int i,j;
	/*
	 * base base energy
	 */
	//cout << "base base energy: " << endl;
	double minDD, dd;
	int sep = abs(nodeA->seqID - nodeB->seqID);
	if(sep > 3) sep = 3;
	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->baseAtomCoords[0]) < 81.0) {
		minDD = 999999.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->baseAtomCoords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->cs1, nodeB->cs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}

	for(i=0;i<nodeA->baseAtomNum;i++){
		for(j=0;j<nodeB->baseAtomNum;j++){
			clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, squareDistance(nodeA->baseAtomCoords[i], nodeB->baseAtomCoords[j]), sep);
		}
	}

	if(verbose && bpEnergy  + clashEnergy + connectionEnergy != 0){
		printf("baseA: %d baseB: %d bpEnergy: %7.3f clashEnergy: %7.3f connectionEnergy: %7.3f\n", nodeA->seqID, nodeB->seqID, bpEnergy, clashEnergy, connectionEnergy);
	}
	return bpEnergy  + clashEnergy + connectionEnergy;
}

inline double getNeighborBaseBaseEnergyTmpBasic(BRNodeBasic* nodeA, BRNodeBasic* nodeB, EnergyTable* et, bool verbose){


	double bpEnergy = 0.0;
	double clashEnergy = 0.0;
	double connectionEnergy = 0.0;

	int i,j;
	/*
	 * base base energy
	 */
	double minDD, dd;
	int sep = abs(nodeA->seqID - nodeB->seqID);
	if(sep > 3) sep = 3;
	if(squareDistance(nodeA->baseAtomCoords[0], nodeB->baseAtomCoords[0]) < 81.0) {
		minDD = 999999.9;
		for(i=0;i<nodeA->baseAtomNum;i++){
			for(j=0;j<nodeB->baseAtomNum;j++){
				dd = squareDistance(nodeA->baseAtomCoords[i], nodeB->baseAtomCoords[j]);
				if(dd < minDD){
					minDD = dd;
				}
				clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, dd, sep);
			}
		}
		if(minDD < 20.25){
			bpEnergy = et->bpET->getEnergy(nodeA->tmpCs1, nodeB->tmpCs1, nodeA->baseType, nodeB->baseType, sep, sqrt(minDD));
		}
	}

	for(i=0;i<nodeA->baseAtomNum;i++){
		for(j=0;j<nodeB->baseAtomNum;j++){
			clashEnergy += et->atET->getBaseBaseEnergy(nodeA->baseType, i, nodeB->baseType, j, squareDistance(nodeA->baseAtomCoordsTmp[i], nodeB->baseAtomCoordsTmp[j]), sep);
		}
	}


	if(verbose && bpEnergy  + clashEnergy != 0){
		printf("baseA: %d baseB: %d bpEnergy: %7.3f clashEnergy: %7.3f\n", nodeA->seqID, nodeB->seqID, bpEnergy, clashEnergy);
	}
	return bpEnergy  + clashEnergy + connectionEnergy;
}

}

#endif /* PRED_BRFOLDINGTREEBASIC_H_ */
