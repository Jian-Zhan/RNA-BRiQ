/*
 * BackboneModelingTemplate.h
 *
 *  Created on: Sep 30, 2019
 *      Author: s2982206
 */

#ifndef PRED_BACKBONEMODELINGTEMPLATE_H_
#define PRED_BACKBONEMODELINGTEMPLATE_H_

#include "pred/BRNode.h"
#include "pred/BRConnection.h"
#include "pred/BRFoldingTree.h"
#include "forcefield/RiboConnectHashMap.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/RiboConnectToPO3.h"
#include "model/RNABaseName.h"
#include "model/BaseRotamerLib.h"
#include "model/PhophateGroup.h"
#include "forcefield/EnergyTable2.h"
#include "tools/InputParser.h"
#include "geometry/RMSD.h"
#include "geometry/Angles.h"
#include "pred/BaseMoveLibrary.h"
#include <iostream>
#include <set>

namespace NSPpred {

using namespace NSPgeometry;
using namespace std;
using namespace NSPforcefield;
using namespace NSPpara;
using namespace NSPmodel;

class BackboneModelingTemplate {
public:

	int seqLen;
	int* seq;
	bool* connectToDownstream;
	BRNode** nodes;
	EnergyTable2* et;
	PO3Builder* pb;
	BaseRotamerLib* rotLib;
	Parameter* para;

	vector<vector<double>> initDihedsList;
	vector<vector<double>> predDihedsList;

	vector<XYZ> initBackboneAtomList;

	int* sepTable;
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


	BackboneModelingTemplate(const string& inputPDB, const string& paraFile);

	void updateDiheds();

	double calPhoEnergy(double len, double xang3, double xang4, int dihed1, int dihed2, int dihed3, int dihed4, int dihed5, XYZ& p, XYZ& o5, XYZ& op1, XYZ op2, int seqID){
		double e = 0;
		double u;

		double len1 = 1.605;
		double len2 = 1.592;
		double len3 = 1.422;
		double ang1 = 120.1;
		double ang2 = 103.5;
		double ang3 = 120.7;
		double ang4 = 111.1;

		u = (len-len3)*para->kBond;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang3-ang3)*para->kAng;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		u = (xang4-ang4)*para->kAng;
		//e += u*u;

		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += pb->eDihed1[dihed1]*para->wtDihed;
		e += pb->eDihed2[dihed2]*para->wtDihed;
		e += pb->eDihed3[dihed3]*para->wtDihed;
		e += pb->eDihed4[dihed4]*para->wtDihed;
		e += pb->eDihed5[dihed5]*para->wtDihed;

		double baseOxyEnergy = 0;
		double clashEnergy = 0;
		/*

		BRNode* node;
		for(int i=0;i<seqLen;i++){
			if(i == seqID || i == seqID + 1) continue;
			node = nodes[i];
			if(squareDistance(node->cs1.origin_, p) < 144){
				// base pho energy
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 3, global2local(node->cs1, o5));
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 4, global2local(node->cs1, op1));
				baseOxyEnergy += et->roET.getEnergy(node->baseType, 4, global2local(node->cs1, op2));
				for(int j=0;j<node->baseAtomNum;j++) {
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 1, squareDistance(node->baseAtomCoords[j], o5));
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 2, squareDistance(node->baseAtomCoords[j], op1));
					clashEnergy += et->atET.getBasePhoEnergy(node->baseType, j, 3, squareDistance(node->baseAtomCoords[j], op2));
				}
				// ribose pho energy

				for(int j=0;j<8;j++){
					clashEnergy += et->atET.getRibosePhoEnergy(j,1,squareDistance(node->riboAtomCoords[j], o5));
					clashEnergy += et->atET.getRibosePhoEnergy(j,2,squareDistance(node->riboAtomCoords[j], op1));
					clashEnergy += et->atET.getRibosePhoEnergy(j,3,squareDistance(node->riboAtomCoords[j], op2));
				}

				// pho pho energy
				for(int j=1;j<4;j++){
					clashEnergy += et->atET.getPhoPhoEnergy(j,1,squareDistance(node->pho.tList[j], o5));
					clashEnergy += et->atET.getPhoPhoEnergy(j,2,squareDistance(node->pho.tList[j], op1));
					clashEnergy += et->atET.getPhoPhoEnergy(j,3,squareDistance(node->pho.tList[j], op2));
				}
			}
		}
		*/

		return e+baseOxyEnergy+clashEnergy;
	}


	//double fastBuildPho(int seqID);

	//double phoEnergy(int seqID, bool verbose); //old method

	double phoEnergy2(int seqID, bool verbose); //updated method, improper dependent dihedral energy

	//double phoEnergy3(int seqID); //pho energy is backbone dependent

	void rebuildPho(bool verbose){
		for(int i=0;i<seqLen;i++){
			phoEnergy2(i, verbose);
			this->nodes[i]->pho = this->nodes[i]->phoTmp;
		}
	}

	void randomInit();

	void updateEnergy(){
		bool verbose = false;
		int i,j, pi, pj, sep;
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<seqLen;j++){
				pi = i*seqLen+j;
				allBaseRiboseE[pi] = 0;
				allBasePhoE[pi] = 0;
				allRiboseRiboseE[pi] = 0;
				allRibosePhoE[pi] = 0;
				allPhoPhoE[pi] = 0;
			}
		}

		for(i=0;i<seqLen;i++){
			BRNode* nodeA = nodes[i];
			for(j=i+1;j<seqLen;j++){
				BRNode* nodeB = nodes[j];
				pi = nodeA->seqID*seqLen+nodeB->seqID;
				pj = nodeB->seqID*seqLen+nodeA->seqID;


				allBaseRiboseE[pi] = getBaseRiboseEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allBaseRiboseE[pj] = getBaseRiboseEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allBasePhoE[pi] = getBasePhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allBasePhoE[pj] = getBasePhoEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allRiboseRiboseE[pi] = getRiboseRiboseEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allRiboseRiboseE[pj] = allRibosePhoE[pi];

				allRibosePhoE[pi] = getRibosePhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allRibosePhoE[pj] = getRibosePhoEnergyBM(nodeB, nodeA, sepTable[pj], et, verbose);

				allPhoPhoE[pi] = getPhoPhoEnergyBM(nodeA, nodeB, sepTable[pi], et, verbose);
				allPhoPhoE[pj] = allPhoPhoE[pi];
			}
		}

		for(i=0;i<seqLen;i++){
			allRotE[i] = nodes[i]->rot->energy;
		}

		for(i=0;i<seqLen;i++){
			allRcE[i] = nodes[i]->pho.e;
		}

		for(i=0;i<seqLen;i++){
			tmpRotE[i] = allRotE[i];
			tmpRcE[i] = allRcE[i];
			for(j=0;j<seqLen;j++){
				tmpBaseRiboseE[i*seqLen+j] = allBaseRiboseE[i*seqLen+j];
				tmpBasePhoE[i*seqLen+j] = allBasePhoE[i*seqLen+j];
				tmpRiboseRiboseE[i*seqLen+j] = allRiboseRiboseE[i*seqLen+j];
				tmpRibosePhoE[i*seqLen+j] = allRibosePhoE[i*seqLen+j];
				tmpPhoPhoE[i*seqLen+j] = allPhoPhoE[i*seqLen+j];
			}
		}
	}

	double mutEnergy(BRNode* node, bool verbose);
	void updateTmpRotamer(BRNode* node, BaseRotamer* rot, bool verbose);
	void clearTmpRotamer(BRNode* node, bool verbose);
	void acceptTmpRotamer(BRNode* node, bool verbose);


	void printPhoEnergy();
	void printPhoTmpEnergy();
	BRTreeInfo* toTreeInfo();
	double totEnergy(bool verbose);

	void debug();

	double rms(){
		vector<XYZ> backboneCoords;
		for(int i=0;i<seqLen;i++){
			for(int j=0;j<8;j++){
				backboneCoords.push_back(nodes[i]->riboAtomCoords[j]);
			}
			if(connectToDownstream[i]){
				backboneCoords.push_back(nodes[i]->pho.tList[0]);
				backboneCoords.push_back(nodes[i]->pho.tList[1]);
			}
		}

		double rms = 0;
		for(int i=0;i<initBackboneAtomList.size();i++){
			rms += squareDistance(initBackboneAtomList[i], backboneCoords[i]);
		}
		return sqrt(rms/initBackboneAtomList.size());
	}

	void printDiheds(const string& outfile);
	void printEnergyDetail();
	void runMC();
	double fragMC();

	void checkTotalEnergy();
	void checkTotalEnergyTmp();

	inline double getBaseRiboseEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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
		return baseOxyEnergy*para->wtBaseOxygen + clashEnergy;
	}


	inline double getBaseRiboseEnergyBMTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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
		return baseOxyEnergy*para->wtBaseOxygen + clashEnergy;
	}

	inline double getBasePhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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
		return baseOxyEnergy + clashEnergy;

	}

	inline double getBasePhoEnergyBMTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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
		return baseOxyEnergy + clashEnergy;
	}

	inline double getRiboseRiboseEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
		if(nodeA->seqID > nodeB->seqID){
			return getRiboseRiboseEnergyBM(nodeB, nodeA, -sep, et, verbose);
		}
		int i,j;
		if(sep == 0){
			cout << "sep0" << endl;
			return 0;
		}
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

	inline double getRiboseRiboseEnergyBMTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
		int i,j;
		if(nodeA->seqID > nodeB->seqID) {
			return getRiboseRiboseEnergyBMTmp(nodeB, nodeA, -sep, et, verbose);
		}
		if(sep == 0){
			cout << "sep0" << endl;
			return 0;
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

	inline double getRibosePhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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

	inline double getRibosePhoEnergyBMTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
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


	inline double getPhoPhoEnergyBM(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){

		if(!nodeA->connectToNeighbor || !nodeB->connectToNeighbor) return 0;

		if(nodeA->seqID > nodeB->seqID){
			return getPhoPhoEnergyBM(nodeB, nodeA, -sep, et, verbose);
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

	inline double getPhoPhoEnergyBMTmp(BRNode* nodeA, BRNode* nodeB, int sep, EnergyTable2* et, bool verbose){
		if(nodeA->seqID > nodeB->seqID) {
			return getPhoPhoEnergyBMTmp(nodeB, nodeA, -sep, et, verbose);
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

	virtual ~BackboneModelingTemplate();
};

} /* namespace NSPforcefield */

#endif /* PRED_BACKBONEMODELINGTEMPLATE_H_ */
