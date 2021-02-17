/*
 * BackboneModelingTemplate.cpp
 *
 *  Created on: Sep 30, 2019
 *      Author: s2982206
 */

#include <pred/BackboneModelingTemplate.h>

namespace NSPpred {

BackboneModelingTemplate::BackboneModelingTemplate(const string& inputPDB, const string& paraFile) {

	RNAPDB pdb(inputPDB, "xxxx");
	vector<RNABase*> baseList0 = pdb.getBaseList();
	vector<RNABase*> baseList;
	for(int i=0;i<baseList0.size();i++){
		if(baseList0[i]->backboneComplete())
			baseList.push_back(baseList0[i]);
	}

	for(int i=0;i<baseList.size();i++){
		baseList[i]->baseSeqID = i;
	}

	this->seqLen = baseList.size();
	this->nodes = new BRNode*[seqLen];
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	cout << "para" << endl;

	this->para = new Parameter(paraFile);

	cout << "et" << endl;
	this->et = new EnergyTable2(paraFile);

	cout << "rotLib" << endl;
	this->rotLib = new BaseRotamerLib();

	cout << "pb" << endl;
	this->pb = new PO3Builder(para, et->atET, et->roET);

	cout << "a" << endl;
	for(int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		BaseRotamer* rot = rotLib->getNearestRotamer(baseList[i]);
		this->nodes[i] = new BRNode(baseList[i], rot);
		this->nodes[i]->updateRotamer(rotLib->getRandomRotamerLv1(seq[i]));
		this->connectToDownstream[i] = false;
		nodes[i]->connectToNeighbor = false;
	}

	connectToDownstream[seqLen-1] = false;
	nodes[seqLen-1]->connectToNeighbor = false;

	this->sepTable = new int[seqLen*seqLen];
	this->allBaseRiboseE = new double[seqLen*seqLen];
	this->tmpBaseRiboseE = new double[seqLen*seqLen];
	this->allBasePhoE = new double[seqLen*seqLen];
	this->tmpBasePhoE = new double[seqLen*seqLen];
	this->allRiboseRiboseE = new double[seqLen*seqLen];
	this->tmpRiboseRiboseE = new double[seqLen*seqLen];
	this->allRibosePhoE = new double[seqLen*seqLen];
	this->tmpRibosePhoE = new double[seqLen*seqLen];
	this->allPhoPhoE = new double[seqLen*seqLen];
	this->tmpPhoPhoE = new double[seqLen*seqLen];
	this->allRotE = new double[seqLen];
	this->tmpRotE = new double[seqLen];
	this->allRcE = new double[seqLen];
	this->tmpRcE = new double[seqLen];

	cout << "b" << endl;

	for(int i=0;i<seqLen-1;i++){
		if(baseList[i]->connectToNeighbor(*baseList[i+1])){
			this->connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
	}

	cout << "c" << endl;
	vector<BaseRotamer*> rotList;
	for(int i=0;i<seqLen;i++){
		BaseRotamer* rot = new BaseRotamer(baseList[i]);
		rotList.push_back(rot);
		LocalFrame cs = nodes[i]->cs1;


		for(int j=0;j<8;j++){
			this->initBackboneAtomList.push_back(local2global(cs, rot->tList1[j]));
		}
		if(connectToDownstream[i]){
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("P")->coord);
			this->initBackboneAtomList.push_back(baseList[i+1]->getAtom("O5'")->coord);
		}
	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			int ij = i*seqLen+j;
			if(i==j) sepTable[ij] = 0;
			else if(j == i+1 && connectToDownstream[i]) sepTable[ij] = 1;
			else if(j == i-1 && connectToDownstream[j]) sepTable[ij] = -1;
			else if(j == i+2 && connectToDownstream[i] && connectToDownstream[i+1]) sepTable[ij] = 2;
			else if(j == i-2 && connectToDownstream[j] && connectToDownstream[j+1]) sepTable[ij] = -2;
			else sepTable[ij] = 3;
		}
	}

	for(int i=0;i<seqLen;i++){
		nodes[i]->updateChildInfo(nodes, seqLen);

		for(int j=0;j<nodes[i]->phoGroupC.size();j++){
			cout << "base: " << i << " pho: " << nodes[i]->phoGroupC[j] << endl;
		}
		for(int j=0;j<nodes[i]->riboseGroupC.size();j++){
			cout << "base: " << i << " ribo: " << nodes[i]->riboseGroupC[j] << endl;
		}
	}

	cout << "d" << endl;

	XYZ c2,c3,o3,p,o5,c5,c4,o4;
	for(int i=0;i<seqLen;i++){

		if(connectToDownstream[i]){
			vector<double> diheds;
			diheds.push_back(rotList[i]->chi);
			diheds.push_back(rotList[i]->improper);
			c2 = baseList[i]->getAtom("C2'")->coord;
			c3 = baseList[i]->getAtom("C3'")->coord;
			o3 = baseList[i]->getAtom("O3'")->coord;
			p = baseList[i+1]->getAtom("P")->coord;
			o5 = baseList[i+1]->getAtom("O5'")->coord;
			c5 = baseList[i+1]->getAtom("C5'")->coord;
			c4 = baseList[i+1]->getAtom("C4'")->coord;
			o4 = baseList[i+1]->getAtom("O4'")->coord;

			diheds.push_back(dihedral(c2, c3, o3, p));
			diheds.push_back(dihedral(c3, o3, p, o5));
			diheds.push_back(dihedral(o3, p, o5, c5));
			diheds.push_back(dihedral(p, o5, c5, c4));
			diheds.push_back(dihedral(o5, c5, c4, o4));
			diheds.push_back(rotList[i+1]->improper);
			this->initDihedsList.push_back(diheds);
		}

	}
	cout << "e" << endl;
	for(int i=0;i<seqLen;i++){
		//phoEnergy(i);
		phoEnergy2(i, false);
		acceptTmpRotamer(nodes[i], false);
	}


	updateEnergy();
}

void BackboneModelingTemplate::updateDiheds(){
	XYZ c2,c3,o3,p,o5,c5,c4,o4;
	this->predDihedsList.clear();
	for(int i=0;i<seqLen-1;i++){
		if(!connectToDownstream[i]) continue;
		BRNode* nodeA = nodes[i];
		BRNode* nodeB = nodes[i+1];
		vector<double> diheds;
		diheds.push_back(nodeA->rot->chi);
		diheds.push_back(nodeA->rot->improper);
		c2 = nodeA->riboAtomCoords[1];
		c3 = nodeA->riboAtomCoords[2];
		o3 = nodeA->riboAtomCoords[6];
		p = nodeA->pho.tList[0];
		o5 = nodeA->pho.tList[1];
		c5 = nodeB->riboAtomCoords[7];
		c4 = nodeB->riboAtomCoords[3];
		o4 = nodeB->riboAtomCoords[4];
		diheds.push_back(dihedral(c2, c3, o3, p));
		diheds.push_back(dihedral(c3, o3, p, o5));
		diheds.push_back(dihedral(o3, p, o5, c5));
		diheds.push_back(dihedral(p, o5, c5, c4));
		diheds.push_back(dihedral(o5, c5, c4, o4));
		diheds.push_back(nodes[i+1]->rot->improper);
		predDihedsList.push_back(diheds);
	}
}

double BackboneModelingTemplate::phoEnergy2(int seqID, bool verbose) {
	if(seqID < 0) return 0;
	if(seqID >= seqLen-1) return 0;
	if(!connectToDownstream[seqID]) return 0;
	if(verbose) {
		cout << "update pho: " << seqID << endl;
	}
	BRNode* nodeA = nodes[seqID];
	BRNode* nodeB = nodes[seqID+1];
	PhophateGroupLocal pl = this->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
	nodeA->phoLocalTmp = pl;
	nodeA->phoTmp = PhophateGroup(pl, nodeA->tmpCs2);
	return nodeA->phoTmp.e;
}


/*
double BackboneModelingTemplate::phoEnergy3(int seqID) {
	if(seqID < 0) return 0;
	if(seqID >= seqLen-1) return 0;
	if(!connectToDownstream[seqID]) return 0;
	BRNode* nodeA = nodes[seqID];
	BRNode* nodeB = nodes[seqID+1];

	XYZ center = (nodeA->riboAtomCoordsTmp[6] + nodeB->riboAtomCoordsTmp[7])*0.5;

	vector<int> baseTypeList;
	vector<LocalFrame> baseCs1List;
	vector<XYZ*> riboCoordsList;
	vector<XYZ*> phoCoordsList;

	for(int i=0;i<seqLen;i++){
		if(abs(i - seqID) < 2) continue;
		if(center.distance(nodes[i]->baseAtomCoordsTmp[0]) < 10){
			baseTypeList.push_back(nodes[i]->baseType);
			baseCs1List.push_back(nodes[i]->tmpCs1);
		}

		if(center.distance(nodes[i]->tmpCs1.origin_) < 8){
			riboCoordsList.push_back(nodes[i]->riboAtomCoordsTmp);
		}

		if(center.distance(nodes[i]->phoTmp.tList[0]) < 8) {
			phoCoordsList.push_back(nodes[i]->phoTmp.tList);
		}
	}

	PhophateGroupLocal pl = this->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper, baseTypeList, baseCs1List, riboCoordsList, phoCoordsList);
	nodeA->phoLocalTmp = pl;
	nodeA->phoTmp = PhophateGroup(pl, nodeA->tmpCs2);
	return nodeA->phoTmp.e;
}
*/

/*
double BackboneModelingTemplate::phoEnergy(int seqID, bool verbose){
	if(seqID < 0) return 0;
	if(seqID >= seqLen-1) return 0;
	if(!connectToDownstream[seqID]) return 0;
	BRNode* nodeA = nodes[seqID];
	BRNode* nodeB = nodes[seqID+1];
	LocalFrame cs = nodeA->tmpCs2;

	double len1 = 1.605;
	double len2 = 1.592;
	double len3 = 1.422;
	double ang1 = 120.1;
	double ang2 = 103.5;
	double ang3 = 120.7;
	double ang4 = 111.1;

	XYZ c2,o2,c3,o3,p,op1,op2,o5,c5,c4,o4;
	c2 = nodeA->riboAtomCoordsTmp[1];
	o2 = nodeA->riboAtomCoordsTmp[5];
	c3 = nodeA->riboAtomCoordsTmp[2];
	o3 = nodeA->riboAtomCoordsTmp[6];
	c5 = nodeB->riboAtomCoordsTmp[7];
	c4 = nodeB->riboAtomCoordsTmp[3];
	o4 = nodeB->riboAtomCoordsTmp[4];

	int d1d2LibSize = pb->d1d2Lib1A.size();

	int d1Index,d2Index, d1d2Index;
	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;
	double u, e;

	double minE = 99999999.9;
	int bestDihed1, bestDihed2, bestIndex1, bestIndex2, localBestD1, localBestD2;
	int dihed1, dihed2;

	for(int i=0;i<d1d2LibSize;i++){
		d1Index = (int)pb->d1d2Lib1A[i].x_;
		d2Index = (int)pb->d1d2Lib1A[i].y_;
		d1d2Index = d1Index*360+d2Index;
		p = local2global(cs, pb->pList[d1d2Index]);
		op1 = local2global(cs,pb->op1List[d1d2Index]);
		op2 = local2global(cs,pb->op2List[d1d2Index]);
		o5 = local2global(cs,pb->o5List[d1d2Index]);
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = calPhoEnergy(xd3, xang3, xang4, d1Index, d2Index, (int)xdihed3, (int)xdihed4, (int)xdihed5, p, o5, op1, op2, seqID);
		if(e < minE){
			bestDihed1 = d1Index;
			bestDihed2 = d2Index;
			bestIndex1 = i;
			minE = e;
		}
	}

	for(int i=0;i<d1d2LibSize;i++){
		d1Index = (int)pb->d1d2Lib2A[bestIndex1][i].x_;
		d2Index = (int)pb->d1d2Lib2A[bestIndex1][i].y_;
		d1d2Index = d1Index*360+d2Index;
		p = local2global(cs, pb->pList[d1d2Index]);
		op1 = local2global(cs,pb->op1List[d1d2Index]);
		op2 = local2global(cs,pb->op2List[d1d2Index]);
		o5 = local2global(cs,pb->o5List[d1d2Index]);
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = calPhoEnergy(xd3, xang3, xang4, d1Index, d2Index, (int)xdihed3, (int)xdihed4, (int)xdihed5, p, o5, op1, op2, seqID);
		if(e < minE){
			bestDihed1 = d1Index;
			bestDihed2 = d2Index;
			bestIndex2 = i;
			minE = e;
		}
	}

	if(bestDihed1 < 220){
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<9;i++){
			dihed1 = localBestD1 + (i-4);
			if(dihed1 < 0) dihed1 = dihed1 +360;
			if(dihed1 > 359) dihed1 = dihed1 - 360;

			for(int j=0;j<9;j++){
				dihed2 = localBestD2 + (j-4);
				if(dihed2 < 0) dihed2 = dihed2 +360;
				if(dihed2 > 359) dihed2 = dihed2 - 360;
				d1d2Index = dihed1*360+dihed2;
				p = local2global(cs, pb->pList[d1d2Index]);
				op1 = local2global(cs,pb->op1List[d1d2Index]);
				op2 = local2global(cs,pb->op2List[d1d2Index]);
				o5 = local2global(cs,pb->o5List[d1d2Index]);
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = calPhoEnergy(xd3, xang3, xang4, dihed1, dihed2, (int)xdihed3, (int)xdihed4, (int)xdihed5, p, o5, op1, op2, seqID);
				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
	}
	else {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<11;i++){
			dihed1 = localBestD1 + (i-5)*2;
			if(dihed1 < 0) dihed1 = dihed1 +360;
			if(dihed1 > 359) dihed1 = dihed1 - 360;

			for(int j=0;j<11;j++){
				dihed2 = localBestD2 + (j-5)*2;
				if(dihed2 < 0) dihed2 = dihed2 +360;
				if(dihed2 > 359) dihed2 = dihed2 - 360;
				d1d2Index = dihed1*360+dihed2;
				p = local2global(cs, pb->pList[d1d2Index]);
				op1 = local2global(cs,pb->op1List[d1d2Index]);
				op2 = local2global(cs,pb->op2List[d1d2Index]);
				o5 = local2global(cs,pb->o5List[d1d2Index]);
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = calPhoEnergy(xd3, xang3, xang4, dihed1, dihed2, (int)xdihed3, (int)xdihed4, (int)xdihed5, p, o5, op1, op2, seqID);
				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}
	}

	nodeA->phoTmp.e = minE;
	d1d2Index = bestDihed1*360+bestDihed2;
	nodeA->phoTmp.tList[0] = local2global(cs, pb->pList[d1d2Index]);
	nodeA->phoTmp.tList[1] = local2global(cs,pb->o5List[d1d2Index]);
	nodeA->phoTmp.tList[2] = local2global(cs,pb->op1List[d1d2Index]);
	nodeA->phoTmp.tList[3] = local2global(cs,pb->op2List[d1d2Index]);

	//cout << minE << endl;
	return minE;
}
*/


/*
double BackboneModelingTemplate::fastBuildPho(int seqID){
	if(seqID < 0) return 0;
	if(seqID >= seqLen-1) return 0;
	if(!connectToDownstream[seqID]) return 0;
	BRNode* nodeA = nodes[seqID];
	BRNode* nodeB = nodes[seqID+1];
	double len = nodeA->riboAtomCoordsTmp[6].distance(nodeB->riboAtomCoordsTmp[7]);
	int index = pb->move2index.toIndex(nodeA->tmpCs2, nodeB->tmpCs3);
	int d1d2Index = pb->dihed1List[index]*360+pb->dihed2List[index];
	LocalFrame cs = nodeA->tmpCs2;
	nodeA->phoTmp.tList[0] = local2global(cs, pb->pList[d1d2Index]);
	nodeA->phoTmp.tList[1] = local2global(cs,pb->o5List[d1d2Index]);
	nodeA->phoTmp.tList[2] = local2global(cs,pb->op1List[d1d2Index]);
	nodeA->phoTmp.tList[3] = local2global(cs,pb->op2List[d1d2Index]);
	double e = pb->e2List[index]*para->wtConnect;
	if(len < 2.4)
		e = e + (2.4-len)*(2.4-len)*20;
	if(len > 4.4)
		e = e + (len-4.4)*(len-4.4)*20;
	nodeA->phoTmp.e = e;
	return e;
}
*/

void BackboneModelingTemplate::updateTmpRotamer(BRNode* node, BaseRotamer* rot, bool verbose){
	node->tmpRot = rot;
	for(int i=0;i<11;i++){
		node->riboAtomCoordsTmp[i] = local2global(node->tmpCs1, node->tmpRot->tList1[i]);
	}
	node->tmpCs2 = node->tmpCs1 + rot->mv12;

	int chainBreakNum = node->phoGroupC.size();
	int indexA, indexB;

	BRNode* nodeA;
	BRNode* nodeB;
	for(int i=0;i<chainBreakNum;i++){
		if(verbose){
			cout << "connection: " << indexA << " " << indexB << endl;
		}
		indexA = node->phoGroupC[i];
		indexB = indexA+1;
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];


		nodeA->phoLocalTmp = pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
	}
}

void BackboneModelingTemplate::clearTmpRotamer(BRNode* nodeSelect, bool verbose){
	int i,j;
	nodeSelect->tmpRot = nodeSelect->rot;
	nodeSelect->tmpCs2 = nodeSelect->cs2;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoordsTmp[j] = nodeSelect->riboAtomCoords[j];
	}
	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	/*
	int pi, pj;

	int id = nodeSelect->seqID;
	tmpRotE[id] = allRotE[id];
	tmpRcE[id] = allRcE[id];
	for(int j=id-2;j<id+3;j++){
		if(id < 0) continue;
		if(id >= seqLen) continue;
		pi = id*seqLen+j;
		pj = j*seqLen+id;
		tmpBaseRiboseE[pj] = allBaseRiboseE[pj];
		tmpBasePhoE[pj] = allBasePhoE[pj];
		tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
		tmpRiboseRiboseE[pj] = allRiboseRiboseE[pj];
		tmpRibosePhoE[pi] = allRibosePhoE[pi];
		tmpRibosePhoE[pj] = allRibosePhoE[pj];
		tmpPhoPhoE[pi] = allPhoPhoE[pi];
		tmpPhoPhoE[pj] = allPhoPhoE[pj];
	}

	if(id > 0 && connectToDownstream[id-1]){
		id = id-1;
		tmpRcE[id] = allRcE[id];
		for(int j=id-2;j<id+3;j++){
			if(id < 0) continue;
			if(id >= seqLen) continue;
			pi = id*seqLen+j;
			pj = j*seqLen+id;

			tmpBasePhoE[pj] = allBasePhoE[pj];
			tmpRibosePhoE[pj] = allRibosePhoE[pj];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
			tmpPhoPhoE[pj] = allPhoPhoE[pj];
		}
	}
	*/
}

void BackboneModelingTemplate::acceptTmpRotamer(BRNode* nodeSelect, bool verbose){
	int i,j;
	nodeSelect->rot = nodeSelect->tmpRot;
	nodeSelect->cs2 = nodeSelect->tmpCs2;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoords[j] = nodeSelect->riboAtomCoordsTmp[j];
	}
	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}
	/*
	int pi, pj;

	int id = nodeSelect->seqID;
	allRotE[id] = tmpRotE[id];
	allRcE[id] = tmpRcE[id];
	for(int j=id-2;j<id+3;j++){
		if(id < 0) continue;
		if(id >= seqLen) continue;
		pi = id*seqLen+j;
		pj = j*seqLen+id;
		allBaseRiboseE[pj] = tmpBaseRiboseE[pj];
		allBasePhoE[pj] = tmpBasePhoE[pj];
		allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
		allRiboseRiboseE[pj] = tmpRiboseRiboseE[pj];
		allRibosePhoE[pi] = tmpRibosePhoE[pi];
		allRibosePhoE[pj] = tmpRibosePhoE[pj];
		allPhoPhoE[pi] = tmpPhoPhoE[pi];
		allPhoPhoE[pj] = tmpPhoPhoE[pj];
	}

	if(id > 0 && connectToDownstream[id-1]){
		id = id-1;
		allRcE[id] = tmpRcE[id];
		for(int j=id-2;j<id+3;j++){
			if(id < 0) continue;
			if(id >= seqLen) continue;
			pi = id*seqLen+j;
			pj = j*seqLen+id;

			allBasePhoE[pj] = tmpBasePhoE[pj];
			allRibosePhoE[pj] = tmpRibosePhoE[pj];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
			allPhoPhoE[pj] = tmpPhoPhoE[pj];
		}
	}
	*/
}

double BackboneModelingTemplate::mutEnergy(BRNode* node, bool verbose){

	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	double mutE = 0;

	//int baseGroupANum = node->baseGroupA.size();


	int riboseGroupANum = node->riboseGroupA.size();
	int riboseGroupCNum = node->riboseGroupC.size();

	int phoGroupANum = node->phoGroupA.size();
	int phoGroupCNum = node->phoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		mutE += (nodes[indexX]->tmpRot->energy - nodes[indexX]->rot->energy);

	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		mutE += (nodes[indexX]->phoLocalTmp.e - nodes[indexX]->pho.e);
	}

	/*
	 * base-ribose energy
	 */

	//A-C
	for(i=0;i<seqLen;i++){
		indexX = i;
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getBaseRiboseEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getBaseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * base-pho energy
	 */

	//A-C
	for(i=0;i<seqLen;i++){
		indexX = i;
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getBasePhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getBasePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = node->riboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
			//	printf("mut: %-2d %-2d %2d %7.3f %7.3f %7.3f\n", indexX, indexY, sepTable[pi], getRiboseRiboseEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose), getRiboseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose), getRiboseRiboseEnergyBM(nodeY, nodeX, -sepTable[pi], et, verbose));

				mutE += (getRiboseRiboseEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getRiboseRiboseEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}


	/*
	 * ribose-pho energy
	 */

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = node->riboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-A
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = node->phoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getRibosePhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getRibosePhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	/*
	 * pho-pho energy
	 */

	//A-C
	for(i=0;i<phoGroupANum;i++){
		indexX = node->phoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getPhoPhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getPhoPhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	//C-C
	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			if(squareDistance(nodeX->cs1.origin_, nodeY->cs1.origin_) < 144){
				mutE += (getPhoPhoEnergyBMTmp(nodeX, nodeY, sepTable[pi], et, verbose) - getPhoPhoEnergyBM(nodeX, nodeY, sepTable[pi], et, verbose));
			}
		}
	}

	return mutE;
}



/*
void BackboneModelingTemplate::printPhoEnergy(){
	for(int i=0;i<seqLen-1;i++){

		BRNode* nodeA = nodes[i];
		BRNode* nodeB = nodes[i+1];
		LocalFrame cs = nodeA->cs2;

		double len1 = 1.605;
		double len2 = 1.592;
		double len3 = 1.422;
		double ang1 = 120.1;
		double ang2 = 103.5;
		double ang3 = 120.7;
		double ang4 = 111.1;

		XYZ c2,o2,c3,o3,p,op1,op2,o5,c5,c4,o4;
		c2 = nodeA->riboAtomCoords[1];
		o2 = nodeA->riboAtomCoords[5];
		c3 = nodeA->riboAtomCoords[2];
		o3 = nodeA->riboAtomCoords[6];
		p = nodeA->pho.tList[0];
		o5 = nodeA->pho.tList[1];
		op1 = nodeA->pho.tList[2];
		op2 = nodeA->pho.tList[3];
		c5 = nodeB->riboAtomCoords[7];
		c4 = nodeB->riboAtomCoords[3];
		o4 = nodeB->riboAtomCoords[4];

		double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5, xdihed1, xdihed2;
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed1 = dihedral(c2,c3,o3,p);
		xdihed2 = dihedral(c3,o3,p, o5);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed1 < 0) xdihed1 += 360;
		if(xdihed2 < 0) xdihed2 += 360;
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		int x1 = (int)xdihed1;
		int x2 = (int)xdihed2;
		int x3 = (int)xdihed3;
		int x4 = (int)xdihed4;
		int x5 = (int)xdihed5;
		double e = calPhoEnergy(xd3, xang3, xang4, x1, x2, x3, x4, x5, p, o5, op1, op2, i);
		printf("%4.2f %5.1f %5.1f %6.2f %6.2f %6.2f %6.2f %6.2f ", xd3, xang3, xang4, xdihed1, xdihed2, xdihed3, xdihed4, xdihed5);

		printf("%8.3f %8.3f\n", e, this->nodes[i]->pho.e);
	}
	cout << endl;
}
*/

/*
void BackboneModelingTemplate::printPhoTmpEnergy(){
	for(int i=0;i<seqLen-1;i++){

		BRNode* nodeA = nodes[i];
		BRNode* nodeB = nodes[i+1];
		LocalFrame cs = nodeA->cs2;

		double len1 = 1.605;
		double len2 = 1.592;
		double len3 = 1.422;
		double ang1 = 120.1;
		double ang2 = 103.5;
		double ang3 = 120.7;
		double ang4 = 111.1;

		XYZ c2,o2,c3,o3,p,op1,op2,o5,c5,c4,o4;
		c2 = nodeA->riboAtomCoordsTmp[1];
		o2 = nodeA->riboAtomCoordsTmp[5];
		c3 = nodeA->riboAtomCoordsTmp[2];
		o3 = nodeA->riboAtomCoordsTmp[6];
		p = nodeA->phoTmp.tList[0];
		o5 = nodeA->phoTmp.tList[1];
		op1 = nodeA->phoTmp.tList[2];
		op2 = nodeA->phoTmp.tList[3];
		c5 = nodeB->riboAtomCoordsTmp[7];
		c4 = nodeB->riboAtomCoordsTmp[3];
		o4 = nodeB->riboAtomCoordsTmp[4];

		double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5, xdihed1, xdihed2;
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed1 = dihedral(c2,c3,o3,p);
		xdihed2 = dihedral(c3,o3,p, o5);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed1 < 0) xdihed1 += 360;
		if(xdihed2 < 0) xdihed2 += 360;
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		int x1 = (int)xdihed1;
		int x2 = (int)xdihed2;
		int x3 = (int)xdihed3;
		int x4 = (int)xdihed4;
		int x5 = (int)xdihed5;

		double e = calPhoEnergy(xd3, xang3, xang4, x1, x2, x3, x4, x5, p, o5, op1, op2, i);
		printf("%4.2f %5.1f %5.1f %6.2f %6.2f %6.2f %6.2f %6.2f ", xd3, xang3, xang4, xdihed1, xdihed2, xdihed3, xdihed4, xdihed5);

		printf("%8.3f %8.3f\n", e, this->nodes[i]->phoTmp.e);
	}
	cout << endl;
}
*/

BRTreeInfo* BackboneModelingTemplate::toTreeInfo(){
	vector<BRNode*> flexNodes;
	for(int i=0;i<seqLen;i++){
		flexNodes.push_back(nodes[i]);
	}
	return new BRTreeInfo(seqLen, seq, connectToDownstream, nodes, flexNodes, 0.0);
}

double BackboneModelingTemplate::totEnergy(bool verbose){
	double tot = 0;
	int i,j;

	//cout << "total energy:" << endl;
	//ribose-ribose, pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];

			//printf("%-2d %-2d %8.3f\n", i, j, getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose));
			tot += getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			tot += getPhoPhoEnergyBM(nodeA, nodeB, sep, et, verbose);
		}
	}


	//base-ribose, base-pho, ribose-pho

	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			tot += getBaseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			tot += getBasePhoEnergyBM(nodeA, nodeB, sep,  et, verbose);
			tot += getRibosePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
		}
	}



	for(int i=0;i<seqLen;i++){
		if(verbose){
			cout << "rotEnergy: " << i << " " << nodes[i]->rot->energy << endl;
		}
		tot += nodes[i]->rot->energy;
	}


	for(int i=0;i<seqLen;i++){
		if(verbose){
			cout << "phoEnergy: " << i << " " <<  nodes[i]->pho.e << endl;
		}
		tot += nodes[i]->pho.e;
	}

	return tot;
}

void BackboneModelingTemplate::printEnergyDetail(){
	double e = 0;

	for(int i=0;i<seqLen;i++){
		e += nodes[i]->rot->energy;
		e += nodes[i]->pho.e;
//		printf("Node: %2d RotamerE: %9.3f PhoE: %9.3f\n", i, nodes[i]->rot->energy, nodes[i]->pho.e);
	}
	bool verbose = false;

	for(int i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(int j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int pi = i*seqLen+j;
			int pj = j*seqLen+i;
			int sepi = sepTable[pi];
			int sepj = sepTable[pj];
			{
				double tot = 0;
				double bb = 0;
				double br = getBaseRiboseEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double rb = getBaseRiboseEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double bp = getBasePhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double pb = getBasePhoEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double rr = getRiboseRiboseEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double rp = getRibosePhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				double pr = getRibosePhoEnergyBM(nodeB, nodeA, sepj, et, verbose);
				double pp = getPhoPhoEnergyBM(nodeA, nodeB, sepi, et, verbose);
				tot = bb + br + rb + bp + pb + rr + rp + pr + pp;
				if(abs(tot) > 0.001)
					printf("Pair: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
				e += tot;
			}

			{
				double tot = 0;
				double bb = 0;
				double br = getBaseRiboseEnergyBMTmp(nodeA, nodeB, sepi, et, verbose);
				double rb = getBaseRiboseEnergyBMTmp(nodeB, nodeA, sepj, et, verbose);
				double bp = getBasePhoEnergyBMTmp(nodeA, nodeB, sepi, et, verbose);
				double pb = getBasePhoEnergyBMTmp(nodeB, nodeA, sepj, et, verbose);
				double rr = getRiboseRiboseEnergyBMTmp(nodeA, nodeB, sepi, et, verbose);
				double rp = getRibosePhoEnergyBMTmp(nodeA, nodeB, sepi, et, verbose);
				double pr = getRibosePhoEnergyBMTmp(nodeB, nodeA, sepj, et, verbose);
				double pp = getPhoPhoEnergyBMTmp(nodeA, nodeB, sepi, et, verbose);
				tot = bb + br + rb + bp + pb + rr + rp + pr + pp;
				if(abs(tot) > 0.001)
					printf("TmPr: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
			}

		}
	}

	printf("total energy: %9.3f\n",  e);
	cout << endl;
}

void BackboneModelingTemplate::printDiheds(const string& outfile){
	ofstream out;
	out.open(outfile, ios::out);
	if(!out.is_open()) {
		cout << "fail to open: " << outfile << endl;
		exit(1);
	}

	for(int i=0;i<initDihedsList.size();i++){
		char xx[20];
		for(int j=0;j<8;j++){
			sprintf(xx,"%7.2f ", initDihedsList[i][j]);
			out << string(xx);
		}
		out << "   ";
		for(int j=0;j<8;j++){
			sprintf(xx,"%7.2f ", predDihedsList[i][j]);
			out << string(xx);
		}
		out << endl;
	}
	out << "RMSD " << rms() << endl;
	out.close();
}

double BackboneModelingTemplate::fragMC(){

	int stepNum = 2000;
	double curEne = totEnergy(false);
	double lastEne = curEne;
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	BaseRotamer* rotMut;

	double minEne = 9999999999.9;
	double rmsd = 999.9;
	int fragLen = 7;
	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		cout << startIndex << endl;
		for(double T=20.0;T>0.05;T=T*0.8){
			for(int k=0;k<stepNum;k++){
				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->getRandomRotamerLv1(node->baseType);
				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	stepNum = 3000;
	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		double T0 = 10.0;
		double T1 = 0.01;
		cout << startIndex << endl;
		for(double T=T0;T>T1;T=T*0.9){
			for(int k=0;k<stepNum;k++){

				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->getRandomRotamerLv1(node->baseType);

				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						updateDiheds();
						printf("R1: ene: %8.3f rmsd: %6.4f\n", curEne,rmsd);
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	for(int startIndex=0;startIndex<seqLen-fragLen+4;startIndex+=3){
		double T0 = 10.0;
		double T1 = 0.01;
		cout << startIndex << endl;
		for(double T=T0;T>T1;T=T*0.9){
			for(int k=0;k<stepNum;k++){

				randIndex = rand()%fragLen + startIndex;
				if(randIndex >= seqLen) continue;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->getRandomRotamerLv1(node->baseType);

				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						updateDiheds();
						printf("R2: ene: %8.3f rmsd: %6.4f\n", curEne,rmsd);
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		node = nodes[i];
		for(int j=0;j<200;j++){
			int type1 = node->rot->rotTypeLv1;
			rotMut = rotLib->rotLib[node->baseType][type1*200+j];
			updateTmpRotamer(node, rotMut, false);
			mutE = mutEnergy(node, false);
			if(mutE < 0) {
				curEne += mutE;
				acceptTmpRotamer(node, false);
				if(curEne < minEne){
					minEne = curEne;
					updateDiheds();
				}
			}
			else{
				clearTmpRotamer(node, false);
			}
		}
	}
	return rmsd;
}

void BackboneModelingTemplate::debug(){
	int stepNum = 100;
	double curEne = totEnergy(false);
	double lastEne = curEne;
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	BaseRotamer* rotMut;
	double T = 2.0;
	double minEne = 99999.9;

	double e = totEnergy(true);
	cout << "curE: " << curEne << " totE: " << e << endl;

	for(int k=0;k<1000;k++){
		cout << "step: " << k << endl;
		randIndex = rand()%seqLen;
		node = nodes[randIndex];
		randType = rand()%2;
		rotMut = rotLib->getRandomRotamerLv1(node->baseType);

		printEnergyDetail();

		updateTmpRotamer(node, rotMut,false);
		mutE = mutEnergy(node,false);


		if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
			cout << "accept: " << mutE << endl;
			curEne += mutE;
			acceptTmpRotamer(node, false);
			if(curEne < minEne){
				minEne = curEne;
				updateDiheds();
			}
		}
		else{
			cout <<"reject" << mutE <<  endl;
			clearTmpRotamer(node, false);
		}

		printEnergyDetail();

		/*
		checkTotalEnergy();
		checkTotalEnergyTmp();
		*/
		double totE = totEnergy(false);
		if(abs(curEne - totE) > 0.001){
			cout << "error" << endl;
		}
		printf(" T=%8.5f  addEne: %9.3f totEne: %9.3f rms: %5.3f\n", T, curEne, totE, rms());
	}
}


void BackboneModelingTemplate::runMC(){

	int stepNum = 1000*seqLen;
	double curEne = totEnergy(false);
	double lastEne = curEne;
	double mutE;
	srand(time(0));
	int randIndex, randType;
	BRNode* node;
	BaseRotamer* rotMut;

	double minEne = 9999999999.9;


	for(int n=0;n<1;n++){
		for(double T=10.0;T>0.01;T=T*0.9){
			int acNum = 0;
			for(int k=0;k<stepNum;k++){
				randIndex = rand()%seqLen;
				node = nodes[randIndex];
				randType = rand()%2;
				rotMut = rotLib->getRandomRotamerLv1(node->baseType);
				updateTmpRotamer(node, rotMut, false);
				mutE = mutEnergy(node, false);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					acceptTmpRotamer(node, false);
					acNum ++;
					if(curEne < minEne){
						minEne = curEne;
						updateDiheds();
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
			double totE = totEnergy(false);
			printf("R: %d T=%8.5f accept=%5d addEne: %9.3f totEne: %9.3f minEne: %9.3f rms: %5.3f\n", n, T, acNum, curEne, totE,minEne, rms());
		}

		for(int i=0;i<seqLen;i++){
			node = nodes[i];
			for(int j=0;j<1500;j++){
				int type1 = node->rot->rotTypeLv1;
				rotMut = rotLib->rotLib[node->baseType][j];
				updateTmpRotamer(node, rotMut,false);
				mutE = mutEnergy(node,false);
				if(mutE < 0) {
					curEne += mutE;
					acceptTmpRotamer(node, false);
					if(curEne < minEne){
						minEne = curEne;
						updateDiheds();
					}
				}
				else{
					clearTmpRotamer(node, false);
				}
			}
		}

	}
}

void BackboneModelingTemplate::checkTotalEnergy(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;



	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allBaseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBaseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//base-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBasePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allBasePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBasePhoE[i*seqLen+j], e);
			}
		}
	}

	//ribose-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRiboseRiboseEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allRiboseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allRiboseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//ribose-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRibosePhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allRibosePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allRibosePhoE[i*seqLen+j], e);
			}
		}
	}

	//pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getPhoPhoEnergyBM(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->allPhoPhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: pho: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allPhoPhoE[i*seqLen+j], e);
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->rot->energy;
		if(abs(allRotE[i] - e) > 0.001)
			cout << "DIFF rot: " << i << " " << allRotE[i] << " " << e << endl;

	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->pho.e;
		if(abs(allRcE[i] -e) > 0.001)
			cout << "DIFF connect: " << i << " " << allRcE[i] << " " << e << endl;
	}
}

void BackboneModelingTemplate::checkTotalEnergyTmp(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;

	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergyBMTmp(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpBaseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBaseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//base-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBasePhoEnergyBMTmp(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpBasePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBasePhoE[i*seqLen+j], e);
			}
		}
	}

	//ribose-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRiboseRiboseEnergyBMTmp(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpRiboseRiboseE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d ribose: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpRiboseRiboseE[i*seqLen+j], e);
			}
		}
	}

	//ribose-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getRibosePhoEnergyBMTmp(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpRibosePhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: ribose: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpRibosePhoE[i*seqLen+j], e);
			}
		}
	}

	//pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getPhoPhoEnergyBMTmp(nodeA, nodeB, sep, et, verbose);
			if(abs(e - this->tmpPhoPhoE[i*seqLen+j]) > 0.001){
				printf("DIFF: pho: %2d pho: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpPhoPhoE[i*seqLen+j], e);
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->tmpRot->energy;
		if(abs(tmpRotE[i] - e) > 0.001)
			cout << "DIFF rot: " << i << " " << tmpRotE[i] << " " << e << endl;
	}

	for(int i=0;i<seqLen;i++){
		e = nodes[i]->phoTmp.e;
		if(abs(tmpRcE[i] -e) > 0.001)
			cout << "DIFF connect: " << i << " " << tmpRcE[i] << " " << e << endl;
	}
}


BackboneModelingTemplate::~BackboneModelingTemplate() {
	for(int i=0;i<seqLen;i++){
		delete nodes[i];
	}
	delete [] nodes;
	delete [] seq;
	delete [] connectToDownstream;
	delete para;
	delete et;
	delete rotLib;

	delete [] allBaseRiboseE;
	delete [] tmpBaseRiboseE;
	delete [] allBasePhoE;
	delete [] tmpBasePhoE;
	delete [] allRiboseRiboseE;
	delete [] tmpRiboseRiboseE;
	delete [] allRibosePhoE;
	delete [] tmpRibosePhoE;
	delete [] allPhoPhoE;
	delete [] tmpPhoPhoE;
	delete [] allRotE;
	delete [] tmpRotE;
	delete [] allRcE;
	delete [] tmpRcE;
}

} /* namespace NSPforcefield */
