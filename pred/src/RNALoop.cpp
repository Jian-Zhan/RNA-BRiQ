/*
 * RNALoop.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: s2982206
 */

#include "pred/RNALoop.h"

namespace NSPpred {

RNALoop::RNALoop(const string& pdbFile, EnergyTable* et){
	this->ncg = new NeighborConnectionGenerator(1);
	this->et = et;

	RNAPDB pdb(pdbFile, "XXXX");
	RNAChain* rc = pdb.getFirstChain();
	vector<RNABase*> baseList = rc->getBaseList();
	this->length = baseList.size();
	this->typeList = new int[length];
	this->nodeList = new BaseNode*[length];
	this->moveList = new NeighborMove*[length-1];
	RNABaseName rn;
	for(int i=0;i<length;i++){
		this->typeList[i] = baseList[i]->baseTypeInt;
		LocalFrame cs = baseList[i]->getCoordSystem();
		this->nodeList[i] = new BaseNode(cs, typeList[i], i);
	}
	for(int i=0;i<length-1;i++) {
		LocalFrame cs1 = baseList[i]->getCoordSystem();
		LocalFrame cs2 = baseList[i+1]->getCoordSystem();
		BaseDistanceMatrix dm(cs1, cs2);
		CsMove move = cs1.getMove(cs2);
		XYZ p1 = baseList[i]->getAtom("P")->coord;
		XYZ p2 = baseList[i+1]->getAtom("P")->coord;
		double ene = et->getBaseBaseEnergyNeighbor(dm, typeList[i], typeList[i+1]);
		NeighborMove* nbMove = new NeighborMove(move, p1, p2, ene);
		this->moveList[i] = nbMove;
	}
}

RNALoop::RNALoop(const string& seq, NeighborConnectionGenerator* ncg, EnergyTable* et) {
	this->length = seq.length();
	map<char, int> typeMap;
	typeMap['A'] = 0;
	typeMap['U'] = 1;
	typeMap['G'] = 2;
	typeMap['C'] = 3;

	this->typeList = new int[length];
	this->nodeList = new BaseNode*[length];
	this->moveList = new NeighborMove*[length-1];

	for(int i=0;i<length;i++){
		typeList[i] = typeMap[seq[i]];
	}

	LocalFrame cs;
	NeighborMove* mv = ncg->standardMove();
	nodeList[0] = new BaseNode(cs, typeList[0], 0);
	for(int i=1;i<length;i++){
		LocalFrame csi = cs.add(mv->move);
		nodeList[i] = new BaseNode(csi, typeList[i], i);
	}
	for(int i=0;i<length-1;i++){
		moveList[i] = mv;
	}
	this->ncg = ncg;
	this->et = et;
}

void RNALoop::sampleLoopWithTerminalConstrain(BaseDistanceMatrix& dm, RNALoop& ref) {
	double minD = 999.99;
	double minE = 999.0;
	srand(time(NULL));
	int moveNum = length-1;
	int i,j;

	int ofIndex = 1;
	for(i=0;i<100000000;i++){
		int randPos = rand()%moveNum;
		int pairType = typeList[randPos]*4 + typeList[randPos+1];
		NeighborMove* nm = ncg->getRandomMove(pairType);
		moveList[randPos] = nm;
		for(j=randPos+1;j<length;j++){
			LocalFrame lf = nodeList[j-1]->cs.add(moveList[j-1]->move);
			nodeList[j]->updateCs(lf);
		}
		LocalFrame lf1 = nodeList[0]->cs;
		LocalFrame lf2 = nodeList[length-1]->cs;
		BaseDistanceMatrix x(lf1, lf2);
		double d = x.distanceTo(dm);
		if(d < minD) {
			minD = d;
			printf("step: %9d minD: %5.3f\n", i, d);
		}



		if(d < 1.0) {
			double rms = rmsdTo(ref);
			double e = getTotalEnergy();
			if(e < minE){
				minE = e;
				printf("d: %5.3f r: %5.3f e: %7.3f\n",d,rms,minE);
				char ss[200];
				sprintf(ss, "/export/home/s2982206/aspen/cppWorkplace/RNAModeling/demo/pred-%d.pdb",ofIndex);
				ofIndex++;
				string fileName = string(ss);
				printLoop(fileName);
			}
		}
	}
}

double RNALoop::getTotalEnergy(){
	double nbEnergy = 0.0;
	for(int i=0;i<length-1;i++) {
		nbEnergy += moveList[i]->energy;
	}

	double bpEnergy = 0.0;
	for(int i=0;i<length-1;i++) {
		for(int j=i+2;j<length-1;j++) {
			BaseDistanceMatrix dm(nodeList[i]->cs, nodeList[j]->cs);
			double centerDist = nodeList[i]->baseCenter.distance(nodeList[j]->baseCenter);
			bpEnergy += et->getBaseBaseEnergy(dm, typeList[i], typeList[j], centerDist);
		}
	}
	return nbEnergy + bpEnergy;
}

double RNALoop::rmsdTo(RNALoop& other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;
	RnaAtomLib atLib;
	if(length != other.length) return -1;

	for(int i=0;i<length;i++){

		vector<Atom*> aList = nodeList[i]->toAtomList(atLib);
		vector<Atom*> bList = other.nodeList[i]->toAtomList(atLib);
		if(aList.size() != bList.size()) return -1;
		for(int j=0;j<aList.size();j++){
			XYZ a = aList[j]->coord;
			XYZ b = bList[j]->coord;
			tList1.push_back(a);
			tList2.push_back(b);
		}
	}




	double r = rmsd(tList1, tList2);
	if(r > 100) {
		return 10.0;
		double r2 = rmsd2(tList1, tList2);
		cout << r << " " << r2 << endl;
		exit(1);
	}
	return rmsd(tList1, tList2);
}

void RNALoop::printLoop(const string& file) {
	cout << "print loop: " << endl;
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	for(int i=0;i<length;i++) {
		sprintf(ss, "%d", i+1);
		RNABase* base = new RNABase(string(ss), 'A', s[typeList[i]]);
		vector<Atom*> aList = nodeList[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
	}
	ofstream of;
	of.open(file, ios::out);
	rc.printPDBFormat(of, 1);
	of.close();
}

RNALoop::~RNALoop() {
	delete typeList;
	delete nodeList;
	delete moveList;
}

} /* namespace NSPpred */
