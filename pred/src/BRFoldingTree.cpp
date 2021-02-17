/*
 * BRFoldingTree.cpp
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#include <pred/BRFoldingTree.h>

namespace NSPpred {

BRTreeInfo::BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, double ene) {
	this->seqLen = seqLen;
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->nodes = new BRNode*[seqLen];
	this->ene = ene;
	for(int i=0;i<seqLen;i++){
		this->seq[i] = seq[i];
		this->connectToDownstream[i] = con[i];
	}
	for(int i=0;i<seqLen;i++) {
		BRNode* br = new BRNode(nodes[i]->baseType, nodes[i]->seqID);
		br->copyValueFrom(*nodes[i]);
		this->nodes[i] = br;
	}
	this->ene = ene;
	for(int i=0;i<seqLen;i++){
		this->freeNodeIDs.push_back(i);
	}
}

BRTreeInfo::BRTreeInfo(int seqLen, int* seq, bool* con, BRNode** nodes, vector<BRNode*>& flexibleNodes, double ene) {
	this->seqLen = seqLen;
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->nodes = new BRNode*[seqLen];
	this->ene = ene;
	for(int i=0;i<seqLen;i++){
		this->seq[i] = seq[i];
		this->connectToDownstream[i] = con[i];
	}
	for(int i=0;i<seqLen;i++) {
		BRNode* br = new BRNode(nodes[i]->baseType, nodes[i]->seqID);
		br->copyValueFrom(*nodes[i]);
		this->nodes[i] = br;
	}
	this->ene = ene;
	for(int i=0;i<flexibleNodes.size();i++){
		this->freeNodeIDs.push_back(flexibleNodes[i]->seqID);
	}
}
/*
BRTreeInfo::BRTreeInfo(const string& pdbFile){


	RNAPDB pdb(pdbFile, "xxxx");
	vector<RNABase*> baseList = pdb.getBaseList();

	this->seqLen = baseList.size();
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->nodes = new BRNode*[seqLen];
	this->ene = 0;

	for(int i=0;i<seqLen;i++) {
		this->nodes[i] = new BRNode(baseList[i]);
		this->seq[i] = this->nodes[i]->baseType;

		if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1])) {
			this->connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
		else {
			this->connectToDownstream[i] = false;
			nodes[i]->connectToNeighbor = false;
		}
	}
}
*/

double BRTreeInfo::rmsd(BRTreeInfo* other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;

	vector<XYZ> subList1;
	vector<XYZ> subList2;

	RnaAtomLib atLib;
	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = nodes[i]->toBaseAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = other->nodes[i]->toBaseAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->freeNodeIDs.size();i++){
		vector<Atom*> aList = nodes[freeNodeIDs[i]]->toBaseAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			subList1.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	for(unsigned int i=0;i<this->freeNodeIDs.size();i++){
		vector<Atom*> aList = other->nodes[freeNodeIDs[i]]->toBaseAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			subList2.push_back(aList[j]->coord);
		for(int k=0;k<aList.size();k++){
			delete aList[k];
		}
	}

	return subRMSD(tList1, tList2, subList1, subList2);
}

void BRTreeInfo::printPDB(ofstream& of, int modelID) {
	RNAChain rc;
	string s = "AUGC";
	char ss[20];


	sprintf(ss, "MODEL%9d", modelID);
	of << string(ss) << endl;

	RnaAtomLib atLib;
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	for(int i=0;i<this->seqLen-1;i++){
		vector<Atom*> aList = nodes[i]->phoAtoms();
		RNABase* base = rc.getBaseList()[i+1];
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
	}

	rc.printPDBFormat(of, 1);
	of << "ENDMDL" << endl;

}

void BRTreeInfo::printPDB(const string& outputFile) {
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}

	for(int i=0;i<this->seqLen-1;i++){
		vector<Atom*> aList = nodes[i]->phoAtoms();
		RNABase* base = rc.getBaseList()[i+1];
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
	}

	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of << "ene: " << this->ene << endl;
	of.close();


}

void BRTreeInfo::printTmpPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {
		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), "A", s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toTmpAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
		atomNum += aList.size();
	}
	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of.close();
}

BRTreeInfo::~BRTreeInfo() {
	delete [] seq;
	delete [] connectToDownstream;
	for(int i=0;i<seqLen;i++){
		delete nodes[i];
	}
	delete [] nodes;
}

BRFoldingTree::BRFoldingTree(const string& inputFile, EnergyTable* et){

	this->et = et;

	this->rotLib = new BaseRotamerLib(); //deleted
	cout << "rot lib" << endl;
	this->fragLib = new FragmentLibrary(rotLib);  //deleted

	cout << "frag lib" << endl;

	RNABaseLib baseLib;

	NSPtools::InputParser input(inputFile);

	input.printOptions();

	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string nwcSec = input.getValue("nwc");

	string flexOption = input.getValue("fixed");
	string fixedGroup = input.getValue("group");
	string chainBreak = input.getValue("break");
	string bulge = input.getValue("bulge");
	string ag = input.getValue("AG");

	vector<string> spt;
	vector<string> spt2;
	vector<int> chainBreakPoints;


	bool freeIndexTag[baseSeq.length()];
	bool riboFreeIndexTag[baseSeq.length()];

	for(int i=0;i<baseSeq.length();i++){
		freeIndexTag[i] = true;
		riboFreeIndexTag[i] = true;
	}

	splitString(flexOption, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		fixedList.push_back(atoi(spt[i].c_str()));
		freeIndexTag[atoi(spt[i].c_str())] = false;
	}

	splitString(fixedGroup, " ", &spt);
	int* groupIndex = new int[baseSeq.length()];
	for(int i=0;i<baseSeq.length();i++)
		groupIndex[i] = -1;

	for(unsigned int i=0;i<spt.size();i++){
		splitString(spt[i], ",", &spt2);
		set<int> group;
		for(unsigned int j=0;j<spt2.size();j++) {
			group.insert(atoi(spt2[j].c_str()));
			groupIndex[atoi(spt2[j].c_str())] = i;
			freeIndexTag[atoi(spt2[j].c_str())] = false;
		}
		this->fixedGroups.push_back(group);
	}


	splitString(chainBreak, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		int pos = atoi(spt[i].c_str());
		if(pos >= baseSeq.length()){
			cout << "invalid chain break position: " << pos << endl;
			exit(0);
		}
		chainBreakPoints.push_back(atoi(spt[i].c_str()));
	}

	splitString(bulge, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		bulgeList.push_back(atoi(spt[i].c_str()));
	}

	splitString(ag, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		agList.push_back(atoi(spt[i].c_str()));
	}

	RNAPDB pdb(pdbFile, "xxxx");
	vector<RNABase*> baseList0 = pdb.getBaseList();
	vector<RNABase*> baseList;
	for(int i=0;i<baseList0.size();i++){
		baseList.push_back(baseLib.toStandardBase(baseList0[i]));
	}

	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size: " << baseList.size() << " not equal to sequence length: " << baseSeq.length() << endl;
		exit(0);
	}
	if(baseSeq.length() != baseSec.length()) {
		cout << "seq length not equal to sec length" << endl;
		exit(0);
	}

	this->seqLen = baseSeq.length();


	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->nwcPairPosID = new int[seqLen];

	this->fixed = new bool[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->loopRevPoints = new bool[seqLen];
	this->chainBreakPoints = new bool[seqLen];

	this->nodes = new BRNode*[seqLen];
	this->initRotList = new BaseRotamer*[seqLen];

	this->sepTable = new int[seqLen*seqLen];

	this->allBaseClashE = new double[seqLen*seqLen];
	this->tmpBaseClashE = new double[seqLen*seqLen];
	this->allBaseBaseE = new double[seqLen*seqLen];
	this->tmpBaseBaseE = new double[seqLen*seqLen];
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

	for(int i=0;i<seqLen;i++){
		allRotE[i] = 0;
		tmpRotE[i] = 0;
		allRcE[i] = 0;
		tmpRcE[i] = 0;
		for(int j=0;j<seqLen;j++){
			int pi = i*seqLen+j;
			allBaseClashE[pi] = 0;
			tmpBaseClashE[pi] = 0;
			allBaseBaseE[pi] = 0;
			tmpBaseBaseE[pi] = 0;
			allBaseRiboseE[pi] = 0;
			tmpBaseRiboseE[pi] = 0;
			allBasePhoE[pi] = 0;
			tmpBasePhoE[pi] = 0;
			allRiboseRiboseE[pi] = 0;
			tmpRiboseRiboseE[pi] = 0;
			allRibosePhoE[pi] = 0;
			tmpRibosePhoE[pi] = 0;
			allPhoPhoE[pi] = 0;
			tmpPhoPhoE[pi] = 0;
		}
	}

	RNABaseName rn;
	for(unsigned int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;
		this->nwcPairPosID[i] = -1;
		this->fixed[i] = false;

		BaseRotamer* rot;
		if(!baseList[i]->backboneComplete())
			rot = rotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		else
			rot = new BaseRotamer(baseList[i]);

		this->initRotList[i] = rot;
		this->nodes[i] = new BRNode(baseList[i], rot);
		this->nodes[i]->groupID = groupIndex[i];
		cout << "node: " << i << " " << nodes[i]->baseType << endl;

		/*
		if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1])) {
			this->connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
		else {
			this->connectToDownstream[i] = false;
			nodes[i]->connectToNeighbor = false;
		}
		*/

		LocalFrame cs = baseList[i]->getCoordSystem();
		this->nodes[i]->tmpRot = this->nodes[i]->rot;

	}
	delete [] groupIndex;


	for(int i=0;i<seqLen-1;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;
	for(int i=0;i<chainBreakPoints.size();i++){
		connectToDownstream[chainBreakPoints[i]] = false;
	}
	for(int i=0;i<seqLen;i++){
		nodes[i]->connectToNeighbor = connectToDownstream[i];
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

	for(int i=0;i<fixedList.size();i++) {
		this->fixed[fixedList[i]] = true;
		this->nodes[fixedList[i]]->fixed = true;
	}

	for(int i=0;i<seqLen;i++){
		if(!fixed[i])
			flexibleNodes.push_back(nodes[i]);
	}

	this->initTreeInfo = new BRTreeInfo(seqLen, seq, connectToDownstream, nodes, this->flexibleNodes, 0.0);

	/*
	 * parse secondary structure
	 */

	char ss[seqLen];
	for(int i=0;i<seqLen;i++){
		ss[i] = baseSec[i];
	}

	for(int i=0;i<seqLen;i++) {
		char c = ss[i];
		if(c == ')') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == '>') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '<') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
	}

	/*
	 * parse non-Watson-Crick pair
	 */

	char nss[seqLen];
	for(int i=0;i<seqLen;i++){
		nss[i] = nwcSec[i];
	}

	for(int i=0;i<seqLen;i++) {
		char c = nss[i];
		if(c == ')') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == '>') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '<') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
	}

	if(fixedList.size() == 0)
		this->rootNode = nodes[0];
	else
		this->rootNode = nodes[fixedList[0]];

	pseudoNode = new BRNode(0,-1);
	this->rootNode->father = pseudoNode;


	for(int i=0;i<seqLen;i++){
		printf("%-2d ", wcPairPosID[i]);
	}
	cout << endl;
	for(int i=0;i<seqLen;i++){
		printf("%-2d ", nwcPairPosID[i]);
	}

	cout << endl;

	this->nodeConnectMatrix = new bool[seqLen*seqLen];
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			if(i == j) this->nodeConnectMatrix[i*seqLen+j] = true;
			else this->nodeConnectMatrix[i*seqLen+j] = false;
		}
	}

	cout << "FIX" << endl;
	buildFixedNodes();
	cout << "FIX pair" << endl;
	buildBasePairOfFixedNodes();
	cout << "AG pair" << endl;
	buildAGPairs();
	cout << "PAIR" << endl;
	buildPairs();
	cout << "BULGE" << endl;
	buildBulged13();
	cout << "MID" << endl;
	findMidPoint();
	cout << "NEIGHBOR" << endl;
	buildNeighbor();
	cout << "REVERSE" << endl;
	buildReverseNodes();

	for(int i=0;i<baseSeq.length();i++){
		riboFreeIndexTag[i] = freeIndexTag[i];
	}

	for(int i=0;i<seqLen;i++){
		if(freeIndexTag[i] && i>0 && connectToDownstream[i-1])
			riboFreeIndexTag[i-1] = true;
		if(freeIndexTag[i] && i<seqLen-1 && connectToDownstream[i])
			riboFreeIndexTag[i+1] = true;
	}

	for(int i=0;i<flexibleConnectionList.size();i++){
		if(flexibleConnectionList[i]->ctType == "wc" || flexibleConnectionList[i]->ctType == "nwc"){
			BRNode* cn = flexibleConnectionList[i]->childNode;
			freeIndexTag[cn->seqID] = false;
		}
	}


	for(int i=0;i<baseSeq.length();i++){
		if(freeIndexTag[i])
		{
			this->freeNodeIDs.push_back(i);
		}
	}

	for(int i=0;i<seqLen;i++){
		if(riboFreeIndexTag[i])
			riboFlexibleNodes.push_back(nodes[i]);
	}

	for(int i=0;i<fixedConnectionList.size();i++) {
		cout << "ct: " << fixedConnectionList[i]->fatherNode->seqID << " " << fixedConnectionList[i]->childNode->seqID << endl;
		fixedConnectionList[i]->childNode->upConnection = fixedConnectionList[i];
	}

	for(int i=0;i<flexibleConnectionList.size();i++) {
		cout << "ct: " << flexibleConnectionList[i]->fatherNode->seqID << " " << flexibleConnectionList[i]->childNode->seqID << endl;
		flexibleConnectionList[i]->childNode->upConnection = flexibleConnectionList[i];
	}

	for(int i=0;i<fixedConnectionList.size();i++) {

		cout << "connect: " << fixedConnectionList[i]->fatherNode->seqID << " " << fixedConnectionList[i]->childNode->seqID << endl;
		fixedConnectionList[i]->updateChildInfo(nodes);
	}

	for(int i=0;i<flexibleConnectionList.size();i++) {
		flexibleConnectionList[i]->updateChildInfo(nodes);
	}

	for(int i=0;i<seqLen;i++){
		cout << "update node child info: " << i << endl;
		nodes[i]->updateChildInfo(nodes, seqLen);
	}

	cout << "check ribose" << endl;
	checkRibose();

	updatePhoGroups();

	cout << "check ribose" << endl;
	checkRibose();

	updateEnergies(0.0);
}

BRFoldingTree::BRFoldingTree(const string& inputFile){
	this->rotLib = new BaseRotamerLib(); //deleted
	cout << "rot lib" << endl;
	this->fragLib = new FragmentLibrary(rotLib);  //deleted

	cout << "frag lib" << endl;

	RNABaseLib baseLib;

	NSPtools::InputParser input(inputFile);

	input.printOptions();

	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string nwcSec = input.getValue("nwc");

	string flexOption = input.getValue("fixed");
	string fixedGroup = input.getValue("group");
	string chainBreak = input.getValue("break");
	string bulge = input.getValue("bulge");

	vector<string> spt;
	vector<string> spt2;
	vector<int> chainBreakPoints;


	bool freeIndexTag[baseSeq.length()];
	bool riboFreeIndexTag[baseSeq.length()];

	for(int i=0;i<baseSeq.length();i++){
		freeIndexTag[i] = true;
		riboFreeIndexTag[i] = true;
	}


	splitString(flexOption, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		fixedList.push_back(atoi(spt[i].c_str()));
		freeIndexTag[atoi(spt[i].c_str())] = false;
	}

	splitString(chainBreak, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		int pos = atoi(spt[i].c_str());
		if(pos >= baseSeq.length()){
			cout << "invalid chain break position: " << pos << endl;
			exit(0);
		}
		chainBreakPoints.push_back(atoi(spt[i].c_str()));
	}

	splitString(fixedGroup, " ", &spt);
	int* groupIndex = new int[baseSeq.length()];
	for(int i=0;i<baseSeq.length();i++)
		groupIndex[i] = -1;

	for(unsigned int i=0;i<spt.size();i++){
		splitString(spt[i], ",", &spt2);
		set<int> group;
		for(unsigned int j=0;j<spt2.size();j++) {
			group.insert(atoi(spt2[j].c_str()));
			groupIndex[atoi(spt2[j].c_str())] = i;
			freeIndexTag[atoi(spt2[j].c_str())] = false;
		}
		this->fixedGroups.push_back(group);
	}

	splitString(bulge, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		bulgeList.push_back(atoi(spt[i].c_str()));
	}



	RNAPDB pdb(pdbFile, "xxxx");
	vector<RNABase*> baseList0 = pdb.getBaseList();
	vector<RNABase*> baseList;
	for(int i=0;i<baseList0.size();i++){
		baseList.push_back(baseLib.toStandardBase(baseList0[i]));
	}

	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size: " << baseList.size() << " not equal to sequence length: " << baseSeq.length() << endl;
		exit(0);
	}
	if(baseSeq.length() != baseSec.length()) {
		cout << "seq length not equal to sec length" << endl;
		exit(0);
	}

	this->seqLen = baseSeq.length();


	this->seq = new int[seqLen];
	this->wcPairPosID = new int[seqLen];
	this->nwcPairPosID = new int[seqLen];

	this->fixed = new bool[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->loopRevPoints = new bool[seqLen];
	this->chainBreakPoints = new bool[seqLen];

	this->nodes = new BRNode*[seqLen];
	this->initRotList = new BaseRotamer*[seqLen];

	this->sepTable = new int[seqLen*seqLen];

	this->allBaseClashE = new double[seqLen*seqLen];
	this->tmpBaseClashE = new double[seqLen*seqLen];
	this->allBaseBaseE = new double[seqLen*seqLen];
	this->tmpBaseBaseE = new double[seqLen*seqLen];
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

	for(int i=0;i<seqLen;i++){
		allRotE[i] = 0;
		tmpRotE[i] = 0;
		allRcE[i] = 0;
		tmpRcE[i] = 0;
		for(int j=0;j<seqLen;j++){
			int pi = i*seqLen+j;
			allBaseClashE[pi] = 0;
			tmpBaseClashE[pi] = 0;
			allBaseBaseE[pi] = 0;
			tmpBaseBaseE[pi] = 0;
			allBaseRiboseE[pi] = 0;
			tmpBaseRiboseE[pi] = 0;
			allBasePhoE[pi] = 0;
			tmpBasePhoE[pi] = 0;
			allRiboseRiboseE[pi] = 0;
			tmpRiboseRiboseE[pi] = 0;
			allRibosePhoE[pi] = 0;
			tmpRibosePhoE[pi] = 0;
			allPhoPhoE[pi] = 0;
			tmpPhoPhoE[pi] = 0;
		}
	}

	RNABaseName rn;
	for(unsigned int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;
		this->nwcPairPosID[i] = -1;
		this->fixed[i] = false;

		BaseRotamer* rot;
		if(!baseList[i]->backboneComplete())
		{
			rot = rotLib->getLowestEnergyRotamer(baseList[i]->baseTypeInt);
		}
		else {
			rot = new BaseRotamer(baseList[i]);
		}
		this->initRotList[i] = rot;
		this->nodes[i] = new BRNode(baseList[i], rot);
		this->nodes[i]->groupID = groupIndex[i];
		cout << "node: " << i << " " << nodes[i]->baseType << endl;

		/*
		if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1])) {
			this->connectToDownstream[i] = true;
			nodes[i]->connectToNeighbor = true;
		}
		else {
			this->connectToDownstream[i] = false;
			nodes[i]->connectToNeighbor = false;
		}
		*/

		LocalFrame cs = baseList[i]->getCoordSystem();
		this->nodes[i]->tmpRot = this->nodes[i]->rot;

	}
	delete [] groupIndex;


	for(int i=0;i<seqLen-1;i++){
		connectToDownstream[i] = true;
	}
	connectToDownstream[seqLen-1] = false;
	for(int i=0;i<chainBreakPoints.size();i++){
		connectToDownstream[chainBreakPoints[i]] = false;
	}
	for(int i=0;i<seqLen;i++){
		nodes[i]->connectToNeighbor = connectToDownstream[i];
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

	for(int i=0;i<fixedList.size();i++) {
		this->fixed[fixedList[i]] = true;
		this->nodes[fixedList[i]]->fixed = true;
	}

	for(int i=0;i<seqLen;i++){
		if(!fixed[i])
			flexibleNodes.push_back(nodes[i]);
	}

	this->initTreeInfo = new BRTreeInfo(seqLen, seq, connectToDownstream, nodes, this->flexibleNodes, 0.0);

	/*
	 * parse secondary structure
	 */

	char ss[seqLen];
	for(int i=0;i<seqLen;i++){
		ss[i] = baseSec[i];
	}

	for(int i=0;i<seqLen;i++) {
		char c = ss[i];
		if(c == ')') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == '>') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(ss[j] == '<') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << baseSec << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
	}

	/*
	 * parse non-Watson-Crick pair
	 */

	char nss[seqLen];
	for(int i=0;i<seqLen;i++){
		nss[i] = nwcSec[i];
	}

	for(int i=0;i<seqLen;i++) {
		char c = nss[i];
		if(c == ')') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
		if(c == '>') {
			int preIndex = -1;
			for(int j=i-1;j>=0;j--) {
				if(nss[j] == '<') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid nwcSeq: " << nwcSec << endl;
				exit(1);
			}
			nss[i] = '.';
			nss[preIndex] = '.';
			nwcPairPosID[i] = preIndex;
			nwcPairPosID[preIndex] = i;
		}
	}

	if(fixedList.size() == 0)
		this->rootNode = nodes[0];
	else
		this->rootNode = nodes[fixedList[0]];

	pseudoNode = new BRNode(0,-1);
	this->rootNode->father = pseudoNode;


	for(int i=0;i<seqLen;i++){
		printf("%-2d ", wcPairPosID[i]);
	}
	cout << endl;
	for(int i=0;i<seqLen;i++){
		printf("%-2d ", nwcPairPosID[i]);
	}

	cout << endl;



	this->nodeConnectMatrix = new bool[seqLen*seqLen];
	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			if(i == j) this->nodeConnectMatrix[i*seqLen+j] = true;
			else this->nodeConnectMatrix[i*seqLen+j] = false;
		}
	}

	printConnections();
	cout << "FIX" << endl;
	buildFixedNodes();
	cout << "FIX pair" << endl;
	buildBasePairOfFixedNodes();
	cout << "PAIR" << endl;
	buildPairs();
	cout << "BULGE" << endl;
	buildBulged13();
	cout << "MID" << endl;
	findMidPoint();
	cout << "NEIGHBOR" << endl;
	buildNeighbor();
	cout << "REVNEIGHBOR" << endl;
	buildReverseNodes();
}

void BRFoldingTree::updateCoordinate(BRNode* node){
	BRConnection* ct = node->upConnection;
	if(ct != NULL && !node->fixed){

		//cout << "update coordinate: " << node->seqID << endl;
		BRNode* father = ct->fatherNode;
		node->cs1 = father->cs1 + ct->cm;
		node->tmpCs1 = node->cs1;
		node->tmpRot = node->rot;
		for(int i=0;i<node->baseAtomNum;i++){
			node->baseAtomCoords[i] = local2global(node->cs1, node->atomCoordLocal[i]);
			node->baseAtomCoordsTmp[i] = node->baseAtomCoords[i];
		}
		node->cs2 = node->cs1 + node->rot->mv12;
		node->cs3 = node->cs1 + node->rot->mv13;
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(int i=0;i<11;i++){
			node->riboAtomCoords[i] = local2global(node->cs1, node->rot->tList1[i]);
			node->riboAtomCoordsTmp[i] = node->riboAtomCoords[i];
		}
	}

	if(node->leftChild != NULL)
		updateCoordinate(node->leftChild);
	if(node->midChild != NULL)
		updateCoordinate(node->midChild);
	if(node->rightChild != NULL)
		updateCoordinate(node->rightChild);
	if(node->reverseChild != NULL)
		updateCoordinate(node->reverseChild);

}

void BRFoldingTree::randInit(){

	for(int i=0;i<riboFlexibleNodes.size();i++){
		cout << "random init node: " << riboFlexibleNodes[i]->seqID << endl;
		riboFlexibleNodes[i]->rot = rotLib->getRandomRotamerLv1(riboFlexibleNodes[i]->baseType);
		riboFlexibleNodes[i]->tmpRot = riboFlexibleNodes[i]->rot;
		riboFlexibleNodes[i]->cs2 = riboFlexibleNodes[i]->cs1 + riboFlexibleNodes[i]->rot->mv12;
		riboFlexibleNodes[i]->tmpCs2 = riboFlexibleNodes[i]->cs2;
		riboFlexibleNodes[i]->cs3 = riboFlexibleNodes[i]->cs1 + riboFlexibleNodes[i]->rot->mv13;
		riboFlexibleNodes[i]->tmpCs3 = riboFlexibleNodes[i]->cs3;


		for(int j=0;j<11;j++){
			riboFlexibleNodes[i]->riboAtomCoords[j] = local2global(riboFlexibleNodes[i]->cs1, riboFlexibleNodes[i]->rot->tList1[j]);
			riboFlexibleNodes[i]->riboAtomCoordsTmp[j] = local2global(riboFlexibleNodes[i]->tmpCs1, riboFlexibleNodes[i]->tmpRot->tList1[j]);
		}
	}

	for(int i=0;i<flexibleConnectionList.size();i++){
		cout << "random init connection: " << flexibleConnectionList[i]->fatherNode->seqID << " " << flexibleConnectionList[i]->childNode->seqID << endl;
		F2Fragment* frag = flexibleConnectionList[i]->f2Lib->getRandomFrag();
		flexibleConnectionList[i]->cm = frag->cm;
		flexibleConnectionList[i]->cmTmp = flexibleConnectionList[i]->cm;
		flexibleConnectionList[i]->f2Frag = frag;
		flexibleConnectionList[i]->f2FragTmp = frag;
	}

	updateCoordinate(rootNode);
	updatePhoGroups();
	updateEnergies(0.0);
}

void BRFoldingTree::buildFixedNodes() {
	if(this->fixedList.size() == 0)
		return;
	for(int i=1;i<fixedList.size();i++){
		cout << "build fixed connect: " << fixedList[i-1] << "->" << fixedList[i] << endl;

		nodes[fixedList[i-1]]->rightChild = nodes[fixedList[i]];
		nodes[fixedList[i]]->father = nodes[fixedList[i-1]];

		int idA = fixedList[i-1];
		int idB = fixedList[i];

		vector<int> ctList;
		for(int a=0;a<seqLen;a++){
			if(nodeConnectMatrix[a*seqLen+idA]){
				ctList.push_back(a);
			}
			else if(nodeConnectMatrix[a*seqLen+idB]){
				ctList.push_back(a);
			}
		}
		for(int a=0;a<ctList.size();a++){
			for(int b=0;b<ctList.size();b++){
				nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
			}
		}
		//printBaseConnection();

		BRConnection* ct = new BRConnection(fragLib, "jump" ,nodes[fixedList[i-1]], nodes[fixedList[i]], seqLen);
		ct->fixed = true;
		this->fixedConnectionList.push_back(ct);
	}
}

void BRFoldingTree::buildGroupNodes(){
	int n = this->fixedGroups.size();
	int pre, cur;
	BRNode* nodeA;
	BRNode* nodeB;
	set<int>::iterator it;

	for(int i=0;i<n;i++){

		it = fixedGroups[i].begin();
		int pre = *it;
		it++;
		while(it != fixedGroups[i].end()){
			cur = *it;
			//build ct

			nodeA = nodes[pre];
			nodeB = nodes[cur];

			cout << "build group connect: " << pre << "->" << cur << endl;

			nodeA->rightChild = nodeB;
			nodeB->father = nodeA;

			int idA = nodeA->seqID;
			int idB = nodeB->seqID;
			vector<int> ctList;
			for(int a=0;a<seqLen;a++){
				if(nodeConnectMatrix[a*seqLen+idA]){
					ctList.push_back(a);
				}
				else if(nodeConnectMatrix[a*seqLen+idB]){
					ctList.push_back(a);
				}
			}
			for(int a=0;a<ctList.size();a++){
				for(int b=0;b<ctList.size();b++){
					nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
				}
			}
			//printBaseConnection();
			BRConnection* ct = new BRConnection(fragLib, "jump" ,nodeA, nodeB, seqLen);
			ct->fixed = true;
			this->fixedConnectionList.push_back(ct);

			pre = cur;
			it++;
		}
	}
}

void BRFoldingTree::buildBasePairOfFixedNodes() {
	if(this->fixedList.size() == 0)
		return;
	for(int x=0;x<fixedList.size();x++){
		BRNode* node = this->nodes[fixedList[x]];
		int seqID = node->seqID;
		if(wcPairPosID[seqID] > 0){
			BRNode* wc = nodes[wcPairPosID[seqID]];
			if(wc->father == NULL){
				node->leftChild = wc;
				wc->father = node;
				BRConnection* ct = new BRConnection(fragLib, "wc", node, wc, seqLen);
				cout << "build wc connection: " << node->seqID << " " << wc->seqID << endl;
				int idA = node->seqID;
				int idB = wc->seqID;
				vector<int> ctList;
				for(int a=0;a<seqLen;a++){
					if(nodeConnectMatrix[a*seqLen+idA]){
						ctList.push_back(a);
					}
					else if(nodeConnectMatrix[a*seqLen+idB]){
						ctList.push_back(a);
					}
				}
				for(int a=0;a<ctList.size();a++){
					for(int b=0;b<ctList.size();b++){
						nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
					}
				}
				//printBaseConnection();
				bool isFixed = false;
				for(int i=0;i<fixedGroups.size();i++) {
					if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
						isFixed = true;
				}
				if(isFixed || fixed[idB]) {
					ct->fixed = true;
					this->fixedConnectionList.push_back(ct);
				}
				else {
					this->flexibleConnectionList.push_back(ct);
				}
			}
		}
		else if(nwcPairPosID[seqID] > 0) {
			BRNode* nwc = nodes[nwcPairPosID[seqID]];
			if(nwc->father == NULL){
				node->leftChild = nwc;
				nwc->father = node;
				BRConnection* ct = new BRConnection(fragLib, "nwc", node, nwc, seqLen);
				cout << "build nwc connection: " << node->seqID << " " << nwc->seqID << endl;
				bool isFixed = false;
				int idA = node->seqID;
				int idB = nwcPairPosID[seqID];
				vector<int> ctList;
				for(int a=0;a<seqLen;a++){
					if(nodeConnectMatrix[a*seqLen+idA]){
						ctList.push_back(a);
					}
					else if(nodeConnectMatrix[a*seqLen+idB]){
						ctList.push_back(a);
					}
				}
				for(int a=0;a<ctList.size();a++){
					for(int b=0;b<ctList.size();b++){
						nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
					}
				}
				//printBaseConnection();
				for(int i=0;i<fixedGroups.size();i++) {

					if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
						isFixed = true;
				}
				if(isFixed || fixed[idB]) {
					ct->fixed = true;
					this->fixedConnectionList.push_back(ct);
				}
				else
					this->flexibleConnectionList.push_back(ct);
			}
		}
	}
}

void BRFoldingTree::buildAGPairs(){
	BRNode* nodeA;
	BRNode* nodeB;
	for(int x=0;x<this->agList.size();x++){
		int k = agList[x];
		if(nwcPairPosID[k] >= 0 && nodes[k]->leftChild == NULL && nodes[nwcPairPosID[k]]->father == NULL){
			nodeA = nodes[k];
			nodeB = nodes[nwcPairPosID[k]];

			nodeA->leftChild = nodeB;
			nodeB->father = nodeA;

			string ctType = "AG";
			if(nodeA->baseType == 0 && nodeB->baseType == 2)
				ctType = "AG";
			else
				ctType = "GA";

			BRConnection* ct = new BRConnection(fragLib, ctType, nodeA, nodeB, seqLen);
			cout << "build " + ctType + " connection: " << nodeA->seqID << " " << nodeB->seqID << endl;
			bool isFixed = false;
			int idA = nodeA->seqID;
			int idB = nodeB->seqID;
			vector<int> ctList;
			for(int a=0;a<seqLen;a++){
				if(nodeConnectMatrix[a*seqLen+idA]){
					ctList.push_back(a);
				}
				else if(nodeConnectMatrix[a*seqLen+idB]){
					ctList.push_back(a);
				}
			}
			for(int a=0;a<ctList.size();a++){
				for(int b=0;b<ctList.size();b++){
					nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
				}
			}
			//printBaseConnection();
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else {
				this->flexibleConnectionList.push_back(ct);
			}
		}
	}
}

void BRFoldingTree::buildPairs(){

	BRNode* nodeA;
	BRNode* nodeB;
	for(int k=0;k<this->seqLen;k++){
		if(wcPairPosID[k] >=0 ){
			if(this->nodes[k]->father != NULL || this->nodes[wcPairPosID[k]]->father != NULL) continue;

			if(this->nodes[k]->leftChild == NULL){
				nodeA = this->nodes[k];
				nodeB = this->nodes[wcPairPosID[k]];
			}
			else
				continue;

			nodeA->leftChild = nodeB;
			nodeB->father = nodeA;

			BRConnection* ct = new BRConnection(fragLib, "wc", nodeA, nodeB, seqLen);
			cout << "build wc connection: " << nodeA->seqID << " " << nodeB->seqID << endl;
			bool isFixed = false;
			int idA = nodeA->seqID;
			int idB = nodeB->seqID;
			vector<int> ctList;
			for(int a=0;a<seqLen;a++){
				if(nodeConnectMatrix[a*seqLen+idA]){
					ctList.push_back(a);
				}
				else if(nodeConnectMatrix[a*seqLen+idB]){
					ctList.push_back(a);
				}
			}
			for(int a=0;a<ctList.size();a++){
				for(int b=0;b<ctList.size();b++){
					nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
				}
			}
			//printBaseConnection();
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else {
				this->flexibleConnectionList.push_back(ct);
			}
		}
	}

	for(int k=0;k<this->seqLen;k++){
		if(nwcPairPosID[k] >=0 ){
			if(this->nodes[k]->father != NULL || this->nodes[nwcPairPosID[k]]->father != NULL) continue;
			if(this->nodes[k]->leftChild == NULL){
				nodeA = this->nodes[k];
				nodeB = this->nodes[nwcPairPosID[k]];
			}
			else
				continue;

			nodeA->leftChild = nodeB;
			nodeB->father = nodeA;

			BRConnection* ct = new BRConnection(fragLib, "nwc", nodeA, nodeB, seqLen);
			cout << "build nwc connection: " << nodeA->seqID << " " << nodeB->seqID << endl;
			bool isFixed = false;
			int idA = nodeA->seqID;
			int idB = nodeB->seqID;
			vector<int> ctList;
			for(int a=0;a<seqLen;a++){
				if(nodeConnectMatrix[a*seqLen+idA]){
					ctList.push_back(a);
				}
				else if(nodeConnectMatrix[a*seqLen+idB]){
					ctList.push_back(a);
				}
			}
			for(int a=0;a<ctList.size();a++){
				for(int b=0;b<ctList.size();b++){
					nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
				}
			}
			//printBaseConnection();
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else {
				this->flexibleConnectionList.push_back(ct);
			}
		}
	}
}



void BRFoldingTree::buildBulged13(){
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

	for(int i=1;i<seqLen-1;i++){
		if(pairID[i-1] >=0 && pairID[i] < 0 && pairID[i+1] >=0 && pairID[i-1] == pairID[i+1]+1 && rand()%2 ==0){
			bulgeList.push_back(i);
		}
	}

	for(int k=0;k<bulgeList.size();k++){
		BRNode* node = nodes[bulgeList[k] -1];
		BRNode* nnb = nodes[node->seqID + 2];
		if(nnb->father != NULL) continue;
		if(nodeConnectMatrix[node->seqID*seqLen+nnb->seqID]) continue;
		BRConnection* ct = new BRConnection(fragLib, "bulge13", node, nnb, seqLen);
		node->bulge13Child = nnb;
		nnb->father = node;
		cout << "build bulge13 connection: " << node->seqID << " " << nnb->seqID << endl;
		bool isFixed = false;
		int idA = node->seqID;
		int idB = idA + 2;
		vector<int> ctList;
		for(int a=0;a<seqLen;a++){
			if(nodeConnectMatrix[a*seqLen+idA]){
				ctList.push_back(a);
			}
			else if(nodeConnectMatrix[a*seqLen+idB]){
				ctList.push_back(a);
			}
		}
		for(int a=0;a<ctList.size();a++){
			for(int b=0;b<ctList.size();b++){
				nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
			}
		}
		//printBaseConnection();
		for(int i=0;i<fixedGroups.size();i++) {

			if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
				isFixed = true;
		}
		if(isFixed || fixed[idB]) {
			ct->fixed = true;
			this->fixedConnectionList.push_back(ct);
		}
		else
			this->flexibleConnectionList.push_back(ct);
	}
}

void BRFoldingTree::buildNeighbor(){

	for(int k=0;k<seqLen-1;k++){
		BRNode* node = nodes[k];
		BRNode* nb = nodes[k+1];
		if(nodeConnectMatrix[node->seqID*seqLen + nb->seqID]) continue;
		if(connectToDownstream[k] && nb->father == NULL) {


			if(loopRevPoints[k] == true) {
				cout << "stop at loop mid: " << k << endl;
				//do nothing
			}
			// helix neighbor pair
			else if(wcPairPosID[node->seqID] > 0 && wcPairPosID[k+1] > 0 &&
			   wcPairPosID[node->seqID] == wcPairPosID[k+1] + 1) {

				node->midChild = nb;
				nb->father = node;
				BRConnection* ct = new BRConnection(fragLib, "wcNb" ,node, nb, seqLen);

				cout << "build wc nb connection: " << node->seqID << " " << nb->seqID << endl;
				bool isFixed = false;
				int idA = node->seqID;
				int idB = k + 1;
				vector<int> ctList;
				for(int a=0;a<seqLen;a++){
					if(nodeConnectMatrix[a*seqLen+idA]){
						ctList.push_back(a);
					}
					else if(nodeConnectMatrix[a*seqLen+idB]){
						ctList.push_back(a);
					}
				}
				for(int a=0;a<ctList.size();a++){
					for(int b=0;b<ctList.size();b++){
						nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
					}
				}
			//	printBaseConnection();
				for(int i=0;i<fixedGroups.size();i++) {

					if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
						isFixed = true;
				}
				if(isFixed || fixed[idB]) {
					ct->fixed = true;
					this->fixedConnectionList.push_back(ct);
				}
				else
					this->flexibleConnectionList.push_back(ct);
			}
			// non-watson crick helix neighbor
			else if((wcPairPosID[node->seqID] > 0 && nwcPairPosID[k+1] > 0  && wcPairPosID[node->seqID] == nwcPairPosID[k+1] + 1) ||
					(nwcPairPosID[node->seqID] > 0 && wcPairPosID[k+1] > 0  && nwcPairPosID[node->seqID] == wcPairPosID[k+1] + 1) ||
					(nwcPairPosID[node->seqID] > 0 && nwcPairPosID[k+1] > 0  && nwcPairPosID[node->seqID] == nwcPairPosID[k+1] + 1)) {
				BRConnection* ct = new BRConnection(fragLib, "nwcNb", node, nb, seqLen);
				node->midChild = nb;
				nb->father = node;
				cout << "build nwc nb connection: " << node->seqID << " " << nb->seqID << endl;
				bool isFixed = false;
				int idA = node->seqID;
				int idB = k + 1;
				vector<int> ctList;
				for(int a=0;a<seqLen;a++){
					if(nodeConnectMatrix[a*seqLen+idA]){
						ctList.push_back(a);
					}
					else if(nodeConnectMatrix[a*seqLen+idB]){
						ctList.push_back(a);
					}
				}
				for(int a=0;a<ctList.size();a++){
					for(int b=0;b<ctList.size();b++){
						nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
					}
				}
				//printBaseConnection();

				for(int i=0;i<fixedGroups.size();i++) {

					if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
						isFixed = true;
				}
				if(isFixed || fixed[idB]) {
					ct->fixed = true;
					this->fixedConnectionList.push_back(ct);
				}
				else
					this->flexibleConnectionList.push_back(ct);
			}
			// loop neighbor
			else {
				BRConnection* ct = new BRConnection(fragLib, "loopNb", node, nb, seqLen);
				node->midChild = nb;
				nb->father = node;
				cout << "build loop nb connection: " << node->seqID << " " << nb->seqID << endl;
				bool isFixed = false;
				int idA = node->seqID;
				int idB = k + 1;
				vector<int> ctList;
				for(int a=0;a<seqLen;a++){
					if(nodeConnectMatrix[a*seqLen+idA]){
						ctList.push_back(a);
					}
					else if(nodeConnectMatrix[a*seqLen+idB]){
						ctList.push_back(a);
					}
				}
				for(int a=0;a<ctList.size();a++){
					for(int b=0;b<ctList.size();b++){
						nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
					}
				}
				//printBaseConnection();
				for(int i=0;i<fixedGroups.size();i++) {
					if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
						isFixed = true;
				}
				if(isFixed || fixed[idB]) {
					ct->fixed = true;
					this->fixedConnectionList.push_back(ct);
				}
				else
					this->flexibleConnectionList.push_back(ct);
			}
		}
	}
}


void BRFoldingTree::buildFrom2(BRNode* node) {

	/*
	 * base first, ribose second
	 */

	cout << "build node: " << node->seqID << endl;
	//if(node->father != NULL && node->seqID < seqLen-1)
	//	buildFrom2(nodes[node->seqID+1]);

	/*
	 * build left child
	 */
	int seqID = node->seqID;

	if(wcPairPosID[seqID] > 0){
		BRNode* wc = nodes[wcPairPosID[seqID]];
		if(wc->father == NULL){
			node->leftChild = wc;
			wc->father = node;
			BRConnection* ct = new BRConnection(fragLib, "wc", node, wc, seqLen);
			cout << "build wc connection: " << node->seqID << " " << wc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = wcPairPosID[seqID];
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else {
				this->flexibleConnectionList.push_back(ct);
			}
		}
	}
	else if(nwcPairPosID[seqID] > 0) {
		BRNode* nwc = nodes[nwcPairPosID[seqID]];
		if(nwc->father == NULL){
			node->leftChild = nwc;
			nwc->father = node;
			BRConnection* ct = new BRConnection(fragLib, "nwc", node, nwc, seqLen);
			cout << "build nwc connection: " << node->seqID << " " << nwc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nwcPairPosID[seqID];
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
	}

	/*
	 * build middle node
	 */
	if(seqID + 1 < this->seqLen && connectToDownstream[seqID] && nodes[seqID+1]->father == NULL) {
		BRNode* nb = nodes[seqID+1];
		// helix neighbor pair
		if(wcPairPosID[node->seqID] > 0 && wcPairPosID[seqID+1] > 0 &&
		   wcPairPosID[node->seqID] == wcPairPosID[seqID+1] + 1) {
			node->midChild = nb;
			nb->father = node;
			BRConnection* ct = new BRConnection(fragLib, "wcNb" ,node, nb, seqLen);
			cout << "build wc nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
		// non-watson crick helix neighbor
		else if((wcPairPosID[node->seqID] > 0 && nwcPairPosID[seqID+1] > 0  && wcPairPosID[node->seqID] == nwcPairPosID[seqID+1] + 1) ||
				(nwcPairPosID[node->seqID] > 0 && wcPairPosID[seqID+1] > 0  && nwcPairPosID[node->seqID] == wcPairPosID[seqID+1] + 1) ||
				(nwcPairPosID[node->seqID] > 0 && nwcPairPosID[seqID+1] > 0  && nwcPairPosID[node->seqID] == nwcPairPosID[seqID+1] + 1)) {
			BRConnection* ct = new BRConnection(fragLib, "nwcNb", node, nb, seqLen);
			node->midChild = nb;
			nb->father = node;
			cout << "build nwc nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
		// loop neighbor
		else {
			BRConnection* ct = new BRConnection(fragLib, "loopNb", node, nb, seqLen);
			node->midChild = nb;
			nb->father = node;
			cout << "build loop nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
	}

	/*
	 * build right node
	 */
	/*
	if(seqID + 1 < this->seqLen && (!connectToDownstream[seqID]) && nodes[seqID+1]->father == NULL){
		BRNode* nb = nodes[seqID+1];
		node->rightChild = nb;
		nb->father = node;
		BRConnection* ct = new BRConnection(fragLib, "jump" ,node, nb, seqLen);
		cout << "build jump connection: " << node->seqID << " " << nb->seqID << endl;
		bool isFixed = false;
		int idA = node->seqID;
		int idB = seqID + 1;
		for(int i=0;i<fixedGroups.size();i++) {
			if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
				isFixed = true;
		}
		if(isFixed || fixed[idB]) {
			ct->fixed = true;
			this->fixedConnectionList.push_back(ct);
		}
		else
			this->flexibleConnectionList.push_back(ct);
	}
	*/

	if(seqID < seqLen-1)
		buildFrom2(nodes[seqID+1]);
}

void BRFoldingTree::buildFrom3(BRNode* node) {
	/*
	 * base first, ribose second
	 */

	cout << "build node: " << node->seqID << endl;
	//if(node->father != NULL && node->seqID < seqLen-1)
	//	buildFrom2(nodes[node->seqID+1]);


	/*
	 * build left child: pairing node
	 */
	int seqID = node->seqID;

	if(wcPairPosID[seqID] > 0){
		BRNode* wc = nodes[wcPairPosID[seqID]];
		if(wc->father == NULL){
			node->leftChild = wc;
			wc->father = node;
			BRConnection* ct = new BRConnection(fragLib, "wc", node, wc, seqLen);
			cout << "build wc connection: " << node->seqID << " " << wc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = wcPairPosID[seqID];
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else {
				this->flexibleConnectionList.push_back(ct);
			}
		}
	}
	else if(nwcPairPosID[seqID] > 0) {
		BRNode* nwc = nodes[nwcPairPosID[seqID]];
		if(nwc->father == NULL){
			node->leftChild = nwc;
			nwc->father = node;
			BRConnection* ct = new BRConnection(fragLib, "nwc", node, nwc, seqLen);
			cout << "build nwc connection: " << node->seqID << " " << nwc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nwcPairPosID[seqID];
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
	}

	/*
	 * build middle node: neighbor pair
	 */
	if(seqID + 1 < this->seqLen && connectToDownstream[seqID] && nodes[seqID+1]->father == NULL) {
		BRNode* nb = nodes[seqID+1];
		if(loopRevPoints[seqID] == true) {
			cout << "stop at loop mid: " << seqID << endl;
			//do nothing
		}
		// helix neighbor pair
		else if(wcPairPosID[node->seqID] > 0 && wcPairPosID[seqID+1] > 0 &&
		   wcPairPosID[node->seqID] == wcPairPosID[seqID+1] + 1) {
			node->midChild = nb;
			nb->father = node;
			BRConnection* ct = new BRConnection(fragLib, "wcNb" ,node, nb, seqLen);
			cout << "build wc nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
		// non-watson crick helix neighbor
		else if((wcPairPosID[node->seqID] > 0 && nwcPairPosID[seqID+1] > 0  && wcPairPosID[node->seqID] == nwcPairPosID[seqID+1] + 1) ||
				(nwcPairPosID[node->seqID] > 0 && wcPairPosID[seqID+1] > 0  && nwcPairPosID[node->seqID] == wcPairPosID[seqID+1] + 1) ||
				(nwcPairPosID[node->seqID] > 0 && nwcPairPosID[seqID+1] > 0  && nwcPairPosID[node->seqID] == nwcPairPosID[seqID+1] + 1)) {
			BRConnection* ct = new BRConnection(fragLib, "nwcNb", node, nb, seqLen);
			node->midChild = nb;
			nb->father = node;
			cout << "build nwc nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
		// loop neighbor
		else {
			BRConnection* ct = new BRConnection(fragLib, "loopNb", node, nb, seqLen);
			node->midChild = nb;
			nb->father = node;
			cout << "build loop nb connection: " << node->seqID << " " << nb->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = seqID + 1;
			for(int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
	}

	/*
	 * build right node: jump
	 */
	/*
	if(seqID + 1 < this->seqLen && (!connectToDownstream[seqID]) && nodes[seqID+1]->father == NULL){
		BRNode* nb = nodes[seqID+1];
		node->rightChild = nb;
		nb->father = node;
		BRConnection* ct = new BRConnection(fragLib, "jump" ,node, nb, seqLen);
		cout << "build jump connection: " << node->seqID << " " << nb->seqID << endl;
		bool isFixed = false;
		int idA = node->seqID;
		int idB = seqID + 1;
		for(int i=0;i<fixedGroups.size();i++) {
			if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
				isFixed = true;
		}
		if(isFixed || fixed[idB]) {
			ct->fixed = true;
			this->fixedConnectionList.push_back(ct);
		}
		else
			this->flexibleConnectionList.push_back(ct);
	}
	*/


	if(seqID < seqLen-1)
		buildFrom3(nodes[seqID+1]);
}

void BRFoldingTree::buildReverseNodes() {
	for(int i=seqLen-2;i>=0;i--){
		if(nodes[i]->father == NULL){
			if(nodeConnectMatrix[nodes[i]->seqID*seqLen + nodes[i+1]->seqID]) continue;
			cout << "build reverse connection: " << i+1 << " " << i << endl;
			nodes[i]->father = nodes[i+1];
			nodes[i+1]->reverseChild = nodes[i];
			BRConnection* ct = new BRConnection(fragLib, "revNb", nodes[i+1], nodes[i], seqLen);
			bool isFixed = false;
			int idA = nodes[i+1]->seqID;
			int idB = nodes[i]->seqID;

			vector<int> ctList;
			for(int a=0;a<seqLen;a++){
				if(nodeConnectMatrix[a*seqLen+idA]){
					ctList.push_back(a);
				}
				else if(nodeConnectMatrix[a*seqLen+idB]){
					ctList.push_back(a);
				}
			}
			for(int a=0;a<ctList.size();a++){
				for(int b=0;b<ctList.size();b++){
					nodeConnectMatrix[ctList[a]*seqLen+ctList[b]] = true;
				}
			}
			//printBaseConnection();

			for(int k=0;k<fixedGroups.size();k++) {
				if(fixedGroups[k].find(idA) != fixedGroups[k].end() && fixedGroups[k].find(idB) != fixedGroups[k].end())
					isFixed = true;
			}
			if(isFixed || fixed[idB]) {
				ct->fixed = true;
				this->fixedConnectionList.push_back(ct);
			}
			else
				this->flexibleConnectionList.push_back(ct);
		}
	}
}

void BRFoldingTree::initFromKey(const string& keyInfo){
	vector<string> spt;
	splitString(keyInfo, " ", &spt);
	vector<int> rotIDs;
	vector<int> mvIDs;
	for(int i=2;i<this->seqLen+2;i++){
		rotIDs.push_back(atoi(spt[i].c_str()));
	}

	for(int i=this->seqLen+2;i<spt.size();i++){
		mvIDs.push_back(atoi(spt[i].c_str()));
	}

	for(int i=0;i<this->seqLen;i++){
		this->nodes[i]->rot = rotLib->rotLib[nodes[i]->baseType][rotIDs[i]];
		this->nodes[i]->tmpRot = this->nodes[i]->rot;

		nodes[i]->cs2 = nodes[i]->cs1 + nodes[i]->rot->mv12;
		nodes[i]->cs3 = nodes[i]->cs1 + nodes[i]->rot->mv13;
		nodes[i]->tmpCs2 = nodes[i]->cs2;
		nodes[i]->tmpCs3 = nodes[i]->cs3;

		for(int j=0;j<11;j++){
			nodes[i]->riboAtomCoords[j] = local2global(nodes[i]->cs1, nodes[i]->rot->tList1[j]);
			nodes[i]->riboAtomCoordsTmp[j] = local2global(nodes[i]->tmpCs1, nodes[i]->tmpRot->tList1[j]);
		}
	}

	for(int i=0;i<flexibleConnectionList.size();i++){
		int idA = mvIDs[i*2];
		int idB = mvIDs[i*2+1];
		F2Fragment* frag = flexibleConnectionList[i]->f2Lib->getFrag(idA, idB);
		flexibleConnectionList[i]->f2Frag = frag;
		flexibleConnectionList[i]->f2FragTmp = frag;
		flexibleConnectionList[i]->cm = frag->cm;
		flexibleConnectionList[i]->cmTmp = frag->cm;
	}

	updateCoordinate(rootNode);
	updatePhoGroups();
	updateEnergies(0.0);
}

void BRFoldingTree::keyInfoToRMS(const string& keyFile, const string& outFile){
	//cout << "key to rms" << endl;

	ifstream file;
	//cout << "open file: " << keyFile << endl;
	file.open(keyFile, ios::in);
	if(file.fail()){
		cout << "can't open file: " << keyFile << endl;
		exit(0);
	}

	ofstream out;
	out.open(outFile, ios::out);

	string line;
	vector<string> spt;

	vector<int> rotIDs;
	vector<int> mvIDs;
	while(getline(file, line)){
		//cout << line << endl;
		mvIDs.clear();
		splitString(line, " ", &spt);

		for(int i=2;i<this->seqLen+2;i++){
			rotIDs.push_back(atoi(spt[i].c_str()));
		}

		for(int i=this->seqLen+2;i<spt.size();i++){
			mvIDs.push_back(atoi(spt[i].c_str()));
		}

		for(int i=0;i<flexibleConnectionList.size();i++){
			int idA = mvIDs[i*2];
			int idB = mvIDs[i*2+1];
			F2Fragment* frag = flexibleConnectionList[i]->f2Lib->getFrag(idA, idB);
			flexibleConnectionList[i]->f2Frag = frag;
			flexibleConnectionList[i]->f2FragTmp = frag;
			flexibleConnectionList[i]->cm = frag->cm;
			flexibleConnectionList[i]->cmTmp = frag->cm;
		}
		updateCoordinate(rootNode);
		BRTreeInfo* info = getTreeInfo();
		double rms = info->rmsd(this->initTreeInfo);
		out << spt[0] << " " << spt[1] << " " << rms << endl;
	}

	file.close();
	out.close();
}

void BRFoldingTree::updatePhoGroups(){
	for(int i=0;i<seqLen-1;i++){
		BRNode* nodeA = nodes[i];
		BRNode* nodeB = nodes[i+1];
		if(!connectToDownstream[i]) continue;

		nodeA->phoLocal = this->et->pb->getPhoLocal(nodeA->cs2, nodeB->riboAtomCoords, nodeA->rot->improper, nodeB->rot->improper);
		nodeA->pho = PhophateGroup(nodeA->phoLocal, nodeA->cs2);
		nodeA->phoLocalTmp = nodeA->phoLocal;
		nodeA->phoTmp = nodeA->pho;
	}
}

void BRFoldingTree::updateEnergies(double shift){
	//cout << "update energy: shift=" << shift << endl;
	bool verbose = false;
	int i,j, pi, pj, sep;
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			pi = nodeA->seqID*seqLen+nodeB->seqID;
			pj = nodeB->seqID*seqLen+nodeA->seqID;

			allBaseClashE[pi] = baseBaseClash(nodeA, nodeB, sepTable[pi], et, shift, verbose);
			allBaseClashE[pj] = allBaseClashE[pi];

			allBaseBaseE[pi] = getBaseBaseEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
			allBaseBaseE[pj] = allBaseBaseE[pi];

			allBaseRiboseE[pi] = getBaseRiboseEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
			allBaseRiboseE[pj] = getBaseRiboseEnergy(nodeB, nodeA, sepTable[pj], et, verbose);

			allBasePhoE[pi] = getBasePhoEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
			allBasePhoE[pj] = getBasePhoEnergy(nodeB, nodeA, sepTable[pj], et, verbose);

			allRiboseRiboseE[pi] = getRiboseRiboseEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
			allRiboseRiboseE[pj] = allRibosePhoE[pi];

			allRibosePhoE[pi] = getRibosePhoEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
			allRibosePhoE[pj] = getRibosePhoEnergy(nodeB, nodeA, sepTable[pj], et, verbose);

			allPhoPhoE[pi] = getPhoPhoEnergy(nodeA, nodeB, sepTable[pi], et, verbose);
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
			tmpBaseClashE[i*seqLen+j] = allBaseClashE[i*seqLen+j];
			tmpBaseBaseE[i*seqLen+j] = allBaseBaseE[i*seqLen+j];
			tmpBaseRiboseE[i*seqLen+j] = allBaseRiboseE[i*seqLen+j];
			tmpBasePhoE[i*seqLen+j] = allBasePhoE[i*seqLen+j];
			tmpRiboseRiboseE[i*seqLen+j] = allRiboseRiboseE[i*seqLen+j];
			tmpRibosePhoE[i*seqLen+j] = allRibosePhoE[i*seqLen+j];
			tmpPhoPhoE[i*seqLen+j] = allPhoPhoE[i*seqLen+j];
		}
	}
}

string BRFoldingTree::toCtKey(){
	int flexCtNum = this->flexibleConnectionList.size();
	if(flexCtNum == 0) return "";
	string key = this->flexibleConnectionList[0]->f2Frag->mvType ;
	for(int i=1;i<flexCtNum;i++){
		key = key + flexibleConnectionList[i]->f2Frag->mvType;
	}
	return key;
}

string BRFoldingTree::toCtDetailString(){
	int flexCtNum = this->flexibleConnectionList.size();
	if(flexCtNum == 0) return "";
	char xx[20];
	vector<string> ids;

	for(int i=0;i<seqLen;i++){
		sprintf(xx, "%d", nodes[i]->rot->rotType);
		ids.push_back(string(xx));
	}

	for(int i=0;i<flexCtNum;i++){
		sprintf(xx, "%d %d", this->flexibleConnectionList[i]->f2Frag->level1ID, this->flexibleConnectionList[i]->f2Frag->level2ID);
		ids.push_back(string(xx));
	}
	string key = ids[0];
	for(int i=1;i<ids.size();i++){
		key = key + " " + ids[i];
	}


	return key;
}

void BRFoldingTree::updateCtChildTmpCs(BRConnection* ct, CsMove& cm, bool verbose){
	int i,j,k;
	BRNode* childNode = ct->childNode;
	BRNode* fatherNode = ct->fatherNode;

	if(verbose){
		cout << "update connection child coordinates: " << fatherNode->seqID << "->" << childNode->seqID << endl;
	}

	ct->cmTmp = cm;
	childNode->tmpCs1 = fatherNode->tmpCs1 + cm;
	childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
	childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;


	for(int i=0;i<11;i++){
		childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
	}

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
	}

	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnection* cct = ct->childConnectionList[j];

		CsMove cmv = cct->cmTmp;
		childNode = cct->childNode;
		fatherNode = cct->fatherNode;
		if(childNode->fixed) continue;

		if(verbose){
			cout << "update connection child coordinates2: " << fatherNode->seqID << "->" << childNode->seqID << endl;
		}

		childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;
		childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
		childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
		}
		for(i=0;i<11;i++){
			childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
		}
	}


	int phoChangeNum = ct->ctPhoGroupC.size();

	BRNode* nodeA;
	BRNode* nodeB;
	int indexA, indexB;
	for(i=0;i<phoChangeNum;i++){
		indexA = ct->ctPhoGroupC[i];
		indexB = indexA+1;

		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);

		if(verbose){
			cout << "update pho: " << indexA << endl;
		}

	}

	int upPhoNum = ct->ctPhoGroupB.size();
	int index;
	for(i=0;i<upPhoNum;i++){
		index = ct->ctPhoGroupB[i];
		nodes[index]->phoTmp = PhophateGroup(nodes[index]->phoLocalTmp, nodes[index]->tmpCs2);
	}
}

void BRFoldingTree::clearCtChildTmpCs(BRConnection* ct, bool verbose){

	if(verbose){
		cout << "clear connection child tmp coordinates: " << ct->fatherNode->seqID << ct->childNode->seqID << endl;
	}

	ct->cmTmp = ct->cm;

	int i,j,k, pi;
	BRNode* node;

	for(i=0;i<ct->ctBaseGroupB.size();i++){
		node = nodes[ct->ctBaseGroupB[i]];
		node->tmpCs1 = node->cs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
	}

	for(i=0;i<ct->ctRiboseGroupB.size();i++){
		node = nodes[ct->ctRiboseGroupB[i]];
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(j=0;j<11;j++){
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
	}

	for(i=0;i<ct->ctPhoGroupB.size();i++){
		node = nodes[ct->ctPhoGroupB[i]];
		node->phoTmp = node->pho;
	}

	for(i=0;i<ct->ctPhoGroupC.size();i++){
		node = nodes[ct->ctPhoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}
}


void BRFoldingTree::acceptCtChildTmpCs(BRConnection* ct, bool verbose){
	if(verbose){
		cout << "accept connection child tmp coordinates: " << ct->fatherNode->seqID << ct->childNode->seqID << endl;
	}

	ct->cm = ct->cmTmp;
	int i,j,k, pi;
	BRNode* node;

	for(i=0;i<ct->ctBaseGroupB.size();i++){
		node = nodes[ct->ctBaseGroupB[i]];
		node->cs1 = node->tmpCs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = node->baseAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct->ctRiboseGroupB.size();i++){
		node = nodes[ct->ctRiboseGroupB[i]];
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct->ctPhoGroupB.size();i++){
		node = nodes[ct->ctPhoGroupB[i]];
		node->pho = node->phoTmp;
	}

	for(i=0;i<ct->ctPhoGroupC.size();i++){
		node = nodes[ct->ctPhoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}

	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::updateF2ChildCs(BRConnection* ct, F2Fragment* frag, bool verbose){
	int i,j,k;
	BRNode* childNode = ct->childNode;
	BRNode* fatherNode = ct->fatherNode;
	ct->cm = frag->cm;
	ct->f2Frag = frag;

	childNode->cs1 = fatherNode->cs1 + frag->cm;
	fatherNode->cs2 = fatherNode->cs1 + fatherNode->rot->mv12;
	fatherNode->cs3 = fatherNode->cs1 + fatherNode->rot->mv13;
	childNode->cs2 = childNode->cs1 + childNode->rot->mv12;
	childNode->cs3 = childNode->cs1 + childNode->rot->mv13;

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
	}

	for(int i=0;i<11;i++){
		childNode->riboAtomCoords[i] = local2global(childNode->cs1, childNode->rot->tList1[i]);
	}

	if(verbose){
		cout << "update connection child coordinates: " << fatherNode->seqID << "->" << childNode->seqID << endl;
	}

	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnection* cct = ct->childConnectionList[j];

		CsMove cmv = cct->cm;
		BRNode* childNode = cct->childNode;
		BRNode* fatherNode = cct->fatherNode;

		if(childNode->fixed) continue;

		childNode->cs1 = fatherNode->cs1 + cmv;
		childNode->cs2 = childNode->cs1 + childNode->rot->mv12;
		childNode->cs3 = childNode->cs1 + childNode->rot->mv13;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
		}
		for(i=0;i<11;i++){
			childNode->riboAtomCoords[i] = local2global(childNode->cs1, childNode->rot->tList1[i]);
		}
	}

	ct->childRot = childNode->rot;
	ct->childRotTmp = childNode->rot;

	int chainBreakNum = ct->f2PhoGroupC.size();

	BRNode* nodeA;
	BRNode* nodeB;
	int indexA, indexB;
	for(i=0;i<chainBreakNum;i++){
		indexA = ct->f2PhoGroupC[i];
		indexB = indexA+1;
		if(verbose){
			cout << "connection: " << indexA << " " << indexB << endl;
		}
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		nodeA->phoLocal = this->et->pb->getPhoLocal(nodeA->cs2, nodeB->riboAtomCoords, nodeA->rot->improper, nodeB->rot->improper);

		nodeA->pho = PhophateGroup(nodeA->phoLocal, nodeA->cs2);
	}

	int upPhoNum = ct->f2PhoGroupB.size();
	int index;
	for(i=0;i<upPhoNum;i++){
		index = ct->f2PhoGroupB[i];
		nodes[index]->pho = PhophateGroup(nodes[index]->phoLocal, nodes[index]->cs2);
	}
}

void BRFoldingTree::updateF2ChildTmpCs(BRConnection* ct, F2Fragment* frag, bool verbose){
	int i,j,k;
	BRNode* childNode = ct->childNode;
	BRNode* fatherNode = ct->fatherNode;

	if(verbose){
		cout << "update connection child coordinates: " << fatherNode->seqID << "->" << childNode->seqID << endl;
	}

	ct->cmTmp = frag->cm;
	childNode->tmpCs1 = fatherNode->tmpCs1 + frag->cm;
	childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
	childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;
	ct->f2FragTmp = frag;

	for(int i=0;i<11;i++){
		childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
	}

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
	}

	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnection* cct = ct->childConnectionList[j];

		CsMove cmv = cct->cmTmp;
		childNode = cct->childNode;
		fatherNode = cct->fatherNode;
		if(childNode->fixed) continue;

		if(verbose){
			cout << "update connection child coordinates2: " << fatherNode->seqID << "->" << childNode->seqID << endl;
		}

		childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;
		childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
		childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
		}
		for(i=0;i<11;i++){
			childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
		}
	}

	BRNode* nodeA;
	BRNode* nodeB;

	if(frag->hasRibose && frag->baseFlip) {
		if(verbose){
			cout << "has flip" << endl;
		}
		BaseRotamer* minEChildRot;
		double minE = 99999.9;
		double e;
		childNode = ct->childNode;
		for(i=0;i<5;i++){
			childNode->tmpRot = rotLib->getFlipRotamer(seq[childNode->baseType]);
			e = childNode->tmpRot->energy;
			for(int i=0;i<11;i++){
				childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
			}
			childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
			childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;

			if(childNode->seqID > 0 && connectToDownstream[childNode->seqID-1]){
				nodeA = nodes[childNode->seqID-1];
				nodeB = childNode;
				e += this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper).e;
				e += getRiboseRiboseEnergyTmp(nodeA, nodeB, 1, this->et, false);
			}


			if(connectToDownstream[childNode->seqID]){
				nodeA = childNode;
				nodeB = nodes[childNode->seqID+1];
				e += this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper).e;
				e += getRiboseRiboseEnergyTmp(nodeA, nodeB, 1, this->et, false);
			}
			if(e < minE){
				minE = e;
				minEChildRot = childNode->tmpRot;
			}

			if(verbose){
				printf("rot imp: %6.2f chi: %6.2f ene: %8.3f\n", childNode->tmpRot->improper, childNode->tmpRot->chi, e);
			}
		}

		childNode->tmpRot = minEChildRot;
		for(int i=0;i<11;i++){
			childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
		}
		childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
		childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;
		ct->childRotTmp = minEChildRot;
	}
	else if(frag->hasRibose && rand()%3 == 0){
		BaseRotamer* minEChildRot;
		double minE = 99999.9;
		double e;
		childNode = ct->childNode;
		for(i=0;i<5;i++){
			childNode->tmpRot = rotLib->getRandomRotamerLv1(seq[childNode->baseType]);
			e = childNode->tmpRot->energy;
			for(int i=0;i<11;i++){
				childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
			}
			childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
			childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;

			if(childNode->seqID > 0 && connectToDownstream[childNode->seqID-1]){
				nodeA = nodes[childNode->seqID-1];
				nodeB = childNode;
				e += this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper).e;
			}


			if(connectToDownstream[childNode->seqID]){
				nodeA = childNode;
				nodeB = nodes[childNode->seqID+1];
				e += this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper).e;
			}
			if(e < minE){
				minE = e;
				minEChildRot = childNode->tmpRot;
			}
		}

		childNode->tmpRot = minEChildRot;
		for(int i=0;i<11;i++){
			childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
		}
		childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
		childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;
		ct->childRotTmp = minEChildRot;
	}



	int phoChangeNum = ct->f2PhoGroupC.size();
	int indexA, indexB;
	for(i=0;i<phoChangeNum;i++){
		indexA = ct->f2PhoGroupC[i];
		indexB = indexA+1;
		if(verbose){
			cout << "update pho: " << indexA << " " << indexB << endl;
		}
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
	}

	int upPhoNum = ct->f2PhoGroupB.size();
	int index;
	for(i=0;i<upPhoNum;i++){
		index = ct->f2PhoGroupB[i];
		nodes[index]->phoTmp = PhophateGroup(nodes[index]->phoLocalTmp, nodes[index]->tmpCs2);
	}
}

void BRFoldingTree::clearF2ChildTmpCs(BRConnection* ct, bool verbose){

	if(verbose){
		cout << "clear connection child tmp coordinates: " << ct->fatherNode->seqID << ct->childNode->seqID << endl;
	}

	ct->cmTmp = ct->cm;
	ct->f2FragTmp = ct->f2Frag;
	ct->childRotTmp = ct->childRot;

	int i,j,k, pi;
	BRNode* node;

	for(i=0;i<ct->f2BaseGroupB.size();i++){
		node = nodes[ct->f2BaseGroupB[i]];
		node->tmpCs1 = node->cs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
	}

	for(i=0;i<ct->f2RiboseGroupB.size();i++){
		node = nodes[ct->f2RiboseGroupB[i]];
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(j=0;j<11;j++){
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
	}

	for(i=0;i<ct->f2RiboseGroupC.size();i++){
		node = nodes[ct->f2RiboseGroupC[i]];
		node->tmpRot = node->rot;
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(j=0;j<11;j++){
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
	}

	for(i=0;i<ct->f2PhoGroupB.size();i++){
		node = nodes[ct->f2PhoGroupB[i]];
		node->phoTmp = node->pho;
	}

	for(i=0;i<ct->f2PhoGroupC.size();i++){
		node = nodes[ct->f2PhoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::acceptF2ChildTmpCs(BRConnection* ct, bool verbose){
	if(verbose){
		cout << "accept connection child tmp coordinates: " << ct->fatherNode->seqID << ct->childNode->seqID << endl;
	}

	ct->cm = ct->cmTmp;
	ct->f2Frag = ct->f2FragTmp;
	ct->childRot = ct->childRotTmp;

	int i,j,k, pi;
	BRNode* node;

	for(i=0;i<ct->f2BaseGroupB.size();i++){
		node = nodes[ct->f2BaseGroupB[i]];
		node->cs1 = node->tmpCs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = node->baseAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct->f2RiboseGroupB.size();i++){
		node = nodes[ct->f2RiboseGroupB[i]];
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct->f2RiboseGroupC.size();i++){
		node = nodes[ct->f2RiboseGroupC[i]];
		node->rot = node->tmpRot;
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct->f2PhoGroupB.size();i++){
		node = nodes[ct->f2PhoGroupB[i]];
		node->pho = node->phoTmp;
	}

	for(i=0;i<ct->f2PhoGroupC.size();i++){
		node = nodes[ct->f2PhoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}

	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::updateF3ChildTmpCs(BRConnection* ct1, F3Fragment* frag, bool verbose){
	BRNode* father = ct1->fatherNode;
	BRNode* child = ct1->childNode;
	BRNode* grandChild = child->midChild;

	if(child->fixed) return;
	if(grandChild->fixed) return;

	ct1->cmTmp = frag->cmAB;
	grandChild->upConnection->cmTmp = frag->cmBC;

	father->tmpRot = frag->rotA;
	child->tmpRot = frag->rotB;
	grandChild->tmpRot = frag->rotC;

	int i,j;

	if(verbose){
		cout << "update base: " << child->seqID << endl;
		cout << "update base: " << grandChild->seqID << endl;
	}

	father->tmpCs2 = father->tmpCs1 + father->tmpRot->mv12;
	father->tmpCs3 = father->tmpCs1 + father->tmpRot->mv13;

	child->tmpCs1 = father->tmpCs1 + frag->cmAB;
	child->tmpCs2 = child->tmpCs1 + child->tmpRot->mv12;
	child->tmpCs3 = child->tmpCs1 + child->tmpRot->mv13;

	grandChild->tmpCs1 = child->tmpCs1 + frag->cmBC;
	grandChild->tmpCs2 = grandChild->tmpCs1 + grandChild->tmpRot->mv12;
	grandChild->tmpCs3 = grandChild->tmpCs1 + grandChild->tmpRot->mv13;



	if(verbose){
		cout << "update ribose: " << father->seqID << endl;
		cout << "update ribose: " << child->seqID << endl;
		cout << "update ribose: " << grandChild->seqID << endl;
	}

	for(i=0;i<11;i++){
		father->riboAtomCoordsTmp[i] = local2global(father->tmpCs1, frag->rotA->tList1[i]);
		child->riboAtomCoordsTmp[i] = local2global(child->tmpCs1, frag->rotB->tList1[i]);
		grandChild->riboAtomCoordsTmp[i] = local2global(grandChild->tmpCs1, frag->rotC->tList1[i]);
	}

	for(i=0;i<child->baseAtomNum;i++){
		child->baseAtomCoordsTmp[i] = local2global(child->tmpCs1, child->atomCoordLocal[i]);
	}
	for(i=0;i<grandChild->baseAtomNum;i++){
		grandChild->baseAtomCoordsTmp[i] = local2global(grandChild->tmpCs1, grandChild->atomCoordLocal[i]);
	}

	for(j=0;j<ct1->childConnectionList.size();j++){

		BRConnection* cct = ct1->childConnectionList[j];
		if(verbose){
					cout << "ct: " << cct->fatherNode->seqID << " " << cct->childNode->seqID << endl;
		}
		CsMove cmv = cct->cm;
		BRNode* childNode = cct->childNode;
		BRNode* fatherNode = cct->fatherNode;
		if(childNode->fixed) continue;

		if(childNode->seqID == grandChild->seqID) continue;
		childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;
		childNode->tmpCs2 = childNode->tmpCs1 + childNode->tmpRot->mv12;
		childNode->tmpCs3 = childNode->tmpCs1 + childNode->tmpRot->mv13;

		if(verbose){
			cout << "update node: " << childNode->seqID << endl;
		}

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
		}
		for(i=0;i<11;i++){
			childNode->riboAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->tmpRot->tList1[i]);
		}
	}

	int chainBreakNum = ct1->f3PhoGroupC.size();
	int indexA, indexB;
	double minE;
	int minIndexA1, minIndexB1, minIndexA2, minIndexB2, d1d2Index;
	int riboKeyIndex;
	double phoE;
	BaseRotamer* rotA;
	BaseRotamer* rotB;
	BRNode* nodeA;
	BRNode* nodeB;
	for(int i=0;i<chainBreakNum;i++){

		indexA = ct1->f3PhoGroupC[i];
		indexB = indexA+1;
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		if(verbose){
			cout << "update pho: " << indexA << endl;
		}

		if(indexA == ct1->fatherNode->seqID){
			nodeA->phoLocalTmp = frag->plA;
			nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
		}
		else if(indexA == ct1->childNode->seqID){
			nodeA->phoLocalTmp = frag->plB;
			nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
		}
		else {
			nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
			nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
		}
	}

	int upPhoNum = ct1->f3PhoGroupB.size();
	int index;

	//cout << "pho without local change: " << upPhoNum << endl;
	for(int i=0;i<upPhoNum;i++){

		index = ct1->f3PhoGroupB[i];
		if(verbose){
					cout << "update pho: " <<  index << " without local change" << endl;
			}
		nodes[index]->phoTmp = PhophateGroup(nodes[index]->phoLocalTmp, nodes[index]->tmpCs2);
	}

}

void BRFoldingTree::clearF3ChildTmpCs(BRConnection* ct1, bool verbose){

	BRNode* father = ct1->fatherNode;
	BRNode* child = ct1->childNode;
	BRNode* grandChild = child->midChild;

	father->tmpRot = father->rot;
	child->tmpRot = child->rot;
	grandChild->tmpRot = grandChild->rot;

	int i,j,k, pi;
	BRNode* node;

	ct1->cmTmp = ct1->cm;
	BRConnection* ct2 = ct1->childNode->midChild->upConnection;
	ct2->cmTmp = ct2->cm;

	for(i=0;i<ct1->f3BaseGroupB.size();i++){
		node = nodes[ct1->f3BaseGroupB[i]];
		node->tmpCs1 = node->cs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
	}

	for(i=0;i<ct1->f3BaseGroupC.size();i++){
		node = nodes[ct1->f3BaseGroupC[i]];
		node->tmpCs1 = node->cs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoordsTmp[j] = node->baseAtomCoords[j];
		}
	}

	for(i=0;i<ct1->f3RiboseGroupB.size();i++){
		node = nodes[ct1->f3RiboseGroupB[i]];
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(j=0;j<11;j++){
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
	}

	for(i=0;i<ct1->f3RiboseGroupC.size();i++){
		node = nodes[ct1->f3RiboseGroupC[i]];
		node->tmpRot = node->rot;
		node->tmpCs2 = node->cs2;
		node->tmpCs3 = node->cs3;
		for(j=0;j<11;j++){
			node->riboAtomCoordsTmp[j] = node->riboAtomCoords[j];
		}
	}

	for(i=0;i<ct1->f3PhoGroupB.size();i++){
		node = nodes[ct1->f3PhoGroupB[i]];
		node->phoTmp = node->pho;
	}

	for(i=0;i<ct1->f3PhoGroupC.size();i++){
		node = nodes[ct1->f3PhoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::acceptF3ChildTmpCs(BRConnection* ct1, bool verbose){
	BRNode* father = ct1->fatherNode;
	BRNode* child = ct1->childNode;
	BRNode* grandChild = child->midChild;

	father->rot = father->tmpRot;
	child->rot = child->tmpRot;
	grandChild->rot = grandChild->tmpRot;

	int i,j,k, pi;
	BRNode* node;

	ct1->cm = ct1->cmTmp;
	BRConnection* ct2 = ct1->childNode->midChild->upConnection;
	ct2->cm = ct2->cmTmp;

	for(i=0;i<ct1->f3BaseGroupB.size();i++){
		node = nodes[ct1->f3BaseGroupB[i]];
		node->cs1 = node->tmpCs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = node->baseAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct1->f3BaseGroupC.size();i++){
		node = nodes[ct1->f3BaseGroupC[i]];
		node->cs1 = node->tmpCs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = node->baseAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct1->f3RiboseGroupB.size();i++){
		node = nodes[ct1->f3RiboseGroupB[i]];
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct1->f3RiboseGroupC.size();i++){
		node = nodes[ct1->f3RiboseGroupC[i]];
		node->rot = node->tmpRot;
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<ct1->f3PhoGroupB.size();i++){
		node = nodes[ct1->f3PhoGroupB[i]];
		node->pho = node->phoTmp;
	}

	for(i=0;i<ct1->f3PhoGroupC.size();i++){
		node = nodes[ct1->f3PhoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}

	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::updateBaseRotamerTmp(BRNode* node, BaseRotamer* rot,  bool verbose){

	node->tmpRot = rot;
	for(int i=0;i<11;i++){
		node->riboAtomCoordsTmp[i] = local2global(node->tmpCs1, node->tmpRot->tList1[i]);
	}
	node->tmpCs2 = node->tmpCs1 + rot->mv12;
	node->tmpCs3 = node->tmpCs1 + rot->mv13;

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

		nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
	}
}

void BRFoldingTree::clearBaseRotamerTmp(BRNode* nodeSelect, bool verbose){

	int i,j;
	nodeSelect->tmpRot = nodeSelect->rot;
	nodeSelect->tmpCs2 = nodeSelect->cs2;
	nodeSelect->tmpCs3 = nodeSelect->cs3;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoordsTmp[j] = nodeSelect->riboAtomCoords[j];
	}
	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}
	int pi;
	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::acceptBaseRotamerTmp(BRNode* nodeSelect, bool verbose){
	int i,j;
	nodeSelect->rot = nodeSelect->tmpRot;
	nodeSelect->cs2 = nodeSelect->tmpCs2;
	nodeSelect->cs3 = nodeSelect->tmpCs3;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoords[j] = nodeSelect->riboAtomCoordsTmp[j];
	}
	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}
	int pi;
	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::updateReverseBaseRotamerTmp(BRNode* node, BaseRotamer* rot, bool verbose){
	//cs2 not changed

	node->tmpCs1 = node->tmpCs2 + rot->mv21;
	node->tmpCs3 = node->tmpCs1 + rot->mv13;
	node->tmpRot = rot;
	for(int i=0;i<node->baseAtomNum;i++){
		node->baseAtomCoordsTmp[i] = local2global(node->tmpCs1, node->atomCoordLocal[i]);
	}
	for(int i=0;i<11;i++){
		node->riboAtomCoordsTmp[i] = local2global(node->tmpCs1, node->tmpRot->tList1[i]);
	}

	if(node->upConnection != NULL)
		node->upConnection->cmTmp = node->tmpCs1 - node->upConnection->fatherNode->tmpCs1;
	if(node->leftChild != NULL){
		node->leftChild->upConnection->cmTmp = node->leftChild->tmpCs1 - node->tmpCs1;
	}
	if(node->midChild != NULL){
		node->midChild->upConnection->cmTmp = node->midChild->tmpCs1 - node->tmpCs1;
	}
	if(node->rightChild != NULL){
		node->rightChild->upConnection->cmTmp = node->rightChild->tmpCs1 - node->tmpCs1;
	}
	if(node->bulge13Child != NULL) {
		node->bulge13Child->upConnection->cmTmp = node->bulge13Child->tmpCs1 - node->tmpCs1;
	}
	if(node->bulge14Child != NULL){
		node->bulge14Child->upConnection->cmTmp = node->bulge14Child->tmpCs1 - node->tmpCs1;
	}
	if(node->reverseChild != NULL){
		node->reverseChild->upConnection->cmTmp = node->reverseChild->tmpCs1 - node->tmpCs1;
	}
	if(node->revBulge13Child != NULL) {
		node->revBulge13Child->upConnection->cmTmp = node->revBulge13Child->tmpCs1 - node->tmpCs1;
	}
	if(node->revBulge14Child != NULL){
		node->revBulge14Child->upConnection->cmTmp = node->revBulge14Child->tmpCs1 - node->tmpCs1;
	}
	int chainBreakNum = node->phoGroupC.size();
	int indexA, indexB;

	BRNode* nodeA;
	BRNode* nodeB;
	for(int i=0;i<chainBreakNum;i++){

		indexA = node->phoGroupC[i];
		indexB = indexA+1;
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		if(verbose){
			cout << "update pho: " << indexA << endl;
		}

		nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
	}

}

void BRFoldingTree::clearReverseBaseRotamerTmp(BRNode* nodeSelect, bool verbose){
	int i,j;
	nodeSelect->tmpRot = nodeSelect->rot;
	nodeSelect->tmpCs1 = nodeSelect->cs1;
	nodeSelect->tmpCs3 = nodeSelect->cs3;

	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoordsTmp[j] = nodeSelect->riboAtomCoords[j];
	}
	for(j=0;j<nodeSelect->baseAtomNum;j++){
		nodeSelect->baseAtomCoordsTmp[j] = nodeSelect->baseAtomCoords[j];
	}

	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	if(nodeSelect->upConnection != NULL)
		nodeSelect->upConnection->cmTmp = nodeSelect->upConnection->cm;
	if(nodeSelect->leftChild != NULL){
		nodeSelect->leftChild->upConnection->cmTmp = nodeSelect->leftChild->upConnection->cm;
	}
	if(nodeSelect->midChild != NULL){
		nodeSelect->midChild->upConnection->cmTmp = nodeSelect->midChild->upConnection->cm;
	}
	if(nodeSelect->rightChild != NULL){
		nodeSelect->rightChild->upConnection->cmTmp = nodeSelect->rightChild->upConnection->cm;
	}
	if(nodeSelect->bulge13Child != NULL){
		nodeSelect->bulge13Child->upConnection->cmTmp = nodeSelect->bulge13Child->upConnection->cm;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cmTmp = nodeSelect->bulge14Child->upConnection->cm;
	}
	if(nodeSelect->reverseChild != NULL){
		nodeSelect->reverseChild->upConnection->cmTmp = nodeSelect->reverseChild->upConnection->cm;
	}
	if(nodeSelect->revBulge13Child != NULL){
		nodeSelect->revBulge13Child->upConnection->cmTmp = nodeSelect->revBulge13Child->upConnection->cm;
	}
	if(nodeSelect->revBulge14Child != NULL){
		nodeSelect->revBulge14Child->upConnection->cmTmp = nodeSelect->revBulge14Child->upConnection->cm;
	}

	int pi;

	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}


}

void BRFoldingTree::acceptReverseBaseRotamerTmp(BRNode* nodeSelect, bool verbose){
	int i,j;
	nodeSelect->rot = nodeSelect->tmpRot;
	nodeSelect->cs1 = nodeSelect->tmpCs1;
	nodeSelect->cs3 = nodeSelect->tmpCs3;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoords[j] = nodeSelect->riboAtomCoordsTmp[j];
	}
	for(j=0;j<nodeSelect->baseAtomNum;j++){
		nodeSelect->baseAtomCoords[j] = nodeSelect->baseAtomCoordsTmp[j];
	}

	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocal = node->phoLocalTmp;
		node->pho = node->phoTmp;
	}

	if(nodeSelect->upConnection != NULL)
		nodeSelect->upConnection->cm = nodeSelect->upConnection->cmTmp;
	if(nodeSelect->leftChild != NULL){
		nodeSelect->leftChild->upConnection->cm = nodeSelect->leftChild->upConnection->cmTmp;
	}
	if(nodeSelect->midChild != NULL){
		nodeSelect->midChild->upConnection->cm = nodeSelect->midChild->upConnection->cmTmp;
	}
	if(nodeSelect->rightChild != NULL){
		nodeSelect->rightChild->upConnection->cm = nodeSelect->rightChild->upConnection->cmTmp;
	}
	if(nodeSelect->bulge13Child != NULL){
		nodeSelect->bulge13Child->upConnection->cm = nodeSelect->bulge13Child->upConnection->cmTmp;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cm = nodeSelect->bulge14Child->upConnection->cmTmp;
	}
	if(nodeSelect->reverseChild != NULL){
		nodeSelect->reverseChild->upConnection->cm = nodeSelect->reverseChild->upConnection->cmTmp;
	}
	if(nodeSelect->revBulge13Child != NULL){
		nodeSelect->revBulge13Child->upConnection->cm = nodeSelect->revBulge13Child->upConnection->cmTmp;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cm = nodeSelect->bulge14Child->upConnection->cmTmp;
	}

	int pi;

	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

void BRFoldingTree::updateSingleBaseCoordTmp(BRNode* node, CsMove& move, bool verbose){

	node->tmpCs1 = node->cs1 + move;
	node->tmpCs2 = node->tmpCs1 + node->tmpRot->mv12;
	node->tmpCs3 = node->tmpCs1 + node->tmpRot->mv13;

	if(node->upConnection != NULL)
		node->upConnection->cmTmp = node->tmpCs1 - node->upConnection->fatherNode->tmpCs1;
	if(node->leftChild != NULL){
		node->leftChild->upConnection->cmTmp = node->leftChild->tmpCs1 - node->tmpCs1;
	}
	if(node->midChild != NULL){
		node->midChild->upConnection->cmTmp = node->midChild->tmpCs1 - node->tmpCs1;
	}
	if(node->rightChild != NULL){
		node->rightChild->upConnection->cmTmp = node->rightChild->tmpCs1 - node->tmpCs1;
	}
	if(node->bulge13Child != NULL) {
		node->bulge13Child->upConnection->cmTmp = node->bulge13Child->tmpCs1 - node->tmpCs1;
	}
	if(node->bulge14Child != NULL){
		node->bulge14Child->upConnection->cmTmp = node->bulge14Child->tmpCs1 - node->tmpCs1;
	}
	if(node->reverseChild != NULL){
		node->reverseChild->upConnection->cmTmp = node->reverseChild->tmpCs1 - node->tmpCs1;
	}
	if(node->revBulge13Child != NULL) {
		node->revBulge13Child->upConnection->cmTmp = node->revBulge13Child->tmpCs1 - node->tmpCs1;
	}
	if(node->revBulge14Child != NULL){
		node->revBulge14Child->upConnection->cmTmp = node->revBulge14Child->tmpCs1 - node->tmpCs1;

	}

	for(int i=0;i<node->baseAtomNum;i++){
		node->baseAtomCoordsTmp[i] = local2global(node->tmpCs1, node->atomCoordLocal[i]);
	}

	for(int i=0;i<11;i++){
		node->riboAtomCoordsTmp[i] = local2global(node->tmpCs1, node->rot->tList1[i]);
	}

	int chainBreakNum = node->phoGroupC.size();
	int indexA, indexB;

	BRNode* nodeA;
	BRNode* nodeB;
	for(int i=0;i<chainBreakNum;i++){

		indexA = node->phoGroupC[i];
		indexB = indexA+1;
		nodeA = nodes[indexA];
		nodeB = nodes[indexB];

		if(verbose){
			cout << "update pho: " << indexA << endl;
		}

		nodeA->phoLocalTmp = this->et->pb->getPhoLocal(nodeA->tmpCs2, nodeB->riboAtomCoordsTmp, nodeA->tmpRot->improper, nodeB->tmpRot->improper);
		nodeA->phoTmp = PhophateGroup(nodeA->phoLocalTmp, nodeA->tmpCs2);
	}
}

void BRFoldingTree::clearSingleBaseCoordTmp(BRNode* nodeSelect, bool verbose){
	nodeSelect->tmpCs1 = nodeSelect->cs1;
	nodeSelect->tmpCs2 = nodeSelect->cs2;
	nodeSelect->tmpCs3 = nodeSelect->cs3;

	if(nodeSelect->upConnection != NULL)
		nodeSelect->upConnection->cmTmp = nodeSelect->upConnection->cm;
	if(nodeSelect->leftChild != NULL){
		nodeSelect->leftChild->upConnection->cmTmp = nodeSelect->leftChild->upConnection->cm;
	}
	if(nodeSelect->midChild != NULL){
		nodeSelect->midChild->upConnection->cmTmp = nodeSelect->midChild->upConnection->cm;
	}
	if(nodeSelect->rightChild != NULL){
		nodeSelect->rightChild->upConnection->cmTmp = nodeSelect->rightChild->upConnection->cm;
	}
	if(nodeSelect->bulge13Child != NULL){
		nodeSelect->bulge13Child->upConnection->cmTmp = nodeSelect->bulge13Child->upConnection->cm;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cmTmp = nodeSelect->bulge14Child->upConnection->cm;
	}
	if(nodeSelect->reverseChild != NULL){
		nodeSelect->reverseChild->upConnection->cmTmp = nodeSelect->reverseChild->upConnection->cm;
	}
	if(nodeSelect->revBulge13Child != NULL){
		nodeSelect->revBulge13Child->upConnection->cmTmp = nodeSelect->revBulge13Child->upConnection->cm;
	}
	if(nodeSelect->revBulge14Child != NULL){
		nodeSelect->revBulge14Child->upConnection->cmTmp = nodeSelect->revBulge14Child->upConnection->cm;
	}

	int i,j;


	nodeSelect->tmpCs1 = nodeSelect->cs1;
	for(j=0;j<nodeSelect->baseAtomNum;j++){
		nodeSelect->baseAtomCoordsTmp[j] = nodeSelect->baseAtomCoords[j];
	}

	nodeSelect->tmpCs2 = nodeSelect->cs2;
	nodeSelect->tmpCs3 = nodeSelect->cs3;
	for(j=0;j<11;j++){
		nodeSelect->riboAtomCoordsTmp[j] = nodeSelect->riboAtomCoords[j];
	}

	BRNode* node;
	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->phoLocalTmp = node->phoLocal;
		node->phoTmp = node->pho;
	}

	int pi;

	for(i=0;i<seqLen;i++){
		tmpRotE[i] = allRotE[i];
		tmpRcE[i] = allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseClashE[pi] = allBaseClashE[pi];
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
			tmpBaseRiboseE[pi] = allBaseRiboseE[pi];
			tmpBasePhoE[pi] = allBasePhoE[pi];
			tmpRiboseRiboseE[pi] = allRiboseRiboseE[pi];
			tmpRibosePhoE[pi] = allRibosePhoE[pi];
			tmpPhoPhoE[pi] = allPhoPhoE[pi];
		}
	}

}

void BRFoldingTree::acceptSingleBaseCoordTmp(BRNode* nodeSelect, bool verbose){

	if(verbose){
		cout << "accept node: " << nodeSelect->seqID << endl;
	}
	nodeSelect->cs1 = nodeSelect->tmpCs1;
	nodeSelect->cs2 = nodeSelect->tmpCs2;
	nodeSelect->cs3 = nodeSelect->tmpCs3;

	if(nodeSelect->upConnection != NULL)
		nodeSelect->upConnection->cm = nodeSelect->upConnection->cmTmp;
	if(nodeSelect->leftChild != NULL){
		nodeSelect->leftChild->upConnection->cm = nodeSelect->leftChild->upConnection->cmTmp;
	}
	if(nodeSelect->midChild != NULL){
		nodeSelect->midChild->upConnection->cm = nodeSelect->midChild->upConnection->cmTmp;
	}
	if(nodeSelect->rightChild != NULL){
		nodeSelect->rightChild->upConnection->cm = nodeSelect->rightChild->upConnection->cmTmp;
	}
	if(nodeSelect->bulge13Child != NULL){
		nodeSelect->bulge13Child->upConnection->cm = nodeSelect->bulge13Child->upConnection->cmTmp;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cm = nodeSelect->bulge14Child->upConnection->cmTmp;
	}
	if(nodeSelect->reverseChild != NULL){
		nodeSelect->reverseChild->upConnection->cm = nodeSelect->reverseChild->upConnection->cmTmp;
	}
	if(nodeSelect->revBulge13Child != NULL){
		nodeSelect->revBulge13Child->upConnection->cm = nodeSelect->revBulge13Child->upConnection->cmTmp;
	}
	if(nodeSelect->bulge14Child != NULL){
		nodeSelect->bulge14Child->upConnection->cm = nodeSelect->bulge14Child->upConnection->cmTmp;
	}

	int i,j;
	BRNode* node;

	for(i=0;i<nodeSelect->baseGroupC.size();i++){
		node = nodes[nodeSelect->baseGroupC[i]];
		node->cs1 = node->tmpCs1;
		for(j=0;j<node->baseAtomNum;j++){
			node->baseAtomCoords[j] = node->baseAtomCoordsTmp[j];
		}
	}

	for(i=0;i<nodeSelect->riboseGroupC.size();i++){
		node = nodes[nodeSelect->riboseGroupC[i]];
		node->rot = node->tmpRot;
		node->cs2 = node->tmpCs2;
		node->cs3 = node->tmpCs3;
		for(j=0;j<11;j++){
			node->riboAtomCoords[j] = node->riboAtomCoordsTmp[j];
		}
	}

	for(i=0;i<nodeSelect->phoGroupC.size();i++){
		node = nodes[nodeSelect->phoGroupC[i]];
		node->pho = node->phoTmp;
		node->phoLocal = node->phoLocalTmp;
	}

	int pi;

	for(i=0;i<seqLen;i++){
		allRotE[i] = tmpRotE[i];
		allRcE[i] = tmpRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseClashE[pi] = tmpBaseClashE[pi];
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
			allBaseRiboseE[pi] = tmpBaseRiboseE[pi];
			allBasePhoE[pi] = tmpBasePhoE[pi];
			allRiboseRiboseE[pi] = tmpRiboseRiboseE[pi];
			allRibosePhoE[pi] = tmpRibosePhoE[pi];
			allPhoPhoE[pi] = tmpPhoPhoE[pi];
		}
	}
}

pair<double,double> BRFoldingTree::ctMutEnergy(BRConnection* selectConnect, double breakCTWT, double connectWT,double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	int baseGroupANum = selectConnect->ctBaseGroupA.size();
	int baseGroupBNum = selectConnect->ctBaseGroupB.size();

	int riboseGroupANum = selectConnect->ctRiboseGroupA.size();
	int riboseGroupBNum = selectConnect->ctRiboseGroupB.size();

	int phoGroupANum = selectConnect->ctPhoGroupA.size();
	int phoGroupBNum = selectConnect->ctPhoGroupB.size();
	int phoGroupCNum = selectConnect->ctPhoGroupC.size();

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->ctPhoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
	}

	/*
	 * base-base energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->ctBaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupBNum;j++){
			indexY = selectConnect->ctBaseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	/*
	 * base-ribose energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->ctBaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->ctRiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->ctBaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = selectConnect->ctRiboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * base-pho energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->ctBaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->ctPhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->ctBaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->ctBaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->ctPhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->ctBaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->ctRiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->ctRiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	/*
	 * ribose-pho energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->ctRiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->ctPhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->ctRiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->ctRiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->ctPhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->ctRiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * pho-pho energy
	 */

	//A-B

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->ctPhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->ctPhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//A-C

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->ctPhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//B-C

	for(i=0;i<phoGroupBNum;i++){
		indexX = selectConnect->ctPhoGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//C-C

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->ctPhoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = selectConnect->ctPhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	double e1 = 0;
	double e2 = 0;
	double eb2 = 0;
	double e3 = 0;

	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += this->tmpRcE[i] - this->allRcE[i];
		else
			e2 += this->tmpRcE[i] - this->allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
			e3 += this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
			e3 += this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
			e1 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
			e3 += this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
		}
	}

	pair<double, double> p(e1+eb2+e2+e3, e1+eb2*breakCTWT+e2*connectWT+e3*clashWT);
	return p;
}

vector<double> BRFoldingTree::getF2MutEnergyDetail(BRConnection* selectConnect, double breakCTWT, double connectWT, double clashWT, double shift, bool verbose){

	vector<double> outList;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	int baseGroupANum = selectConnect->f2BaseGroupA.size();
	int baseGroupBNum = selectConnect->f2BaseGroupB.size();

	int riboseGroupANum = selectConnect->f2RiboseGroupA.size();
	int riboseGroupBNum = selectConnect->f2RiboseGroupB.size();
	int riboseGroupCNum = selectConnect->f2RiboseGroupC.size();



	int phoGroupANum = selectConnect->f2PhoGroupA.size();
	int phoGroupBNum = selectConnect->f2PhoGroupB.size();
	int phoGroupCNum = selectConnect->f2PhoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		this->tmpRotE[indexX] = nodes[indexX]->tmpRot->energy;
	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->f2PhoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
	}

	/*
	 * base-base energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupBNum;j++){
			indexY = selectConnect->f2BaseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	/*
	 * base-ribose energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->f2RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = selectConnect->f2RiboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * base-pho energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->f2RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//A-C

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//B-C

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//C-C

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	/*
	 * ribose-pho energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-B

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);

		}
	}

	/*
	 * pho-pho energy
	 */

	//A-B

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->f2PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//A-C

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->f2PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//B-C

	for(i=0;i<phoGroupBNum;i++){
		indexX = selectConnect->f2PhoGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//C-C

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->f2PhoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	double e1 = 0;
	double eb2 = 0;
	double e2 = 0;
	double e3 = 0;
	double e4 = 0;
	double e5 = 0;
	double e6 = 0;
	double e7 = 0;
	double e8 = 0;

	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += (this->tmpRcE[i] - this->allRcE[i])*breakCTWT;
		else
			e2 += (this->tmpRcE[i] - this->allRcE[i])*connectWT;
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e4 += (this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi])*clashWT;
			e5 += (this->tmpBasePhoE[pi] - this->allBasePhoE[pi])*clashWT;
			e7 += (this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi])*clashWT;
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += (this->tmpBaseClashE[pi] - this->allBaseClashE[pi])*clashWT;
			e6 += (this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi])*clashWT;
			e8 += (this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi])*clashWT;
		}
	}

	outList.push_back(e1);
	outList.push_back(eb2);
	outList.push_back(e2);
	outList.push_back(e3);
	outList.push_back(e4);
	outList.push_back(e5);
	outList.push_back(e6);
	outList.push_back(e7);
	outList.push_back(e8);

	return outList;
}

pair<double,double> BRFoldingTree::f2MutEnergy(BRConnection* selectConnect, double breakCTWT, double connectWT, double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	int baseGroupANum = selectConnect->f2BaseGroupA.size();
	int baseGroupBNum = selectConnect->f2BaseGroupB.size();

	int riboseGroupANum = selectConnect->f2RiboseGroupA.size();
	int riboseGroupBNum = selectConnect->f2RiboseGroupB.size();
	int riboseGroupCNum = selectConnect->f2RiboseGroupC.size();



	int phoGroupANum = selectConnect->f2PhoGroupA.size();
	int phoGroupBNum = selectConnect->f2PhoGroupB.size();
	int phoGroupCNum = selectConnect->f2PhoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		this->tmpRotE[indexX] = nodes[indexX]->tmpRot->energy;
	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->f2PhoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
	}

	/*
	 * base-base energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupBNum;j++){
			indexY = selectConnect->f2BaseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	/*
	 * base-ribose energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->f2RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = selectConnect->f2RiboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * base-pho energy
	 */

	//A-B

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = selectConnect->f2BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = selectConnect->f2BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = selectConnect->f2RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//A-C

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//B-C


	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//C-C

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<riboseGroupCNum;j++){
			indexY = selectConnect->f2RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	/*
	 * ribose-pho energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<riboseGroupANum;i++){
		indexX = selectConnect->f2RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<riboseGroupBNum;i++){
		indexX = selectConnect->f2RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = selectConnect->f2PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-B

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C

	for(i=0;i<riboseGroupCNum;i++){
		indexX = selectConnect->f2RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);

		}
	}

	/*
	 * pho-pho energy
	 */

	//A-B

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->f2PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = selectConnect->f2PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//A-C

	for(i=0;i<phoGroupANum;i++){
		indexX = selectConnect->f2PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//B-C

	for(i=0;i<phoGroupBNum;i++){
		indexX = selectConnect->f2PhoGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//C-C

	for(i=0;i<phoGroupCNum;i++){
		indexX = selectConnect->f2PhoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = selectConnect->f2PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	double e1 = 0;
	double e2 = 0;
	double e3 = 0;
	double eb2 = 0;

	if(verbose){

		double e;
		for(i=0;i<seqLen;i++){
			e = this->tmpRotE[i] - this->allRotE[i];
			if(e > 0.5){
				printf("rotE: %d %7.3f\n", i,e);
			}
			e = this->tmpRcE[i] - this->allRcE[i];
			if(e > 0.5){
				printf("conE: %d %7.3f\n", i, e);
			}
			for(j=0;j<seqLen;j++){
				pi = i*seqLen+j;
				e = this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
				if(e > 0.5){
					printf("baseRibose: %d %d %7.3f\n", i, j, e);
				}
				e = this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
				if(e > 0.5){
					printf("basePho: %d %d %7.3f\n", i, j, e);
				}
				e = this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
				if(e > 0.5){
					printf("ribosePho: %d %d %7.3f\n", i, j, e);
				}
			}
			for(j=i+1;j<seqLen;j++){
				pi = i*seqLen+j;
				e = this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
				if(e > 0.5){
					printf("baseClashE: %d %d %7.3f\n", i, j, e);
				}
				e = this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
				if(e > 0.5){
					printf("baseBaseE: %d %d %7.3f\n", i, j, e);
				}
				e = this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
				if(e > 0.5){
					printf("riboseRiboseE: %d %d %7.3f\n", i, j, e);
				}
				e = this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
				if(e > 0.5){
					printf("phoE phoE: %d %d %7.3f\n", i, j, e);
				}
			}
		}
	}


	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += this->tmpRcE[i] - this->allRcE[i];
		else
			e2 += this->tmpRcE[i] - this->allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
			e3 += this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
			e3 += this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
			e1 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
			e3 += this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
		}
	}

	pair<double, double> p(e1+eb2+e2+e3, e1+eb2*breakCTWT+e2*connectWT+e3*clashWT);
	return p;
}

pair<double,double> BRFoldingTree::f3MutEnergy(BRConnection* ct1,double breakCTWT, double connectWT,double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	int baseGroupANum = ct1->f3BaseGroupA.size();
	int baseGroupBNum = ct1->f3BaseGroupB.size();
	int baseGroupCNum = ct1->f3BaseGroupC.size();

	int riboseGroupANum = ct1->f3RiboseGroupA.size();
	int riboseGroupBNum = ct1->f3RiboseGroupB.size();
	int riboseGroupCNum = ct1->f3RiboseGroupC.size();

	int phoGroupANum = ct1->f3PhoGroupA.size();
	int phoGroupBNum = ct1->f3PhoGroupB.size();
	int phoGroupCNum = ct1->f3PhoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = ct1->f3RiboseGroupC[i];
		this->tmpRotE[indexX] = nodes[indexX]->tmpRot->energy;
	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = ct1->f3PhoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
	}


	if(verbose){
		cout << "base groupA: ";
		for(i=0;i<baseGroupANum;i++){
			cout << ct1->f3BaseGroupA[i] << " ";
		}
		cout << "base groupB: ";
		for(i=0;i<baseGroupBNum;i++){
			cout << ct1->f3BaseGroupB[i] << " ";
		}
		cout << "base groupC: ";
		for(i=0;i<baseGroupCNum;i++){
			cout << ct1->f3BaseGroupC[i] << " ";
		}
		cout << endl;

		cout << "ribo groupA: ";
		for(i=0;i<riboseGroupANum;i++){
			cout << ct1->f3RiboseGroupA[i] << " ";
		}
		cout << "ribo groupB: ";
		for(i=0;i<riboseGroupBNum;i++){
			cout << ct1->f3RiboseGroupB[i] << " ";
		}
		cout << "ribo groupC: ";
		for(i=0;i<riboseGroupCNum;i++){
			cout << ct1->f3RiboseGroupC[i] << " ";
		}
		cout << endl;

		cout << "pho groupA: ";
		for(i=0;i<phoGroupANum;i++){
			cout << ct1->f3PhoGroupA[i] << " ";
		}
		cout << "pho groupB: ";
		for(i=0;i<phoGroupBNum;i++){
			cout << ct1->f3PhoGroupB[i] << " ";
		}
		cout << "pho groupC: ";
		for(i=0;i<phoGroupCNum;i++){
			cout << ct1->f3PhoGroupC[i] << " ";
		}
		cout << endl;
	}


	/*
	 * base-base energy
	 */

	//A-B
	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupBNum;j++){
			indexY = ct1->f3BaseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	//A-C
	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupCNum;j++){
			indexY = ct1->f3BaseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift,  verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}


	//B-C
	for(i=0;i<baseGroupBNum;i++){
		indexX = ct1->f3BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupCNum;j++){
			indexY = ct1->f3BaseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	/*
	 * base-ribose energy
	 */

	//A-B
	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = ct1->f3RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C

	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);

		}
	}

	//B-A

	for(i=0;i<baseGroupBNum;i++){
		indexX = ct1->f3BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = ct1->f3RiboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C

	for(i=0;i<baseGroupBNum;i++){
		indexX = ct1->f3BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A

	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = ct1->f3RiboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-B

	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = ct1->f3RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C

	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * base-pho energy
	 */

	//A-B
	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = ct1->f3PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C
	for(i=0;i<baseGroupANum;i++){
		indexX = ct1->f3BaseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A
	for(i=0;i<baseGroupBNum;i++){
		indexX = ct1->f3BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = ct1->f3PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C
	for(i=0;i<baseGroupBNum;i++){
		indexX = ct1->f3BaseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A
	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = ct1->f3PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-B
	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = ct1->f3PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C

	for(i=0;i<baseGroupCNum;i++){
		indexX = ct1->f3BaseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * ribose-ribose energy
	 */

	//A-B

	for(i=0;i<riboseGroupANum;i++){
		indexX = ct1->f3RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupBNum;j++){
			indexY = ct1->f3RiboseGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = ct1->f3RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}
	//B-C
	for(i=0;i<riboseGroupBNum;i++){
		indexX = ct1->f3RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}
	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = ct1->f3RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<riboseGroupCNum;j++){
			indexY = ct1->f3RiboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;

			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}
	/*
	 * ribose-pho energy
	 */

	//A-B
	for(i=0;i<riboseGroupANum;i++){
		indexX = ct1->f3RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = ct1->f3PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//A-C
	for(i=0;i<riboseGroupANum;i++){
		indexX = ct1->f3RiboseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-A
	for(i=0;i<riboseGroupBNum;i++){
		indexX = ct1->f3RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = ct1->f3PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//B-C
	for(i=0;i<riboseGroupBNum;i++){
		indexX = ct1->f3RiboseGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A
	for(i=0;i<riboseGroupCNum;i++){
		indexX = ct1->f3RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = ct1->f3PhoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-B
	for(i=0;i<riboseGroupCNum;i++){
		indexX = ct1->f3RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = ct1->f3PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = ct1->f3RiboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;

			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * pho-pho energy
	 */

	//A-B
	for(i=0;i<phoGroupANum;i++){
		indexX = ct1->f3PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupBNum;j++){
			indexY = ct1->f3PhoGroupB[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;

			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//A-C
	for(i=0;i<phoGroupANum;i++){
		indexX = ct1->f3PhoGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;

			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//B-C
	for(i=0;i<phoGroupBNum;i++){
		indexX = ct1->f3PhoGroupB[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;

			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	//C-C

	for(i=0;i<phoGroupCNum;i++){
		indexX = ct1->f3PhoGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<phoGroupCNum;j++){
			indexY = ct1->f3PhoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen + indexY;
			pj = indexY*seqLen + indexX;

			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}


	double e1 = 0;
	double eb2 = 0;
	double e2 = 0;
	double e3 = 0;

	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += this->tmpRcE[i] - this->allRcE[i];
		else
			e2 += this->tmpRcE[i] - this->allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
			e3 += this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
			e3 += this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
			e1 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
			e3 += this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
		}
	}

	pair<double, double> p(e1+eb2+e2+e3, e1+eb2*breakCTWT+e2*connectWT+e3*clashWT);
	return p;
}

pair<double,double> BRFoldingTree::singleBaseMutEnergy(BRNode* node,double breakCTWT, double connectWT,double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	int baseGroupANum = node->baseGroupA.size();
	int baseGroupCNum = node->baseGroupC.size();

	int riboseGroupANum = node->riboseGroupA.size();
	int riboseGroupCNum = node->riboseGroupC.size();

	int phoGroupANum = node->phoGroupA.size();
	int phoGroupCNum = node->phoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		this->tmpRotE[indexX] = nodes[indexX]->tmpRot->energy;
	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
	}

	/*
	 * base-base energy
	 */

	//A-C
	for(i=0;i<baseGroupANum;i++){
		indexX = node->baseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<baseGroupCNum;j++){
			indexY = node->baseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpBaseBaseE[pi] = getBaseBaseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpBaseClashE[pi] = baseBaseClashTmp(nodeX, nodeY, sepTable[pi], et, shift, verbose);
			this->tmpBaseBaseE[pj] = this->tmpBaseBaseE[pi];
			this->tmpBaseClashE[pj] = this->tmpBaseClashE[pi];
		}
	}

	/*
	 * base-ribose energy
	 */

	//A-C
	for(i=0;i<baseGroupANum;i++){
		indexX = node->baseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A
	for(i=0;i<baseGroupCNum;i++){
		indexX = node->baseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupANum;j++){
			indexY = node->riboseGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C
	for(i=0;i<baseGroupCNum;i++){
		indexX = node->baseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	/*
	 * base-pho energy
	 */

	//A-C
	for(i=0;i<baseGroupANum;i++){
		indexX = node->baseGroupA[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-A
	for(i=0;i<baseGroupCNum;i++){
		indexX = node->baseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupANum;j++){
			indexY = node->phoGroupA[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
		}
	}

	//C-C
	for(i=0;i<baseGroupCNum;i++){
		indexX = node->baseGroupC[i];
		nodeX = nodes[indexX];
		for(j=0;j<phoGroupCNum;j++){
			indexY = node->phoGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
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
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	double e1 = 0;
	double eb2 = 0;
	double e2 = 0;
	double e3 = 0;

	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += this->tmpRcE[i] - this->allRcE[i];
		else
			e2 += this->tmpRcE[i] - this->allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
			e3 += this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
			e3 += this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
			e1 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
			e3 += this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
		}
	}

	pair<double, double> p(e1+eb2+e2+e3, e1+eb2*breakCTWT+e2*connectWT+e3*clashWT);
	return p;
}

pair<double,double> BRFoldingTree::baseRotamerMutEnergy(BRNode* node, double breakCTWT, double connectWT,double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int indexX, indexY;
	BRNode *nodeX, *nodeY;

	//int baseGroupANum = node->baseGroupA.size();


	int riboseGroupANum = node->riboseGroupA.size();
	int riboseGroupCNum = node->riboseGroupC.size();

	int phoGroupANum = node->phoGroupA.size();
	int phoGroupCNum = node->phoGroupC.size();

	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		this->tmpRotE[indexX] = nodes[indexX]->tmpRot->energy;
	}

	for(i=0;i<phoGroupCNum;i++){
		indexX = node->phoGroupC[i];
		this->tmpRcE[indexX] = nodes[indexX]->phoLocalTmp.e;
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
			this->tmpBaseRiboseE[pi] = getBaseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpBasePhoE[pi] = getBasePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
		}
	}

	//C-C
	for(i=0;i<riboseGroupCNum;i++){
		indexX = node->riboseGroupC[i];
		nodeX = nodes[indexX];
		for(j=i+1;j<riboseGroupCNum;j++){
			indexY = node->riboseGroupC[j];
			nodeY = nodes[indexY];
			pi = indexX*seqLen+indexY;
			pj = indexY*seqLen+indexX;
			this->tmpRiboseRiboseE[pi] = getRiboseRiboseEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpRiboseRiboseE[pj] = this->tmpRiboseRiboseE[pi];
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpRibosePhoE[pi] = getRibosePhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
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
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
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
			this->tmpPhoPhoE[pi] = getPhoPhoEnergyTmp(nodeX, nodeY, sepTable[pi], et, verbose);
			this->tmpPhoPhoE[pj] = this->tmpPhoPhoE[pi];
		}
	}

	double e1 = 0;
	double eb2 = 0;
	double e2 = 0;
	double e3 = 0;

	for(i=0;i<seqLen;i++){
		e1 += this->tmpRotE[i] - this->allRotE[i];
		if(chainBreakPoints[i])
			eb2 += this->tmpRcE[i] - this->allRcE[i];
		else
			e2 += this->tmpRcE[i] - this->allRcE[i];
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseRiboseE[pi] - this->allBaseRiboseE[pi];
			e3 += this->tmpBasePhoE[pi] - this->allBasePhoE[pi];
			e3 += this->tmpRibosePhoE[pi] - this->allRibosePhoE[pi];
		}
		for(j=i+1;j<seqLen;j++){
			pi = i*seqLen+j;
			e3 += this->tmpBaseClashE[pi] - this->allBaseClashE[pi];
			e1 += this->tmpBaseBaseE[pi] - this->allBaseBaseE[pi];
			e3 += this->tmpRiboseRiboseE[pi] - this->allRiboseRiboseE[pi];
			e3 += this->tmpPhoPhoE[pi] - this->allPhoPhoE[pi];
		}
	}

	pair<double, double> p(e1+eb2+e2+e3, e1+eb2*breakCTWT+e2*connectWT+e3*clashWT);
	return p;
}

double BRFoldingTree::totalEnergy(double breakCTWT, double connectWT, double clashWT, double shift, bool verbose){
	double tot = 0;
	int i,j;

	//base-base, ribose-ribose, pho-pho
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			tot += getBaseBaseEnergy(nodeA, nodeB, sep, et, verbose);
			tot += baseBaseClash(nodeA, nodeB, sep, et, shift, verbose)*clashWT;
			tot += getRiboseRiboseEnergy(nodeA, nodeB, sep, et, verbose)*clashWT;
			tot += getPhoPhoEnergy(nodeA, nodeB, sep, et, verbose)*clashWT;
		}
	}

	//base-ribose, base-pho, ribose-pho

	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			tot += getBaseRiboseEnergy(nodeA, nodeB, sep, et, verbose)*clashWT;
			tot += getBasePhoEnergy(nodeA, nodeB, sep,  et, verbose)*clashWT;
			tot += getRibosePhoEnergy(nodeA, nodeB, sep, et, verbose)*clashWT;
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
		if(chainBreakPoints[i])
			tot += nodes[i]->pho.e * breakCTWT;
		else
			tot += nodes[i]->pho.e * connectWT;
	}
	return tot;
}

void BRFoldingTree::printDetailEnergy(){
	double e = 0;

	for(int i=0;i<seqLen;i++){
		e += nodes[i]->rot->energy;
		e += nodes[i]->pho.e;
		printf("Node: %2d RotamerE: %9.3f PhoE: %9.3f\n", i, nodes[i]->rot->energy, nodes[i]->pho.e);
		if(i < seqLen-1){
			double dis = nodes[i]->pho.tList[1].distance(nodes[i+1]->riboAtomCoords[7]);
			cout << dis << endl;
		}
	}
	bool verbose = false;

	double shift = 0.0;

	for(int i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(int j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int pi = i*seqLen+j;
			int pj = j*seqLen+i;
			int sepi = sepTable[pi];
			int sepj = sepTable[pj];
			double tot = 0;
			double bb = getBaseBaseEnergy(nodeA, nodeB, sepi, et, verbose) + baseBaseClash(nodeA, nodeB, sepi, et, shift, verbose);

			double bbClash = baseBaseClash(nodeA, nodeB, sepi, et, shift, verbose);
			double bbTmpClash = baseBaseClashTmp(nodeA, nodeB, sepi, et, shift, verbose);
			double br = getBaseRiboseEnergy(nodeA, nodeB, sepi, et, verbose);
			double rb = getBaseRiboseEnergy(nodeB, nodeA, sepj, et, verbose);
			double bp = getBasePhoEnergy(nodeA, nodeB, sepi, et, verbose);
			double pb = getBasePhoEnergy(nodeB, nodeA, sepj, et, verbose);
			double rr = getRiboseRiboseEnergy(nodeA, nodeB, sepi, et, verbose);
			double rp = getRibosePhoEnergy(nodeA, nodeB, sepi, et, verbose);
			double pr = getRibosePhoEnergy(nodeB, nodeA, sepj, et, verbose);
			double pp = getPhoPhoEnergy(nodeA, nodeB, sepi, et, verbose);
			tot = bb + br + rb + bp + pb + rr + rp + pr + pp;
			printf("Pair: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
			e += tot;
		}
	}

	printf("total energy: %9.3f\n",  e);
	cout << endl;
}

void BRFoldingTree::printDetailEnergy(ofstream& of){
	double e = 0;

	char ss[200];
	for(int i=0;i<seqLen;i++){
		e += nodes[i]->rot->energy;
		e += nodes[i]->pho.e;
		sprintf(ss, "Node: %2d RotamerE: %9.3f PhoE: %9.3f\n", i, nodes[i]->rot->energy, nodes[i]->pho.e);
		of << string(ss);
		//printf("Node: %2d RotamerE: %9.3f PhoE: %9.3f\n", i, nodes[i]->rot->energy, nodes[i]->pho.e);
	}
	bool verbose = false;

	double shift = 0.0;

	for(int i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(int j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int pi = i*seqLen+j;
			int pj = j*seqLen+i;
			int sepi = sepTable[pi];
			int sepj = sepTable[pj];
			double tot = 0;

			double bb = getBaseBaseEnergy(nodeA, nodeB, sepi, et, verbose) + baseBaseClash(nodeA, nodeB, sepi, et, shift, verbose);
			double br = getBaseRiboseEnergy(nodeA, nodeB, sepi, et, verbose);
			double rb = getBaseRiboseEnergy(nodeB, nodeA, sepj, et, verbose);
			double bp = getBasePhoEnergy(nodeA, nodeB, sepi, et, verbose);
			double pb = getBasePhoEnergy(nodeB, nodeA, sepj, et, verbose);
			double rr = getRiboseRiboseEnergy(nodeA, nodeB, sepi, et, verbose);
			double rp = getRibosePhoEnergy(nodeA, nodeB, sepi, et, verbose);
			double pr = getRibosePhoEnergy(nodeB, nodeA, sepj, et, verbose);
			double pp = getPhoPhoEnergy(nodeA, nodeB, sepi, et, verbose);

			sprintf(ss, "Pair: %2d-%-2d bb %7.3f br %7.3f rb %7.3f bp %7.3f pb %7.3f rr %7.3f rp %7.3f pr %7.3f pp %7.3f\n",i,j, bb, br, rb, bp, pb, rr, rp, pr, pp);
			of << string(ss);
			e += tot;
		}
	}
	//printf("total energy: %9.3f\n",  e);
}

void BRFoldingTree::trackCoordinateChangeCt(BRConnection* ct) {
	vector<int> baseListA;
	vector<int> baseListB;
	vector<int> riboseListA;
	vector<int> riboseListB;
	vector<int> phoListA;
	vector<int> phoListB;
	vector<int> phoListC;

	vector<vector<int>> listList;


	vector<string> nameList;
	nameList.push_back("baseListA");
	nameList.push_back("baseListB");
	nameList.push_back("riboseListA");
	nameList.push_back("riboseListB");
	nameList.push_back("phoListA");
	nameList.push_back("phoListB");
	nameList.push_back("phoListC");

	vector<vector<int>> ctListList;
	ctListList.push_back(ct->ctBaseGroupA);
	ctListList.push_back(ct->ctBaseGroupB);
	ctListList.push_back(ct->ctRiboseGroupA);
	ctListList.push_back(ct->ctRiboseGroupB);
	ctListList.push_back(ct->ctPhoGroupA);
	ctListList.push_back(ct->ctPhoGroupB);
	ctListList.push_back(ct->ctPhoGroupC);

	cout << "base listA: ";
	for(int i=0;i<ct->ctBaseGroupA.size();i++){
		cout << " " << ct->ctBaseGroupA[i];
	}
	cout << endl;

	cout << "base listB: ";
	for(int i=0;i<ct->ctBaseGroupB.size();i++){
		cout << " " << ct->ctBaseGroupB[i];
	}
	cout << endl;

	cout << "ribose listA: ";
	for(int i=0;i<ct->ctRiboseGroupA.size();i++){
		cout << " " << ct->ctRiboseGroupA[i];
	}
	cout << endl;
	cout << "ribose listB: ";
	for(int i=0;i<ct->ctRiboseGroupB.size();i++){
		cout << " " << ct->ctRiboseGroupB[i];
	}

	cout << endl;
	cout << "pho listA: ";
	for(int i=0;i<ct->ctPhoGroupA.size();i++){
		cout << " " << ct->ctPhoGroupA[i];
	}
	cout << endl;

	cout << "pho listB: ";
	for(int i=0;i<ct->ctPhoGroupB.size();i++){
		cout << " " << ct->ctPhoGroupB[i];
	}
	cout << endl;

	cout << "pho listC: ";
	for(int i=0;i<ct->ctPhoGroupC.size();i++){
		cout << " " << ct->ctPhoGroupC[i];
	}
	cout << endl;


	for(int i=0;i<seqLen;i++){
		if(nodes[i]->baseConsistent())
			baseListA.push_back(i);
		else
			baseListB.push_back(i);
	}

	for(int i=0;i<seqLen;i++){
		if(nodes[i]->riboConsistent())
			riboseListA.push_back(i);
		else if(nodes[i]->rotamerConsistent())
			riboseListB.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i] && nodes[i]->phoConsistent())
			phoListA.push_back(i);
		else if(connectToDownstream[i] && nodes[i]->phoLocalConsistent())
			phoListB.push_back(i);
		else if(connectToDownstream[i])
			phoListC.push_back(i);
	}

	listList.push_back(baseListA);
	listList.push_back(baseListB);
	listList.push_back(riboseListA);
	listList.push_back(riboseListB);
	listList.push_back(phoListA);
	listList.push_back(phoListB);
	listList.push_back(phoListC);

	for(int k=0;k<nameList.size();k++){
		cout << "k=" << k << endl;
		bool changeConsistent = true;
		if(listList[k].size() != ctListList[k].size())
			changeConsistent = false;
		else {
			for(int i=0;i<listList[k].size();i++){
				if(listList[k][i] != ctListList[k][i])
					changeConsistent = false;
			}
		}

		if(!changeConsistent){
			cout << "connect: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << endl;
			cout << nameList[k] << " not consistent" << endl;
			for(int i=0;i<listList[k].size();i++) {
				cout << listList[k][i] << " ";
			}
			cout << endl;

			for(int i=0;i<ctListList[k].size();i++) {
				cout << ctListList[k][i] << " ";
			}
			cout << endl;
		}
	}

	cout << "track finished" << endl;

}

void BRFoldingTree::trackCoordinateChangeF2(BRConnection* ct){
	/*
	 * base ribose change
	 */
	cout << "COOR changed: ";

	vector<int> baseListA;
	vector<int> baseListB;
	vector<int> riboseListA;
	vector<int> riboseListB;
	vector<int> riboseListC;
	vector<int> phoListA;
	vector<int> phoListB;
	vector<int> phoListC;

	vector<vector<int>> listList;


	vector<string> nameList;
	nameList.push_back("baseListA");
	nameList.push_back("baseListB");
	nameList.push_back("riboseListA");
	nameList.push_back("riboseListB");
	nameList.push_back("riboseListC");
	nameList.push_back("phoListA");
	nameList.push_back("phoListB");
	nameList.push_back("phoListC");

	vector<vector<int>> ctListList;
	ctListList.push_back(ct->f2BaseGroupA);
	ctListList.push_back(ct->f2BaseGroupB);
	ctListList.push_back(ct->f2RiboseGroupA);
	ctListList.push_back(ct->f2RiboseGroupB);
	ctListList.push_back(ct->f2RiboseGroupC);
	ctListList.push_back(ct->f2PhoGroupA);
	ctListList.push_back(ct->f2PhoGroupB);
	ctListList.push_back(ct->f2PhoGroupC);

	cout << "base listA: ";
	for(int i=0;i<ct->f2BaseGroupA.size();i++){
		cout << " " << ct->f2BaseGroupA[i];
	}
	cout << endl;

	cout << "base listB: ";
	for(int i=0;i<ct->f2BaseGroupB.size();i++){
		cout << " " << ct->f2BaseGroupB[i];
	}
	cout << endl;

	cout << "ribose listA: ";
	for(int i=0;i<ct->f2RiboseGroupA.size();i++){
		cout << " " << ct->f2RiboseGroupA[i];
	}
	cout << endl;
	cout << "ribose listB: ";
	for(int i=0;i<ct->f2RiboseGroupB.size();i++){
		cout << " " << ct->f2RiboseGroupB[i];
	}
	cout << endl;
	cout << "ribose listC: ";
	for(int i=0;i<ct->f2RiboseGroupC.size();i++){
		cout << " " << ct->f2RiboseGroupC[i];
	}

	cout << endl;
	cout << "pho listA: ";
	for(int i=0;i<ct->f2PhoGroupA.size();i++){
		cout << " " << ct->f2PhoGroupA[i];
	}
	cout << endl;

	cout << "pho listB: ";
	for(int i=0;i<ct->f2PhoGroupB.size();i++){
		cout << " " << ct->f2PhoGroupB[i];
	}
	cout << endl;

	cout << "pho listC: ";
	for(int i=0;i<ct->f2PhoGroupC.size();i++){
		cout << " " << ct->f2PhoGroupC[i];
	}
	cout << endl;



	for(int i=0;i<seqLen;i++){
		if(nodes[i]->baseConsistent())
			baseListA.push_back(i);
		else
			baseListB.push_back(i);
	}

	for(int i=0;i<seqLen;i++){
		if(nodes[i]->riboConsistent())
			riboseListA.push_back(i);
		else if(nodes[i]->rotamerConsistent())
			riboseListB.push_back(i);
		else
			riboseListC.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i] && nodes[i]->phoConsistent())
			phoListA.push_back(i);
		else if(connectToDownstream[i] && nodes[i]->phoLocalConsistent())
			phoListB.push_back(i);
		else if(connectToDownstream[i])
			phoListC.push_back(i);
	}

	listList.push_back(baseListA);
	listList.push_back(baseListB);
	listList.push_back(riboseListA);
	listList.push_back(riboseListB);
	listList.push_back(riboseListC);
	listList.push_back(phoListA);
	listList.push_back(phoListB);
	listList.push_back(phoListC);

	for(int k=0;k<8;k++){
		bool changeConsistent = true;
		if(listList[k].size() != ctListList[k].size())
			changeConsistent = false;
		else {
			for(int i=0;i<listList[k].size();i++){
				if(listList[k][i] != ctListList[k][i])
					changeConsistent = false;
			}
		}

		if(!changeConsistent){
			cout << "connect: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << endl;
			cout << nameList[k] << " not consistent" << endl;
			for(int i=0;i<listList[k].size();i++) {
				cout << listList[k][i] << " ";
			}
			cout << endl;

			for(int i=0;i<ctListList[k].size();i++) {
				cout << ctListList[k][i] << " ";
			}
			cout << endl;
		}
	}


}

void BRFoldingTree::trackCoordinateChangeF3(BRConnection* ct){

	vector<int> baseListA;
	vector<int> baseListB;
	vector<int> baseListC;
	vector<int> riboseListA;
	vector<int> riboseListB;
	vector<int> riboseListC;
	vector<int> phoListA;
	vector<int> phoListB;
	vector<int> phoListC;

	vector<vector<int>> listList;


	vector<string> nameList;
	nameList.push_back("baseListA");

	nameList.push_back("riboseListA");
	nameList.push_back("riboseListB");
	nameList.push_back("riboseListC");
	nameList.push_back("phoListA");
	nameList.push_back("phoListB");
	nameList.push_back("phoListC");

	vector<vector<int>> ctListList;
	ctListList.push_back(ct->f3BaseGroupA);

	ctListList.push_back(ct->f3RiboseGroupA);
	ctListList.push_back(ct->f3RiboseGroupB);
	ctListList.push_back(ct->f3RiboseGroupC);
	ctListList.push_back(ct->f3PhoGroupA);
	ctListList.push_back(ct->f3PhoGroupB);
	ctListList.push_back(ct->f3PhoGroupC);

	for(int i=0;i<seqLen;i++){
		if(nodes[i]->baseConsistent())
			baseListA.push_back(i);
		else
			baseListB.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(nodes[i]->riboConsistent())
			riboseListA.push_back(i);
		else if(nodes[i]->rotamerConsistent())
			riboseListB.push_back(i);
		else
			riboseListC.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i] && nodes[i]->phoConsistent())
			phoListA.push_back(i);
		else if(connectToDownstream[i] && nodes[i]->phoLocalConsistent())
			phoListB.push_back(i);
		else if(connectToDownstream[i])
			phoListC.push_back(i);
	}

	listList.push_back(baseListA);
	listList.push_back(riboseListA);
	listList.push_back(riboseListB);
	listList.push_back(riboseListC);
	listList.push_back(phoListA);
	listList.push_back(phoListB);
	listList.push_back(phoListC);

	for(int k=0;k<nameList.size();k++){
		bool changeConsistent = true;
		if(listList[k].size() != ctListList[k].size())
			changeConsistent = false;
		else {
			for(int i=0;i<listList[k].size();i++){
				if(listList[k][i] != ctListList[k][i])
					changeConsistent = false;
			}
		}

		if(!changeConsistent){
			cout << "connect: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << "->" << ct->childNode->midChild->seqID << endl;
			cout << nameList[k] << " not consistent" << endl;
			cout << "changed: ";
			for(int i=0;i<listList[k].size();i++) {
				cout << listList[k][i] << " ";
			}
			cout << endl;
			cout << "should be change: ";
			for(int i=0;i<ctListList[k].size();i++) {
				cout << ctListList[k][i] << " ";
			}
			cout << endl;
		}
	}
}

void BRFoldingTree::trackCoordinateChangeSingleBase(BRNode* node){
	vector<int> baseListA;
	vector<int> baseListC;
	vector<int> riboseListA;
	vector<int> riboseListC;
	vector<int> phoListA;
	vector<int> phoListC;

	vector<vector<int>> listList;


	vector<string> nameList;
	nameList.push_back("baseListA");
	nameList.push_back("baseListC");
	nameList.push_back("riboseListA");
	nameList.push_back("riboseListC");
	nameList.push_back("phoListA");
	nameList.push_back("phoListC");

	vector<vector<int>> ctListList;
	ctListList.push_back(node->baseGroupA);
	ctListList.push_back(node->baseGroupC);
	ctListList.push_back(node->riboseGroupA);
	ctListList.push_back(node->riboseGroupC);
	ctListList.push_back(node->phoGroupA);
	ctListList.push_back(node->phoGroupC);


	for(int i=0;i<seqLen;i++){
		if(nodes[i]->baseConsistent())
			baseListA.push_back(i);
		else
			baseListC.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(nodes[i]->riboConsistent())
			riboseListA.push_back(i);
		else
			riboseListC.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i] && nodes[i]->phoConsistent())
			phoListA.push_back(i);
		else if(connectToDownstream[i])
			phoListC.push_back(i);
	}

	listList.push_back(baseListA);
	listList.push_back(baseListC);
	listList.push_back(riboseListA);
	listList.push_back(riboseListC);
	listList.push_back(phoListA);
	listList.push_back(phoListC);

	for(int k=0;k<nameList.size();k++){
		bool changeConsistent = true;
		if(listList[k].size() != ctListList[k].size())
			changeConsistent = false;
		else {
			for(int i=0;i<listList[k].size();i++){
				if(listList[k][i] != ctListList[k][i])
					changeConsistent = false;
			}
		}


		if(!changeConsistent){
			cout << "node: " << node->seqID << endl;
			cout << nameList[k] << " not consistent" << endl;
			cout << "changed: ";
			for(int i=0;i<listList[k].size();i++) {
				cout << listList[k][i] << " ";
			}
			cout << endl;
			cout << "should be change: ";
			for(int i=0;i<ctListList[k].size();i++) {
				cout << ctListList[k][i] << " ";
			}
			cout << endl;
		}
	}
}

void BRFoldingTree::trackCoordinateChangeRotamer(BRNode* node){

	vector<int> riboseListA;
	vector<int> riboseListC;
	vector<int> phoListA;
	vector<int> phoListC;

	vector<vector<int>> listList;


	vector<string> nameList;

	nameList.push_back("riboseListA");
	nameList.push_back("riboseListC");
	nameList.push_back("phoListA");
	nameList.push_back("phoListC");

	vector<vector<int>> ctListList;

	ctListList.push_back(node->riboseGroupA);
	ctListList.push_back(node->riboseGroupC);
	ctListList.push_back(node->phoGroupA);
	ctListList.push_back(node->phoGroupC);



	for(int i=0;i<seqLen;i++){
		if(nodes[i]->riboConsistent())
			riboseListA.push_back(i);
		else
			riboseListC.push_back(i);
	}
	for(int i=0;i<seqLen;i++){
		if(connectToDownstream[i] && nodes[i]->phoConsistent())
			phoListA.push_back(i);
		else if(connectToDownstream[i])
			phoListC.push_back(i);
	}


	listList.push_back(riboseListA);
	listList.push_back(riboseListC);
	listList.push_back(phoListA);
	listList.push_back(phoListC);

	for(int k=0;k<nameList.size();k++){
		bool changeConsistent = true;
		if(listList[k].size() != ctListList[k].size())
			changeConsistent = false;
		else {
			for(int i=0;i<listList[k].size();i++){
				if(listList[k][i] != ctListList[k][i])
					changeConsistent = false;
			}
		}


		if(!changeConsistent){
			cout << "node: " << node->seqID << endl;
			cout << nameList[k] << " not consistent" << endl;
			cout << "changed: ";
			for(int i=0;i<listList[k].size();i++) {
				cout << listList[k][i] << " ";
			}
			cout << endl;
			cout << "should be change: ";
			for(int i=0;i<ctListList[k].size();i++) {
				cout << ctListList[k][i] << " ";
			}
			cout << endl;
		}
	}
}


void BRFoldingTree::checkConnection(){
	for(int i=0;i<flexibleConnectionList.size();i++){
		BRConnection* ct = flexibleConnectionList[i];
		LocalFrame cs = ct->fatherNode->cs1 + ct->cm;
		if(!cs.equalTo(ct->childNode->cs1)){

			cout << "connection error: " << ct->fatherNode->seqID << "-" << ct->childNode->seqID << endl;
			cs.print();
			ct->childNode->cs1.print();
		}
	}

	for(int i=0;i<flexibleConnectionList.size();i++){
		BRConnection* ct = flexibleConnectionList[i];
		LocalFrame cs = ct->fatherNode->tmpCs1 + ct->cmTmp;
		if(!cs.equalTo(ct->childNode->tmpCs1)){
			cout << "tmp connection error: " << ct->fatherNode->seqID << "-" << ct->childNode->seqID << endl;
			cs.print();
			ct->childNode->tmpCs1.print();
		}
	}

}

void BRFoldingTree::checkNode(){
	cout << "check nodes: " << endl;
	for(int i=0;i<seqLen;i++){
		BRNode* node = nodes[i];
		LocalFrame cs = node->cs1;
		XYZ t = local2global(cs, node->atomCoordLocal[1]);
		if(t.distance(node->baseAtomCoords[1]) > 0.001){
			cs.print();
			cout << node->atomCoordLocal[1].toString() << endl;
			cout << t.toString() << endl;
			cout << node->baseAtomCoords[1].toString() << endl;
			cout << "node cs error: " << i << endl;
		}
	}

	for(int i=0;i<seqLen;i++){
		BRNode* node = nodes[i];
		LocalFrame cs = node->tmpCs1;
		XYZ t = local2global(cs, node->atomCoordLocal[1]);
		if(t.distance(node->baseAtomCoordsTmp[1]) > 0.001){
			cout << "node tmp cs error: " << i << endl;
		}
	}
}

void BRFoldingTree::checkRibose(){
	cout << "check ribose: " << endl;
	for(int i=0;i<seqLen;i++){
		BRNode* node = nodes[i];
		LocalFrame cs = node->cs1;
		LocalFrame cs2 = node->cs1 + node->rot->mv12;
		LocalFrame cs3 = node->cs1 + node->rot->mv13;
		XYZ t = local2global(cs, node->rot->tList1[1]);
		if(t.distance(node->riboAtomCoords[1]) > 0.001){
			cout << t.toString() << endl;
			cout << node->riboAtomCoords[1].toString() << endl;
			cout << "ribose cs error: " << i << endl;
		}
		if(!cs2.equalTo(node->cs2)) {
			cout << "cs2 error " << i << endl;
		}
		if(!cs3.equalTo(node->cs3)) {
			cout << "cs3 error " << i << endl;
		}
	}

	for(int i=0;i<seqLen;i++){
		BRNode* node = nodes[i];
		LocalFrame cs = node->tmpCs1;
		LocalFrame cs2 = node->tmpCs1 + node->tmpRot->mv12;
		LocalFrame cs3 = node->tmpCs1 + node->tmpRot->mv13;
		XYZ t = local2global(cs, node->tmpRot->tList1[1]);
		if(t.distance(node->riboAtomCoordsTmp[1]) > 0.001){
			cout << t.toString() << endl;
			cout << node->riboAtomCoordsTmp[1].toString() << endl;
			cout << "ribose tmp cs error: " << i << endl;
		}
		if(!cs2.equalTo(node->tmpCs2)) {
			cout << "tmp cs2 error " << i << endl;
		}
		if(!cs3.equalTo(node->tmpCs3)) {
			cout << "tmp cs3 error " << i << endl;
		}
	}
}

void BRFoldingTree::checkPho(){
	cout << "check pho: " << endl;
	for(int i=0;i<seqLen-1;i++){
		BRNode* node = nodes[i];
		if(!node->connectToNeighbor) continue;
		LocalFrame cs = node->cs2;
		XYZ t = local2global(cs, node->phoLocal.p);
		if(t.distance(node->pho.tList[0]) > 0.001){
			cout << "pho cs error: " << i << endl;
		}
	}

	for(int i=0;i<seqLen-1;i++){
		BRNode* node = nodes[i];
		if(!node->connectToNeighbor) continue;
		LocalFrame cs = node->tmpCs2;
		XYZ t = local2global(cs, node->phoLocalTmp.p);
		if(t.distance(node->phoTmp.tList[0]) > 0.001){
			cout << "pho tmp cs error: " << i << endl;
		}
	}
}

void BRFoldingTree::energyChange(){

	int i,j;
	double e1, e2;
	for(i=0;i<seqLen;i++){
		for(j=i+1;j<seqLen;j++){
			if(abs(allBaseBaseE[i*seqLen+j] - tmpBaseBaseE[i*seqLen+j])>0.001){
				printf("EnergyChange: base %2d - base %2d\n", i, j);
			}
			if(abs(allBaseClashE[i*seqLen+j] - tmpBaseClashE[i*seqLen+j])>0.001){
				printf("ClashEnergyChange: base %2d - base %2d\n", i, j);
			}
			if(abs(allRiboseRiboseE[i*seqLen+j] - tmpRiboseRiboseE[i*seqLen+j])>0.001){
				printf("EnergyChange: ribo %2d - ribo %2d\n", i, j);
			}
			if(abs(allPhoPhoE[i*seqLen+j] - tmpPhoPhoE[i*seqLen+j])>0.001){
				printf("EnergyChange: pho %2d - pho %2d\n", i, j);
			}
		}
	}

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			if(abs(allBaseRiboseE[i*seqLen+j] - tmpBaseRiboseE[i*seqLen+j]) > 0.001){
				printf("EnergyChange: base %2d - ribose %2d\n", i, j);
			}
			if(abs(allBasePhoE[i*seqLen+j] - tmpBasePhoE[i*seqLen+j]) > 0.001){
				printf("EnergyChange: base %2d - pho %2d\n", i, j);
			}
			if(abs(allRibosePhoE[i*seqLen+j] - tmpRibosePhoE[i*seqLen+j]) > 0.001){
				printf("EnergyChange: ribose %2d - pho %2d\n", i, j);
			}
		}
	}
}

void BRFoldingTree::checkTotalEnergy(double shift){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;


	//base-base
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e1 = getBaseBaseEnergy(nodeA, nodeB, sep, et, verbose);
			e2 = baseBaseClash(nodeA, nodeB, sep, et,  shift, verbose);
			//e = getBaseBaseEnergy(nodeA, nodeB, sep, et, verbose) + baseBaseClash(nodeA, nodeB, sep, et, verbose);
			if(abs(e1 - this->allBaseBaseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d base: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBaseBaseE[i*seqLen+j], e1);
			}

			if(abs(e2 - this->allBaseClashE[i*seqLen+j]) > 0.001) {
				printf("DIFF-clash: base: %2d base: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,allBaseClashE[i*seqLen+j], e2);
			}
		}
	}

	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergy(nodeA, nodeB, sep, et, verbose);
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
			e = getBasePhoEnergy(nodeA, nodeB, sep, et, verbose);
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
			e = getRiboseRiboseEnergy(nodeA, nodeB, sep, et, verbose);
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
			e = getRibosePhoEnergy(nodeA, nodeB, sep, et, verbose);
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
			e = getPhoPhoEnergy(nodeA, nodeB, sep, et, verbose);
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


void BRFoldingTree::checkTmpTotalEnergy(double shift){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e, e1, e2;

	//base-base
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e1 = getBaseBaseEnergyTmp(nodeA, nodeB, sep, et, verbose);
			e2 = baseBaseClashTmp(nodeA, nodeB, sep, et, shift, verbose);
			//e = getBaseBaseEnergy(nodeA, nodeB, sep, et, verbose) + baseBaseClash(nodeA, nodeB, sep, et, verbose);
			if(abs(e1 - this->tmpBaseBaseE[i*seqLen+j]) > 0.001){
				printf("DIFF: base: %2d base: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBaseBaseE[i*seqLen+j], e1);
			}

			if(abs(e2 - this->tmpBaseClashE[i*seqLen+j]) > 0.001) {
				printf("DIFF-clash: base: %2d base: %2d energyInTable: %7.3f energy: %7.3f\n", i,j,tmpBaseClashE[i*seqLen+j], e2);
			}
		}
	}

	//base-ribose
	for(i=0;i<seqLen;i++){
		BRNode* nodeA = nodes[i];
		for(j=0;j<seqLen;j++){
			BRNode* nodeB = nodes[j];
			int sep = sepTable[i*seqLen+j];
			e = getBaseRiboseEnergyTmp(nodeA, nodeB, sep, et, verbose);
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
			e = getBasePhoEnergyTmp(nodeA, nodeB, sep, et, verbose);
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
			e = getRiboseRiboseEnergyTmp(nodeA, nodeB, sep, et, verbose);
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
			e = getRibosePhoEnergyTmp(nodeA, nodeB, sep, et, verbose);
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
			e = getPhoPhoEnergyTmp(nodeA, nodeB, sep, et, verbose);
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

BRTreeInfo* BRFoldingTree::getTreeInfo() {
	BRTreeInfo* info = new BRTreeInfo(seqLen, seq, connectToDownstream, nodes, this->flexibleNodes, totalEnergy(1.0, 1.0, 1.0, 0.0, false));
	return info;
}

void BRFoldingTree::printTree(int index) {
	cout << "Node: " << index;
	BRNode* node = nodes[index];
	if(node->leftChild != NULL) {
		BRConnection* ct = node->leftChild->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout<< " Left: " << node->leftChild->seqID << " " << ctType;
	}
	if(node->midChild != NULL) {
		BRConnection* ct = node->midChild->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout << " Mid: " << node->midChild->seqID << " " << ctType;
	}
	if(node->rightChild != NULL) {
		BRNode* rn = node->rightChild;
		BRConnection* ct = rn->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout << " Right: " << rn->seqID << " " << ctType;
	}
	if(node->reverseChild != NULL) {
		BRNode* rn = node->reverseChild;
		BRConnection* ct = rn->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout << " Reverse: " << rn->seqID << " " << ctType;
	}
	cout << endl;

	if(node->leftChild != NULL) {
		printTree(node->leftChild->seqID);
	}
	if(node->midChild != NULL) {
		printTree(node->midChild->seqID);
	}
	if(node->rightChild != NULL) {
		printTree(node->rightChild->seqID);
	}
	if(node->reverseChild != NULL){
		printTree(node->reverseChild->seqID);
	}
}

int BRFoldingTree::printConnections() {

	/*
	for(int i=0;i<this->seqLen;i++){
		for(int j=0;j<this->seqLen;j++) {
			printf("%-2d ",this->sepTable[i*seqLen+j]);
		}
		cout << endl;
	}
	*/

	int x = 0;

	for(int i=0;i<this->seqLen;i++){
		if(nodes[i]->father == NULL) x++;
		//nodes[i]->printPartition();
	}

	for(int i=0;i<fixedConnectionList.size();i++) {
		BRConnection* ct = fixedConnectionList[i];
		cout << ct->fatherNode->seqID << " " << ct->childNode->seqID  << endl;

		//fixedConnectionList[i]->printPartition();
		//fixedConnectionList[i]->printF3Partition();
	}
	for(int i=0;i<flexibleConnectionList.size();i++) {
		BRConnection* ct = flexibleConnectionList[i];
		cout << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
		//flexibleConnectionList[i]->printPartition();
		//flexibleConnectionList[i]->printF3Partition();
	}
	return x;
}



BRFoldingTree::~BRFoldingTree() {


	delete [] seq;
	delete [] wcPairPosID;
	delete [] nwcPairPosID;


	delete [] fixed;
	delete [] loopRevPoints;
	delete [] chainBreakPoints;
	delete [] connectToDownstream;

	delete [] nodeConnectMatrix;

	delete [] sepTable;

	delete [] allBaseClashE;
	delete [] tmpBaseClashE;
	delete [] allBaseBaseE;
	delete [] tmpBaseBaseE;
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

	delete fragLib;
	delete rotLib;

	if(pseudoNode != NULL)
		delete pseudoNode;

	if(seqLen > 0) {
		for(int i=0;i<seqLen;i++)
		{
			delete initRotList[i];
			delete nodes[i];
		}
	}

	delete [] initRotList;
	delete [] nodes;


	delete initTreeInfo;

	for(int i=0;i<this->fixedConnectionList.size();i++){
		delete this->fixedConnectionList[i];
	}



	for(int i=0;i<this->flexibleConnectionList.size();i++){
		delete this->flexibleConnectionList[i];
	}
}

} /* namespace NSPpred */
