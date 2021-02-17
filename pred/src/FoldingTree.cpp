/*
 * FoldingTree.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#include "pred/FoldingTree.h"

namespace NSPpred {

void treeInfo::printPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		if(connectToDownstream[i])
			seqID ++;
		else
			seqID += 5;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), 'A', s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(unsigned int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
	}
	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of.close();
}

double treeInfo::rmsd(treeInfo* other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;
	RnaAtomLib atLib;
	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		if(nodes[i]->fixed) continue;
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		if(nodes[i]->fixed) continue;
		vector<Atom*> aList = other->nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
	}
	return NSPgeometry::rmsd(tList1, tList2);

}

treeInfo::~treeInfo() {
	delete[] seq;
	delete[] connectToDownstream;
	for(unsigned int i=0;i<seqLen;i++){
		delete nodes[i];
	}
	delete[] nodes;
}

FoldingTree::FoldingTree(const string& baseSeq, const string& ssSeq) {

	/*
	 * input format:
	 * single chain
	 * baseSeq: UGGCAAUAGCCGGAUG
	 * ssSeq:   .(((.[[.)))..]].
	 * two chains:
	 * baseSeq: GGGCCAAG&CUUGGCCC
	 * ssSeq:   ((((((((&))))))))
	 */
	if(baseSeq.length() != ssSeq.length()) {
		cout << "invalid input: " << baseSeq << " " << ssSeq << endl;
	}
	int n = 0;
	RNABaseName rn;
	for(unsigned int i=0;i<baseSeq.length();i++) {
		int baseType = rn.sinToInt(baseSeq[i]);
		if(baseType >=0  && baseType <= 3) {
			n++;
		}
	}

	this->seqLen = n;
	this->seq = new int[n];
	this->nodes = new BaseNode*[n];
	this->fixed = new bool[n];
	this->connectToDownstream = new bool[n];
	for(unsigned int i=0;i<n;i++) {
		connectToDownstream[i] = false;
		fixed[i] = false;
	}

	char ss[n];
	int index = 0;
	LocalFrame cs;
	for(unsigned int i=0;i<baseSeq.length();i++) {
		int baseType = rn.sinToInt(baseSeq[i]);
		if(baseType >= 0 && baseType <= 3) {
			this->seq[index] = rn.sinToInt(baseSeq[i]);
			this->nodes[index] = new BaseNode(cs, seq[index], index);
			ss[index] = ssSeq[i];

			if(i<baseSeq.length()-1 && baseType >= 0 && baseType <= 3)
				connectToDownstream[index] = true;
			index++;
		}
	}

	this->initTreeInfo = new treeInfo(seqLen, seq, connectToDownstream, nodes, 0.0);

	/*
	 * parse secondary structure
	 */
	this->wcPairPosID = new int[n];
	this->nwcPairPosID = new int[n];
	for(unsigned int i=0;i<n;i++) {
		wcPairPosID[i] = -1;
		nwcPairPosID[i] = -1;
	}
	for(unsigned int i=0;i<n;i++) {
		char c = ss[i];
		if(c == ')') {
			int preIndex = -1;
			for(unsigned int j=i-1;j>=0;j--) {
				if(ss[j] == '(') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << ssSeq << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == ']') {
			int preIndex = -1;
			for(unsigned int j=i-1;j>=0;j--) {
				if(ss[j] == '[') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << ssSeq << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
		if(c == '}') {
			int preIndex = -1;
			for(unsigned int j=i-1;j>=0;j--) {
				if(ss[j] == '{') {
					preIndex = j;
					break;
				}
			}
			if(preIndex < 0) {
				cout << "invalid ssSeq: " << ssSeq << endl;
				exit(1);
			}
			ss[i] = '.';
			ss[preIndex] = '.';
			wcPairPosID[i] = preIndex;
			wcPairPosID[preIndex] = i;
		}
	}



	/*
	 * build tree
	 */
	cout << "tree" << endl;
	buildFrom(nodes[0]);


//	this->nbConnectionOfNode = new BaseConnection*[seqLen];
//	this->wcConnectionOfNode = new BaseConnection*[seqLen];
	this->connectionToFatherNode = new BaseConnection*[seqLen];
	connectionToFatherNode[0] = NULL;
	/*
	 * update base connection
	 */
	cout << "update connection" << endl;
	for(unsigned int i=0;i<flexNbConnectList.size();i++){
		flexNbConnectList[i]->updateChildInfo();
		this->connectionToFatherNode[flexNbConnectList[i]->childNode->seqID] = flexNbConnectList[i];
//		this->nbConnectionOfNode[nbConnectList[i]->fatherNode->seqID] = nbConnectList[i];
	}
	for(unsigned int i=0;i<flexWcConnectList.size();i++){
		flexWcConnectList[i]->updateChildInfo();
		this->connectionToFatherNode[flexWcConnectList[i]->childNode->seqID] = flexWcConnectList[i];
//		this->wcConnectionOfNode[wcConnectList[i]->fatherNode->seqID] = wcConnectList[i];
	}





	/*
	 * init connect map
	 */
	this->connectMap = new int*[seqLen];
	for(unsigned int i=0;i<seqLen;i++){
		this->connectMap[i] = new int[seqLen];
		for(unsigned int j=0;j<seqLen;j++) {
			this->connectMap[i][j] = 0;
		}
	}

	if(nodes[0]->leftChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->leftChild);
	}
	if(nodes[0]->rightChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->rightChild);
	}

	et = new EnergyTable();

	/*
	 * update tPho
	 */
	cout << "update tPho" << endl;
	for(int i=1;i<seqLen;i++){
		if(!connectToDownstream[i-1])
			continue;
		BaseNode* nodePre = this->nodes[i-1];
		BaseNode* node = this->nodes[i];
		XYZ tPho = this->et->getPhoCoord(nodePre->cs, node->cs, seq[i-1], seq[i]);
		node->updatePho(tPho);
		node->updateTmpPho(tPho);
	}

	updatePsudoConnection();
	pseudoConnectEne = 0.0;
	for(int i=0;i<pseudoConnectList.size();i++) {
		BaseDistanceMatrix dm(pseudoConnectList[i]->fatherNode->cs, pseudoConnectList[i]->childNode->cs);
		double d = et->nearestDistanceNeighbor(dm, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
		this->pseudoConnectEne += d*d;
	}

}

FoldingTree::FoldingTree(const string& inputFile) {
	NSPtools::InputParser input(inputFile);

	input.printOptions();

	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string nwcSec = input.getValue("nwc");

	string flexOption = input.getValue("fixed");
	string fixedGroup = input.getValue("group");

	vector<string> spt;
	vector<string> spt2;
	vector<int> fixedList;
	splitString(flexOption, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		fixedList.push_back(atoi(spt[i].c_str()));
	}

	splitString(fixedGroup, " ", &spt);
	for(unsigned int i=0;i<spt.size();i++){
		splitString(spt[i], ",", &spt2);
		set<int> group;
		for(unsigned int j=0;j<spt2.size();j++) {
			group.insert(atoi(spt2[j].c_str()));
		}
		this->fixedGroups.push_back(group);
	}

	RNAPDB pdb(pdbFile, "xxxx");
	vector<RNABase*> baseList = pdb.getBaseList();
	if(baseList.size() != baseSeq.length()) {
		cout << "pdb size not equal to sequence length" << endl;
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
	this->nodes = new BaseNode*[seqLen];

	RNABaseName rn;
	for(unsigned int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;
		this->nwcPairPosID[i] = -1;
		this->fixed[i] = false;
		if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1]))
			this->connectToDownstream[i] = true;
		LocalFrame cs = baseList[i]->getCoordSystem();
		Atom* p = baseList[i]->getAtom("P");
		XYZ tp;
		if(p != NULL)
			tp = p->getCoord();

		this->nodes[i] = new BaseNode(cs, tp, seq[i], i);
	}

	for(unsigned int i=0;i<fixedList.size();i++) {
		this->fixed[fixedList[i]] = true;
		this->nodes[fixedList[i]]->fixed = true;
	}

	this->initTreeInfo = new treeInfo(seqLen, seq, connectToDownstream, nodes, 0.0);

	/*
	 * parse secondary structure
	 */

	char ss[seqLen];
	for(unsigned int i=0;i<seqLen;i++){
		ss[i] = baseSec[i];
	}

	for(unsigned int i=0;i<seqLen;i++) {
		char c = ss[i];
		if(c == ')') {
			int preIndex = -1;
			for(unsigned int j=i-1;j>=0;j--) {
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
			for(unsigned int j=i-1;j>=0;j--) {
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
			for(unsigned int j=i-1;j>=0;j--) {
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
	}

	/*
	 * parse non-Watson-Crick pair
	 */
	char nss[seqLen];
	for(unsigned int i=0;i<seqLen;i++){
		nss[i] = nwcSec[i];
	}

	for(unsigned int i=0;i<seqLen;i++) {
		char c = nss[i];
		if(c == ')') {
			int preIndex = -1;
			for(unsigned int j=i-1;j>=0;j--) {
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
			for(unsigned int j=i-1;j>=0;j--) {
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
			for(unsigned int j=i-1;j>=0;j--) {
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
	}


	if(input.specifiedOption("wtConnect")){
		double cutoff = atof(input.getValue("wtConnect").c_str());
		cout << "set wild type connect cutoff " << cutoff << endl;
		buildConnectionNearWildType(nodes[0], cutoff);
	}
	else {
		cout << "no wild type connect constraint" << endl;
		buildFrom(nodes[0]);
	}

	this->connectionToFatherNode = new BaseConnection*[seqLen];
	for(unsigned int i=0;i<seqLen;i++) {
		if(this->fixed[i])
			this->connectionToFatherNode[i] = NULL;
	}



	for(unsigned int i=0;i<moveFixedConnectList.size();i++){
		moveFixedConnectList[i]->updateChildInfo();
		this->connectionToFatherNode[moveFixedConnectList[i]->childNode->seqID] = moveFixedConnectList[i];
	}

	for(unsigned int i=0;i<childFixedConnecList.size();i++){
		this->connectionToFatherNode[childFixedConnecList[i]->childNode->seqID] = childFixedConnecList[i];
	}

	for(unsigned int i=0;i<flexNbConnectList.size();i++){
		flexNbConnectList[i]->updateChildInfo();
		this->connectionToFatherNode[flexNbConnectList[i]->childNode->seqID] = flexNbConnectList[i];
	}


	for(unsigned int i=0;i<flexWcConnectList.size();i++){
		flexWcConnectList[i]->updateChildInfo();
		this->connectionToFatherNode[flexWcConnectList[i]->childNode->seqID] = flexWcConnectList[i];
	}


	for(unsigned int i=0;i<moveFixedConnectList.size();i++){
		moveFixedConnectList[i]->updateChildConnectionInfo(connectionToFatherNode);
	}

	for(unsigned int i=0;i<flexNbConnectList.size();i++){
		flexNbConnectList[i]->updateChildConnectionInfo(connectionToFatherNode);
	}
	for(unsigned int i=0;i<flexWcConnectList.size();i++){
		flexWcConnectList[i]->updateChildConnectionInfo(connectionToFatherNode);
	}


	//printTreeInfo();

	if(nodes[0]->leftChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->leftChild);
	}
	if(nodes[0]->rightChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->rightChild);
	}

	/*
	cout << "update cs" << endl;
	for(int i=1;i<seqLen;i++) {
		updateChildCs(nodes[i]);
		updateChildTmpCs(nodes[i]);
	}
	*/


	et = new EnergyTable();


	/*
	for(int i=1;i<seqLen;i++){
		if(!connectToDownstream[i-1])
			continue;
		BaseNode* nodePre = this->nodes[i-1];
		BaseNode* node = this->nodes[i];
		XYZ tPho = this->et->getPhoCoord(nodePre->cs, node->cs, seq[i-1], seq[i]);
		node->updatePho(tPho);
		node->updateTmpPho(tPho);
	}
	*/

	/*
	 * init connect map
	 */
	cout << "connect map" << endl;


	this->connectMap = new int*[seqLen];
	for(unsigned int i=0;i<seqLen;i++){
		this->connectMap[i] = new int[seqLen];
		for(unsigned int j=0;j<seqLen;j++) {
			this->connectMap[i][j] = 0;
		}
	}

	/*
	updatePsudoConnection();
	pseudoConnectEne = 0.0;
	for(int i=0;i<pseudoConnectList.size();i++) {
		BaseDistanceMatrix dm(pseudoConnectList[i]->fatherNode->cs, pseudoConnectList[i]->childNode->cs);
		double d = et->nearestDistanceNeighbor(dm, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
		this->pseudoConnectEne += d*d;
	}
	*/
	cout << "finish init" << endl;
}


void FoldingTree::buildFrom(BaseNode* node) {

	int nextIndex = node->seqID + 1;
	if(wcPairPosID[node->seqID] > 0) {
		BaseNode* wc = nodes[wcPairPosID[node->seqID]];
		if(wc ->father  == NULL) {
			node->rightChild = wc;
			wc->father = node;
			string type = "wc";
			BaseConnection* connect = new BaseConnection(type, node, wc, seqLen);
			bool isFixed = false;
			int idA = node->seqID;
			int idB = wcPairPosID[idA];
			for(unsigned int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[idB]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexWcConnectList.push_back(connect);
		}
	}
	else if(nwcPairPosID[node->seqID] > 0) {
		BaseNode* nwc = nodes[nwcPairPosID[node->seqID]];
		if(nwc ->father  == NULL) {
			node->rightChild = nwc;
			nwc->father = node;
			string type = "nwc";
			BaseConnection* connect = new BaseConnection(type, node, nwc, seqLen);
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nwcPairPosID[idA];
			for(unsigned int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[idB]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexWcConnectList.push_back(connect);
		}
	}

	if(nextIndex < this->seqLen) {
		BaseNode* nb = nodes[nextIndex];
		if(nb->father == NULL ) {
			node->leftChild = nb;
			nb->father = node;
			string type = "loopNb";
			if(!connectToDownstream[node->seqID] || nextIndex > node->seqID+1)
				type = "jump";
			else if(wcPairPosID[node->seqID] > 0 && wcPairPosID[nextIndex] > 0 &&
					wcPairPosID[node->seqID] == wcPairPosID[nextIndex] + 1)
				type = "helixNb";
			else if((wcPairPosID[node->seqID] > 0 && nwcPairPosID[nextIndex] > 0  && wcPairPosID[node->seqID] == nwcPairPosID[nextIndex] + 1) ||
					(nwcPairPosID[node->seqID] > 0 && wcPairPosID[nextIndex] > 0  && nwcPairPosID[node->seqID] == wcPairPosID[nextIndex] + 1) ||
					(nwcPairPosID[node->seqID] > 0 && nwcPairPosID[nextIndex] > 0  && nwcPairPosID[node->seqID] == nwcPairPosID[nextIndex] + 1))
				type = "nwcNb";
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nextIndex;
			for(unsigned int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}

			BaseConnection* connect = new BaseConnection(type, node, nb, seqLen);
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[nextIndex]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexNbConnectList.push_back(connect);
		}
	}

	if(nextIndex >= this->seqLen)
		return;
	buildFrom(nodes[nextIndex]);

}

void FoldingTree::updatePsudoConnection() {
	this->pseudoConnectList.clear();
	bool* tmp = new bool[this->seqLen];
	for(int i=0;i<seqLen;i++){
		tmp[i] = false;
	}
	for(int i=0;i<this->moveFixedConnectList.size();i++) {
		BaseConnection* bc = this->moveFixedConnectList[i];
		if(bc->isNeighborConnect())
			tmp[bc->fatherNode->seqID] = true;
		if(bc->connectTye == "wc")
			tmp[bc->childNode->seqID] = true;
	}
	for(int i=0;i<this->childFixedConnecList.size();i++) {
		BaseConnection* bc = this->childFixedConnecList[i];
		if(bc->isNeighborConnect())
			tmp[bc->fatherNode->seqID] = true;
		if(bc->connectTye == "wc")
					tmp[bc->childNode->seqID] = true;
	}
	for(int i=0;i<this->flexNbConnectList.size();i++) {
		tmp[flexNbConnectList[i]->fatherNode->seqID] = true;
	}

	for(int i=0;i<this->flexWcConnectList.size();i++) {
		tmp[flexWcConnectList[i]->childNode->seqID] = true;
	}

	for(int i=0;i<this->seqLen-1;i++) {
		if(connectToDownstream[i] && !tmp[i]) {
			bool isFixed = false;
			int idA = i;
			int idB = i+1;
			for(unsigned int k=0;k<fixedGroups.size();k++) {

				if(fixedGroups[k].find(idA) != fixedGroups[k].end() && fixedGroups[k].find(idB) != fixedGroups[k].end())
					isFixed = true;
			}
			if(isFixed) continue;
			BaseConnection* bc = new BaseConnection("loopNb", this->nodes[i], nodes[i+1], seqLen);
			bc->setNativeMove();
			this->pseudoConnectList.push_back(bc);
		}
	}
}

void FoldingTree::buildConnectionNearWildType(BaseNode* node, double cutoff) {

	int nextIndex = node->seqID + 1;
	if(wcPairPosID[node->seqID] > 0) {
		BaseNode* wc = nodes[wcPairPosID[node->seqID]];
		if(wc ->father  == NULL) {
			node->rightChild = wc;
			wc->father = node;
			string type = "wc";
			CsMove mv = node->cs.getMove(wc->cs);
			BaseConnection* connect = new BaseConnection(type, &mv, cutoff, node, wc, seqLen);
			//cout << "new connect: " << node->seqID << "-" << wc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = wcPairPosID[idA];
			for(unsigned int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[idB]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexWcConnectList.push_back(connect);
		}
	}
	else if(nwcPairPosID[node->seqID] > 0) {
		BaseNode* nwc = nodes[nwcPairPosID[node->seqID]];
		if(nwc ->father  == NULL) {
			node->rightChild = nwc;
			nwc->father = node;
			string type = "nwc";
			CsMove mv = node->cs.getMove(nwc->cs);
			BaseConnection* connect = new BaseConnection(type, &mv, cutoff, node, nwc, seqLen);
			//cout << "new connect: " << node->seqID << "-" << nwc->seqID << endl;
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nwcPairPosID[idA];
			for(unsigned int i=0;i<fixedGroups.size();i++) {

				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[idB]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexWcConnectList.push_back(connect);
		}
	}

	if(nextIndex < this->seqLen) {
		BaseNode* nb = nodes[nextIndex];
		if(nb->father == NULL ) {
			node->leftChild = nb;
			nb->father = node;
			string type = "loopNb";
			if(!connectToDownstream[node->seqID] || nextIndex > node->seqID+1)
				type = "jump";
			else if(wcPairPosID[node->seqID] > 0 && wcPairPosID[nextIndex] > 0 &&
					wcPairPosID[node->seqID] == wcPairPosID[nextIndex] + 1)
				type = "helixNb";
			bool isFixed = false;
			int idA = node->seqID;
			int idB = nextIndex;
			for(unsigned int i=0;i<fixedGroups.size();i++) {
				if(fixedGroups[i].find(idA) != fixedGroups[i].end() && fixedGroups[i].find(idB) != fixedGroups[i].end())
					isFixed = true;
			}

			CsMove mv = node->cs.getMove(nb->cs);
			BaseConnection* connect = new BaseConnection(type, &mv, cutoff, node, nb, seqLen);
			//cout << "new connect: " << node->seqID << "-" << nb->seqID << endl;
			if(isFixed) {
				connect->setMoveFixed();
				this->moveFixedConnectList.push_back(connect);
			}
			else if(fixed[nextIndex]) {
				connect->setChildFixed();
				this->childFixedConnecList.push_back(connect);
			}
			else
				this->flexNbConnectList.push_back(connect);
		}
	}

	if(nextIndex >= this->seqLen)
		return;
	buildConnectionNearWildType(nodes[nextIndex], cutoff);

}

void FoldingTree::randomInit() {
	int nbConnectNum = this->flexNbConnectList.size();
	int wcConnectNum = this->flexWcConnectList.size();
	for(int i=0;i<nbConnectNum;i++) {
		this->flexNbConnectList[i]->initRandom();
	}
	for(int i=0;i<wcConnectNum;i++) {
		this->flexWcConnectList[i]->initRandom();
	}
	if(nodes[0]->leftChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->leftChild);
	}
	if(nodes[0]->rightChild != NULL) {
		updateChildCsAndTmpCs(nodes[0]->rightChild);
	}


	for(int i=1;i<seqLen;i++){
		if(!connectToDownstream[i-1])
			continue;
		BaseNode* nodePre = this->nodes[i-1];
		BaseNode* node = this->nodes[i];
		XYZ tPho = this->et->getPhoCoord(nodePre->cs, node->cs, seq[i-1], seq[i]);
		node->updatePho(tPho);
		node->updateTmpPho(tPho);
	}
}

void FoldingTree::initFromTreeInfo(treeInfo* ti){
	for(int i=0;i<seqLen;i++) {
		nodes[i]->copyValueFrom(ti->nodes[i]);
	}

	int fixConnectNum = this->moveFixedConnectList.size();
	int childFixedConnectNum = this->childFixedConnecList.size();
	int nbConnectNum = this->flexNbConnectList.size();
	int wcConnectNum = this->flexWcConnectList.size();

	for(int i=0;i<fixConnectNum;i++) {
		this->moveFixedConnectList[i]->setNativeMove();
	}
	for(int i=0;i<childFixedConnectNum;i++) {
		this->childFixedConnecList[i]->setNativeMove();
	}
	for(int i=0;i<nbConnectNum;i++) {
		this->flexNbConnectList[i]->setNativeMove();
	}
	for(int i=0;i<wcConnectNum;i++) {
		this->flexWcConnectList[i]->setNativeMove();
	}

}

void FoldingTree::updateConnectMap(){

	for(unsigned int i=0;i<seqLen;i++) {
		for(unsigned int j=0;j<seqLen;j++) {
			connectMap[i][j] = 0;
		}
	}

	for(unsigned int i=0;i<seqLen-1;i++) {
		connectMap[i][i+1] = 1;
		connectMap[i+1][i] = 1;
	}

	for(unsigned int i=0;i<seqLen;i++) {
		XYZ ci = nodes[i]->baseCenter;
		for(unsigned int j=i+2;j<seqLen;j++){
			XYZ cj = nodes[j]->baseCenter;

			if(squareDistance(ci,cj) < 100.0) {
				connectMap[i][j] = 1;
				connectMap[j][i] = 1;
			}
		}
	}

}


double FoldingTree::partialEnergy(BaseConnection* selectConnect) {

	int idA, idB;
	double totE, e, centerDist;
	totE = 0;

	int childNodeNum = selectConnect->childNodeIndexList.size();
	int nonChildNodeNum = selectConnect->nonChildNodeIndexList.size();

	int updatedPhoNum = selectConnect->updatedPhoIndexList.size();
	int unchangedPhoNum = selectConnect->unchangedPhoIndexList.size();

	/*
	 * base-pho energy change
	 */


	for(int i=0;i<updatedPhoNum;i++) {
		idA = selectConnect->updatedPhoIndexList[i];
		XYZ localA = this->nodes[idA]->cs.global2localcrd(nodes[idA]->tPho);
		totE += et->getBasePhoEnergy(localA, 0);
		if(idA > 0 && this->connectToDownstream[idA-1]) {
			XYZ localB = this->nodes[idA-1]->cs.global2localcrd(nodes[idA]->tPho);
			totE += et->getBasePhoEnergy(localB, 1);
		}
	}



	/*
	 * pho-pho energy change
	 */


	for(int i=0;i<updatedPhoNum;i++){
		idA = selectConnect->updatedPhoIndexList[i];
		for(int j=i+1;j<updatedPhoNum;j++){
			idB = selectConnect->updatedPhoIndexList[j];
			double d = nodes[idA]->tPho.distance(nodes[idB]->tPho);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;
			if(sep == 1)
				totE += et->getPhoPhoEnergy(d, 1);
			else if(d < 5)
				totE += et->getPhoPhoEnergy(d, 2);
		}

		for(int j=0;j<unchangedPhoNum;j++){
			idB = selectConnect->unchangedPhoIndexList[j];
			double d = nodes[idA]->tPho.distance(nodes[idB]->tPho);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;
			if(sep == 1)
				totE += et->getPhoPhoEnergy(d, 1);
			else if(d < 5)
				totE += et->getPhoPhoEnergy(d, 2);
		}
	}


	/*
	 * base-base energy change
	 */


	for(unsigned int i=0;i<childNodeNum;i++){
		idA = selectConnect->childNodeIndexList[i];
		for(unsigned int j=0;j<nonChildNodeNum;j++){
			idB = selectConnect->nonChildNodeIndexList[j];
			if(idB - idA == 1 && connectToDownstream[idA]) {
				BaseDistanceMatrix dm(this->nodes[idA]->cs, this->nodes[idB]->cs);
				e = et->getBaseBaseEnergyNeighbor(dm, seq[idA], seq[idB]);
				//printf("Pnb1 IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
				totE += e;
			}
			else if(idA - idB == 1 && connectToDownstream[idB]) {
				BaseDistanceMatrix dm(this->nodes[idB]->cs, this->nodes[idA]->cs);
				e = et->getBaseBaseEnergyNeighbor(dm, seq[idB], seq[idA]);
				totE += e;
				//printf("Pnb2 IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
			}
			else if(idA < idB){
				BaseDistanceMatrix dm(this->nodes[idA]->cs, this->nodes[idB]->cs);
				centerDist = sqrt(squareDistance(nodes[idA]->baseCenter, nodes[idB]->baseCenter));
				if(centerDist > 10.0) continue;
				e = et->getBaseBaseEnergy(dm, seq[idA], seq[idB], centerDist);
				totE += e;
				//printf("P IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
			}
			else if(idA > idB){
				BaseDistanceMatrix dm(this->nodes[idB]->cs, this->nodes[idA]->cs);
				centerDist = sqrt(squareDistance(nodes[idA]->baseCenter, nodes[idB]->baseCenter));
				if(centerDist > 10.0) continue;
				e = et->getBaseBaseEnergy(dm, seq[idB], seq[idA], centerDist);
				totE += e;
				//printf("P IDA: %2d IDB: %2d ENE: %7.3f\n",idB, idA,e);
			}
		}
	}



	return totE;
}

void FoldingTree::updateChildTmpCs(BaseNode* node) {

	if(node->fixed) {

		if(node->leftChild != NULL){
			updateChildTmpCs(node->leftChild);
		}
		if(node->rightChild != NULL) {
			updateChildTmpCs(node->rightChild);
		}
		return;
	}

	LocalFrame cs = node->father->tmpCs.add(connectionToFatherNode[node->seqID]->cm);
	node->updateTmpCs(cs);
	cout << "updateA: " << node->seqID << endl;
	if(node->leftChild != NULL){
		updateChildTmpCs(node->leftChild);
	}
	if(node->rightChild != NULL) {
		updateChildTmpCs(node->rightChild);
	}
}

void FoldingTree::updateChildTmpPho(BaseConnection* selectConnect) {
	int n = selectConnect->updatedPhoIndexList.size();

	for(int i=0;i<n;i++){
		int id = selectConnect->updatedPhoIndexList[i];
		if(!connectToDownstream[id-1])
			continue;

		bool preCsUpdated = selectConnect->childOrNotChild[id-1];
		bool curCsUpdated = selectConnect->childOrNotChild[id];

		BaseNode* nodePre = this->nodes[id-1];
		BaseNode* node = this->nodes[id];
		if(preCsUpdated && curCsUpdated) {
			XYZ tPho = this->et->getPhoCoord(nodePre->tmpCs, node->tmpCs, seq[id-1], seq[id]);
			node->updateTmpPho(tPho);
			//cout << "both tmp: " << id << " " << tPho.toString() << endl;
		}
		else if(curCsUpdated){
			XYZ tPho = this->et->getPhoCoord(nodePre->cs, node->tmpCs, seq[id-1], seq[id]);
			node->updateTmpPho(tPho);
			//cout << "cur tmp: " << id << " " << tPho.toString() << endl;
		}
		else if(preCsUpdated){
			XYZ tPho = this->et->getPhoCoord(nodePre->tmpCs, node->cs, seq[id-1], seq[id]);
			node->updateTmpPho(tPho);
			//cout << "pre tmp: " << id << " " << tPho.toString() << endl;
		}

	}
}

void FoldingTree::updateChildCs(BaseNode* node) {

	if(node->fixed) {

		if(node->leftChild != NULL){
			updateChildCs(node->leftChild);
		}
		if(node->rightChild != NULL) {
			updateChildCs(node->rightChild);
		}
		return;
	}
	LocalFrame cs = node->father->cs.add(connectionToFatherNode[node->seqID]->cm);
	node->updateCs(cs);
	//cout << "update cs: " << node->seqID << endl;

	if(node->leftChild != NULL){
		updateChildCs(node->leftChild);
	}
	if(node->rightChild != NULL) {
		updateChildCs(node->rightChild);
	}
}

void FoldingTree::updateChildCsAndTmpCs(BaseNode* node) {
	if(node->fixed) {
		if(node->leftChild != NULL){
			updateChildCsAndTmpCs(node->leftChild);
		}
		if(node->rightChild != NULL) {
			updateChildCsAndTmpCs(node->rightChild);
		}
		return;
	}

	LocalFrame cs = node->father->cs.add(connectionToFatherNode[node->seqID]->cm);

	node->updateCs(cs);
	node->updateTmpCs(cs);


	if(node->leftChild != NULL){
		updateChildCsAndTmpCs(node->leftChild);
	}
	if(node->rightChild != NULL) {
		updateChildCsAndTmpCs(node->rightChild);
	}
}

void FoldingTree::updateChildPho(BaseConnection* selectConnect) {
	int n = selectConnect->updatedPhoIndexList.size();
	for(int i=0;i<n;i++){
		int id = selectConnect->updatedPhoIndexList[i];

		//if(!connectToDownstream[id-1])
		//	continue;
		//BaseNode* nodePre = this->nodes[id-1];
		BaseNode* node = this->nodes[id];
		node->updatePho(node->tmpPho);
		//XYZ tPho = this->et->getPhoCoord(nodePre->cs, node->cs, seq[id-1], seq[id]);
		//node->updatePho(tPho);
		//cout << "child " << tPho.toString() << endl;
	}
}


void FoldingTree::updateChildTmpCs(BaseConnection* selectConnect, CsMove* mutMove) {


	int childNum = selectConnect->childConnectionList.size();
	LocalFrame cs0 = selectConnect->fatherNode->cs.add(*mutMove);
	selectConnect->childNode->updateTmpCs(cs0);
	//cout << "updateB: " << selectConnect->childNode->seqID << endl;
	BaseConnection* ct;
	for(int i=0;i<childNum;i++) {
		ct = selectConnect->childConnectionList[i];
		LocalFrame cs = ct->fatherNode->tmpCs.add(ct->cm);
		ct->childNode->updateTmpCs(cs);
		//cout << "updateB: " << ct->childNode->seqID << endl;
	}
}

void FoldingTree::updateChildCs(BaseConnection* selectConnect, CsMove* mutMove) {

	int childNum = selectConnect->childConnectionList.size();
	//LocalFrame cs0 = selectConnect->fatherNode->cs.add(*mutMove);
	selectConnect->childNode->updateCs(selectConnect->childNode->tmpCs);
	BaseConnection* ct;
	for(int i=0;i<childNum;i++) {
		ct = selectConnect->childConnectionList[i];
		//LocalFrame cs = ct->fatherNode->cs.add(ct->cm);
		ct->childNode->updateCs(ct->childNode->tmpCs);
	}
}

double FoldingTree::partialMutEnergy(BaseConnection* selectConnect, CsMove* mutMove) {

	/*
	BaseNode* node = selectConnect->childNode;
	LocalFrame cs = node->father->cs.add(*mutMove);
	cout << "updateA: " << node->seqID << endl;
	node->updateTmpCs(cs);
	if(node->leftChild != NULL)
		updateChildTmpCs(node->leftChild);
	if(node->rightChild != NULL)
		updateChildTmpCs(node->rightChild);
	*/



	updateChildTmpCs(selectConnect, mutMove);
	updateChildTmpPho(selectConnect);


	/*
	 * tmp coordinates of all child nodes of "node" have been updated
	 * node(idA): tmpCs
	 * node(idB): cs
	 */
	int idA, idB;
	double totE, e, centerDist;
	totE = 0;

	int childNodeNum = selectConnect->childNodeIndexList.size();
	int nonChildNodeNum = selectConnect->nonChildNodeIndexList.size();

	int updatedPhoNum = selectConnect->updatedPhoIndexList.size();
	int unchangedPhoNum = selectConnect->unchangedPhoIndexList.size();


	/*
	 * base-pho energy change
	 */


	XYZ localA, localB;

	for(int i=0;i<updatedPhoNum;i++) {
		idA = selectConnect->updatedPhoIndexList[i];
		if(selectConnect->childOrNotChild[idA])
			localA = this->nodes[idA]->tmpCs.global2localcrd(nodes[idA]->tmpPho);
		else
			localA = this->nodes[idA]->cs.global2localcrd(nodes[idA]->tmpPho);

		totE += et->getBasePhoEnergy(localA, 0);
		if(idA > 0 && this->connectToDownstream[idA-1]) {
			if(selectConnect->childOrNotChild[idA-1])
				localB = this->nodes[idA-1]->tmpCs.global2localcrd(nodes[idA]->tmpPho);
			else
				localB = this->nodes[idA-1]->cs.global2localcrd(nodes[idA]->tmpPho);
			totE += et->getBasePhoEnergy(localB, 1);
		}
	}



	/*
	 * pho-pho energy change
	 */


	for(int i=0;i<updatedPhoNum;i++){
		idA = selectConnect->updatedPhoIndexList[i];
		for(int j=i+1;j<updatedPhoNum;j++){
			idB = selectConnect->updatedPhoIndexList[j];
			double d = nodes[idA]->tmpPho.distance(nodes[idB]->tmpPho);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;
			if(sep == 1)
				totE += et->getPhoPhoEnergy(d, 1);
			else if(d < 5)
				totE += et->getPhoPhoEnergy(d, 2);
		}

		for(int j=0;j<unchangedPhoNum;j++){
			idB = selectConnect->unchangedPhoIndexList[j];
			double d = nodes[idA]->tmpPho.distance(nodes[idB]->tPho);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;
			if(sep == 1)
				totE += et->getPhoPhoEnergy(d, 1);
			else if(d < 5)
				totE += et->getPhoPhoEnergy(d, 2);
		}
	}


	/*
	 * base-base energy change
	 */


	for(int i=0;i<childNodeNum;i++){
		idA = selectConnect->childNodeIndexList[i];
		for(int j=0;j<nonChildNodeNum;j++){
			idB = selectConnect->nonChildNodeIndexList[j];
			if(idB - idA == 1 && connectToDownstream[idA]) {
				BaseDistanceMatrix dm(this->nodes[idA]->tmpCs, this->nodes[idB]->cs);
				e = et->getBaseBaseEnergyNeighbor(dm, seq[idA], seq[idB]);
				//printf("Tnb1 IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
				totE += e;
			}
			else if(idA - idB == 1 && connectToDownstream[idB]) {
				BaseDistanceMatrix dm(this->nodes[idB]->cs, this->nodes[idA]->tmpCs);
				e = et->getBaseBaseEnergyNeighbor(dm, seq[idB], seq[idA]);
				totE += e;
				//printf("Tnb2 IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
			}
			else if(idA < idB) {
				BaseDistanceMatrix dm(this->nodes[idA]->tmpCs, this->nodes[idB]->cs);
				centerDist = sqrt(squareDistance(nodes[idA]->tmpBaseCenter, nodes[idB]->baseCenter));
				if(centerDist > 10.0) continue;
				e = et->getBaseBaseEnergy(dm, seq[idA], seq[idB], centerDist);
				totE += e;
				//printf("T IDA: %2d IDB: %2d ENE: %7.3f\n",idA, idB,e);
			}
			else if(idA > idB) {
				BaseDistanceMatrix dm(this->nodes[idB]->cs, this->nodes[idA]->tmpCs);
				centerDist = sqrt(squareDistance(nodes[idA]->tmpBaseCenter, nodes[idB]->baseCenter));
				if(centerDist > 10.0) continue;
				e = et->getBaseBaseEnergy(dm, seq[idB], seq[idA], centerDist);
				totE += e;
				//printf("T IDA: %2d IDB: %2d ENE: %7.3f\n",idB, idA,e);
			}
		}
	}


	return totE;
}

double FoldingTree::dPartialMutPseudoConnectEnergy(BaseConnection* selectConnect, CsMove* mutMove){
	bool connectChange = false;
	for(int i=0;i<pseudoConnectList.size();i++) {
		int idA = pseudoConnectList[i]->fatherNode->seqID;
		int idB = idA + 1;
		if(selectConnect->childOrNotChild[idA] || selectConnect->childOrNotChild[idB]) {
			connectChange = true;
		}
	}
	if(!connectChange) return 0.0;

	BaseNode* node = selectConnect->childNode;
	LocalFrame cs = node->father->cs.add(*mutMove);
	node->updateTmpCs(cs);
	if(node->leftChild != NULL)
		updateChildTmpCs(node->leftChild);
	if(node->rightChild != NULL)
		updateChildTmpCs(node->rightChild);
	double e = 0.0;

	double d1,d2;

	for(int i=0;i<pseudoConnectList.size();i++) {
		int idA = pseudoConnectList[i]->fatherNode->seqID;
		int idB = idA + 1;
		if(selectConnect->childOrNotChild[idA] && !selectConnect->childOrNotChild[idB]) {
			BaseDistanceMatrix dm(nodes[idA]->tmpCs, nodes[idB]->cs);
			BaseDistanceMatrix dm2(nodes[idA]->cs, nodes[idB]->cs);
			d1 = et->nearestDistanceNeighbor(dm, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
			d2 = et->nearestDistanceNeighbor(dm2, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
			e += d1*d1 - d2*d2;
		}
		else if(selectConnect->childOrNotChild[idB] && !selectConnect->childOrNotChild[idA]) {
			BaseDistanceMatrix dm(nodes[idA]->cs, nodes[idB]->tmpCs);
			BaseDistanceMatrix dm2(nodes[idA]->cs, nodes[idB]->cs);

			d1 = et->nearestDistanceNeighbor(dm, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
			d2 = et->nearestDistanceNeighbor(dm2, pseudoConnectList[i]->fatherNode->baseType, pseudoConnectList[i]->childNode->baseType);
			e += d1*d1 - d2*d2;
		}
	}
	return e;
}

void FoldingTree::acceptMove(BaseConnection* selectConnect, CsMove* mutMove){

	/*
	selectConnect->cm = *mutMove;
	BaseNode* node = selectConnect->childNode;
	LocalFrame cs = node->father->cs.add(*mutMove);
	node->updateCs(cs);

	if(node->leftChild != NULL)
		updateChildCs(node->leftChild);
	if(node->rightChild != NULL)
		updateChildCs(node->rightChild);

	updateChildPho(selectConnect);
	*/

	selectConnect->cm = *mutMove;
	updateChildCs(selectConnect, mutMove);
	updateChildPho(selectConnect);

}

void FoldingTree::accetpMoveWithoutPho(BaseConnection* selectConnect, CsMove* mutMove) {
	/*
	selectConnect->cm = *mutMove;
	BaseNode* node = selectConnect->childNode;
	LocalFrame cs = node->father->cs.add(*mutMove);
	node->updateCs(cs);

	if(node->leftChild != NULL)
		updateChildCs(node->leftChild);
	if(node->rightChild != NULL)
		updateChildCs(node->rightChild);
	*/
	selectConnect->cm = *mutMove;
	updateChildCs(selectConnect, mutMove);
}

void FoldingTree::printTree(unsigned int index) {
	cout << "Node: " << index;
	BaseNode* node = nodes[index];

	if(node->leftChild != NULL)
	{
		BaseConnection* bc = connectionToFatherNode[node->leftChild->seqID];
		string bcType = "flex";
		if(bc->moveFixed)
			bcType = "fixed";
		if(bc->childFixed)
			bcType = "cFixed";
		cout << " Left: " << node->leftChild->seqID << " " << bc->connectTye <<  " " << bcType;
	}
	if(node->rightChild != NULL) {
		BaseConnection* bc = connectionToFatherNode[node->rightChild->seqID];
		string bcType = "flex";
		if(bc->moveFixed)
			bcType = "fixed";
		if(bc->childFixed)
			bcType = "cFixed";
		cout << " Right: " << node->rightChild->seqID << " " << bc->connectTye <<  " " << bcType;
	}
	cout << endl;

	if(node->leftChild != NULL) {
		printTree(node->leftChild->seqID);
	}
	if(node->rightChild != NULL) {
		printTree(node->rightChild->seqID);
	}
}


void FoldingTree::printConnection(unsigned int index) {
	//cout << "node: " << index << endl;
	BaseNode* node = nodes[index];
	if(node->leftChild != NULL) {
		BaseConnection* ct = connectionToFatherNode[node->leftChild->seqID];
		LocalFrame csA = node->cs;

		LocalFrame csTmpB = csA.add(ct->cm);
		LocalFrame csB = node->leftChild->cs;

		CsMove m1 = csA.getMove(csTmpB);
		CsMove m2 = ct->cm;
		BaseDistanceMatrix dm(csA, csB);
		BaseDistanceMatrix dm2(ct->cm);
		double p = ct->mutator->proportion;
		cout << "connect: " << node->seqID << " " << node->leftChild->seqID << " " << ct->connectTye << " p: " << p << endl;
		//ct->printChildConnection();
		//cout << m1.toString() << endl;
		//cout << m2.toString() << endl;
	}
	if(node->rightChild != NULL) {
		BaseConnection* ct = connectionToFatherNode[node->rightChild->seqID];
		LocalFrame csA = node->cs;
		LocalFrame csTmpB = csA.add(ct->cm);
		LocalFrame csB = node->rightChild->cs;
		CsMove m1 = csA.getMove(csTmpB);
		CsMove m2 = ct->cm;
		BaseDistanceMatrix dm(csA, csB);
		BaseDistanceMatrix dm2(ct->cm);
		double p = ct->mutator->proportion;
		cout << "connect: " << node->seqID << " " << node->rightChild->seqID << " " << ct->connectTye << " p: " << p << endl;
		//ct->printChildConnection();
		//cout << m1.toString() << endl;
		//cout << m2.toString() << endl;
	}
	if(node->leftChild != NULL)
		printConnection(node->leftChild->seqID);
	if(node->rightChild != NULL)
		printConnection(node->rightChild->seqID);

}

treeInfo* FoldingTree::getTreeInfo(){
	double e = totalEnergy();
	return new treeInfo(seqLen, seq, connectToDownstream, nodes, e);
}

void FoldingTree::printPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	for(unsigned int i=0;i<this->seqLen;i++) {
		sprintf(ss, "%d", i+1);
		RNABase* base = new RNABase(string(ss), 'A', s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(unsigned int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
	}
	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of.close();
}

void FoldingTree::printTmpPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	for(unsigned int i=0;i<this->seqLen;i++) {
		sprintf(ss, "%d", i+1);
		RNABase* base = new RNABase(string(ss), 'A', s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toTmpAtomList(atLib);
		for(unsigned int j=0;j<aList.size();j++)
			base->addAtom(aList[j]);
		rc.addBase(base);
	}
	ofstream of;
	of.open(outputFile, ios::out);
	rc.printPDBFormat(of, 1);
	of.close();
}

double FoldingTree::totalEnergy() {

	double totE = 0.0;
	updateConnectMap();
	double e, centerDist;
	for(int idA=0;idA<seqLen;idA++){
		XYZ pA = nodes[idA]->tPho;


		XYZ localA = this->nodes[idA]->cs.global2localcrd(nodes[idA]->tPho);
		totE += et->getBasePhoEnergy(localA, 0);
		if(idA > 0 && this->connectToDownstream[idA-1]) {
			XYZ localB = this->nodes[idA-1]->cs.global2localcrd(nodes[idA]->tPho);
			totE += et->getBasePhoEnergy(localB, 1);
		}



		for(int idB=idA+1;idB<seqLen;idB++){
			XYZ pB = nodes[idB]->tPho;
			double d = pA.distance(pB);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;


			if(sep == 1)
				totE += et->getPhoPhoEnergy(d, 1);
			else if(d < 5)
				totE += et->getPhoPhoEnergy(d, 2);



			if(connectMap[idA][idB] > 0) {

				centerDist = nodes[idA]->baseCenter.distance(nodes[idB]->baseCenter);
				BaseDistanceMatrix dm(this->nodes[idA]->cs, this->nodes[idB]->cs);
				if(idB - idA == 1 && connectToDownstream[idA]) {
					e = et->getBaseBaseEnergyNeighbor(dm, seq[idA], seq[idB]);
					totE += e;
					//printf("tot: %d %d %7.3f\n",idA, idB, e);
				}
				else {
					centerDist = nodes[idA]->baseCenter.distance(nodes[idB]->baseCenter);
					e = et->getBaseBaseEnergy(dm, seq[idA], seq[idB], centerDist);
					totE += e;
					//printf("tot: %d %d %7.3f\n",idA, idB, e);
				}
			}


		}
	}
	return totE;
}

void FoldingTree::printEnergyDetail(){
	double totE = 0.0;
	updateConnectMap();
	double e, centerDist;
	for(int idA=0;idA<seqLen;idA++){
		for(int idB=idA+1;idB<seqLen;idB++){
			if(connectMap[idA][idB] > 0) {
				centerDist = nodes[idA]->baseCenter.distance(nodes[idB]->baseCenter);
				BaseDistanceMatrix dm(this->nodes[idA]->cs, this->nodes[idB]->cs);
				if(idB - idA == 1 && connectToDownstream[idA]) {
					e = et->getBaseBaseEnergyNeighbor(dm, seq[idA], seq[idB]);
					printf("%d %d %7.3f\n",idA,idB,e);
					totE += e;
				}
				else {
					centerDist = nodes[idA]->baseCenter.distance(nodes[idB]->baseCenter);
					e = et->getBaseBaseEnergy(dm, seq[idA], seq[idB], centerDist);
					printf("%d %d %7.3f\n",idA,idB,e);
					totE += e;
				}
			}
		}
	}

	for(int idA=0;idA<seqLen;idA++){
		XYZ pA = nodes[idA]->tPho;
		for(int idB=idA+1;idB<seqLen;idB++){
			XYZ pB = nodes[idB]->tPho;
			double d = pA.distance(pB);
			int sep;
			if(idA == idB + 1 && this->connectToDownstream[idB])
				sep = 1;
			else if(idB == idA + 1 && this->connectToDownstream[idA])
				sep = 1;
			else
				sep = 2;
			if(sep == 1) {
				e = et->getPhoPhoEnergy(d, 1);
				totE += e;
				printf("%d %d %7.3f\n",idA,idB,e);
			}
			else if(d < 5) {
				e = et->getPhoPhoEnergy(d, 2);
				totE += e;
				printf("%d %d %7.3f\n",idA,idB,e);
			}
		}
	}


	for(int idA=0;idA<seqLen;idA++){
		XYZ pA = nodes[idA]->tPho;
		//printf("%d %7.3f %7.3f %7.3f\n",idA, pA.x_, pA.y_, pA.z_);
		XYZ localA = this->nodes[idA]->cs.global2localcrd(nodes[idA]->tPho);
		e = et->getBasePhoEnergy(localA, 0);
		printf("p0: %d %7.3f\n",idA,e);
		totE += e;
		if(idA > 0 && this->connectToDownstream[idA-1]) {
			XYZ localB = this->nodes[idA-1]->cs.global2localcrd(nodes[idA]->tPho);
			e = et->getBasePhoEnergy(localB, 1);
			printf("p1: %d %7.3f\n", idA, e);
			totE += e;
		}
	}
	printf("totE: %7.3f\n", totE);
}

void FoldingTree::printPho() {
	for(int id=0;id<seqLen;id++) {
		XYZ p = nodes[id]->tPho;
		XYZ tmp = nodes[id]->tmpPho;
		printf("%d P: %7.3f %7.3f %7.3f TmpP: %7.3f %7.3f %7.3f\n",id, p.x_, p.y_, p.z_, tmp.x_, tmp.y_, tmp.z_);
	}
}

FoldingTree::~FoldingTree() {
	// TODO Auto-generated destructor stub

	delete[] seq;
	delete[] wcPairPosID;
	delete[] nwcPairPosID;


	delete[] connectToDownstream;
	delete[] connectionToFatherNode;


	for(unsigned int i=0;i<moveFixedConnectList.size();i++){
		delete moveFixedConnectList[i];
	}

	for(unsigned int i=0;i<childFixedConnecList.size();i++) {
		delete childFixedConnecList[i];
	}


	for(unsigned int i=0;i<flexNbConnectList.size();i++) {
		delete flexNbConnectList[i];
	}


	for(unsigned int i=0;i<flexWcConnectList.size();i++){
		delete flexWcConnectList[i];
	}


	for(unsigned int i=0;i<seqLen;i++){
		delete nodes[i];
	}


	for(unsigned int i=0;i<seqLen;i++) {
		delete[] connectMap[i];
	}

	delete[] connectMap;

	if(initTreeInfo != NULL)
		delete initTreeInfo;

	if(et != NULL)
		delete et;

}

} /* namespace NSPpred */
