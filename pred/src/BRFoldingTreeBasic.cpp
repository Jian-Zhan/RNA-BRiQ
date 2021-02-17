/*
 * BRFoldingTreeBasic.cpp
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#include <pred/BRFoldingTreeBasic.h>
namespace NSPpred {

BRTreeInfoBasic::BRTreeInfoBasic(int seqLen, int* seq, bool* con, BRNodeBasic** nodes, double ene) {
	this->seqLen = seqLen;
	this->seq = new int[seqLen];
	this->connectToDownstream = new bool[seqLen];
	this->nodes = new BRNodeBasic*[seqLen];
	this->ene = ene;
	for(int i=0;i<seqLen;i++){
		this->seq[i] = seq[i];
		this->connectToDownstream[i] = con[i];
	}
	for(int i=0;i<seqLen;i++) {
		BRNodeBasic* br = new BRNodeBasic(nodes[i]->baseType, nodes[i]->seqID);
		br->copyValueFrom(*nodes[i]);
		this->nodes[i] = br;
	}
	this->ene = ene;
}

double BRTreeInfoBasic::rmsd(BRTreeInfoBasic* other){
	vector<XYZ> tList1;
	vector<XYZ> tList2;
	RnaAtomLib atLib;
	int seqID = 0;
	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList1.push_back(aList[j]->coord);
	}

	for(unsigned int i=0;i<this->seqLen;i++) {
		vector<Atom*> aList = other->nodes[i]->toAtomList(atLib);
		for(int j=0;j<aList.size();j++)
			tList2.push_back(aList[j]->coord);
	}
	return NSPgeometry::rmsd(tList1, tList2);
}

void BRTreeInfoBasic::printPDB(const string& outputFile) {
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), 'A', s[seq[i]]);
		vector<Atom*> aList = nodes[i]->toAtomList(atLib);
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

void BRTreeInfoBasic::printTmpPDB(const string& outputFile){
	RNAChain rc;
	string s = "AUGC";
	char ss[20];
	RnaAtomLib atLib;
	int seqID = 0;
	int atomNum = 0;
	for(int i=0;i<this->seqLen;i++) {

		seqID++;
		sprintf(ss, "%d", seqID);
		RNABase* base = new RNABase(string(ss), 'A', s[seq[i]]);
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

BRTreeInfoBasic::~BRTreeInfoBasic() {
	delete [] seq;
	delete [] connectToDownstream;
	for(int i=0;i<seqLen;i++){
		delete nodes[i];
	}
	delete [] nodes;
}


BRFoldingTreeBasic::BRFoldingTreeBasic(const string& inputFile, EnergyTable* et){

	wcLib = new BaseMoveLibrary("wc");
	wcNbLib = new BaseMoveLibrary("wcNb");
	nwcNbLib = new BaseMoveLibrary("nwcNb");
	riboConnectLib = new BaseMoveLibrary("riboConnect");
	revConnectLib = new BaseMoveLibrary("revConnect");
	loopNbLib = new BaseMoveLibrary("loopNb");
	nwcLibs = new BaseMoveLibrary*[16];
	for(int i=0;i<16;i++){
		nwcLibs[i] = new BaseMoveLibrary("nwc", i);
	}

	//this->et = NULL;
	this->et = et;

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
	this->nodes = new BRNodeBasic*[seqLen];
	this->allBaseBaseE = new double[seqLen*seqLen];
	this->tmpBaseBaseE = new double[seqLen*seqLen];


	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			int pi = i*seqLen+j;
			allBaseBaseE[pi] = 0;
			tmpBaseBaseE[pi] = 0;
		}
	}

	cout << "a" << endl;

	RNABaseName rn;
	for(unsigned int i=0;i<seqLen;i++){
		this->seq[i] = baseList[i]->baseTypeInt;
		this->wcPairPosID[i] = -1;
		this->nwcPairPosID[i] = -1;
		this->fixed[i] = false;
		if(i<seqLen-1 && baseList[i]->connectToNeighbor(*baseList[i+1]))
			this->connectToDownstream[i] = true;
		LocalFrame cs = baseList[i]->getCoordSystem();
		this->nodes[i] = new BRNodeBasic(baseList[i]);

	}

	cout << "b" << endl;
	for(int i=0;i<fixedList.size();i++) {
		this->fixed[fixedList[i]] = true;
		this->nodes[fixedList[i]]->fixed = true;
	}

	for(int i=0;i<seqLen;i++){
		if(!fixed[i])
			flexibleNodes.push_back(nodes[i]);
	}

	this->initTreeInfo = new BRTreeInfoBasic(seqLen, seq, connectToDownstream, nodes, 0.0);

	/*
	 * parse secondary structure
	 */
	cout << "c" << endl;

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
	}

	/*
	 * parse non-Watson-Crick pair
	 */
	cout << "d" << endl;
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
	}

	cout << "e" << endl;
	buildFrom2(nodes[0]);
	buildReverseNodes();


	cout << "f" << endl;

	for(int i=0;i<fixedConnectionList.size();i++) {
		cout << "ct: " << fixedConnectionList[i]->fatherNode->seqID << " " << fixedConnectionList[i]->childNode->seqID << endl;
		fixedConnectionList[i]->childNode->upConnection = fixedConnectionList[i];
	}

	for(int i=0;i<flexibleConnectionList.size();i++) {
		cout << "ct: " << flexibleConnectionList[i]->fatherNode->seqID << " " << flexibleConnectionList[i]->childNode->seqID << endl;
		flexibleConnectionList[i]->childNode->upConnection = flexibleConnectionList[i];
	}


	for(int i=0;i<fixedConnectionList.size();i++) {

		fixedConnectionList[i]->updateChildInfo(nodes);
	}



	for(int i=0;i<flexibleConnectionList.size();i++) {
		flexibleConnectionList[i]->updateChildInfo(nodes);
	}



	for(int i=0;i<seqLen;i++){
		int connectionType = 0;
		if(nodes[i]->upConnection != NULL)
			connectionType = nodes[i]->upConnection->connectionType;
		nodes[i]->updateChildInfo(nodes, seqLen);
	}

	for(int i=0;i<seqLen;i++){
		nodes[i]->updateThreeBaseChildInfo(nodes, seqLen);
	}

	cout << "h" << endl;
	if(nodes[0]->leftChild != NULL) {
		BRConnectionBasic* ct = nodes[0]->leftChild->upConnection;
		updateConnectionChildCs(ct, ct->cm, false);
		updateConnectionChildTmpCs(ct, ct->cm, false);
	}

	if(nodes[0]->midChild != NULL) {
		BRConnectionBasic* ct = nodes[0]->midChild->upConnection;
		updateConnectionChildCs(ct, ct->cm, false);
		updateConnectionChildTmpCs(ct, ct->cm, false);
	}

	cout << "i" << endl;

	updateEnergies();


	for(int i=0;i<seqLen;i++){
		BRNodeBasic* node = nodes[i];
		int groupID1 = -1;
		int groupID2 = -1;
		int groupID3 = -1;
		if(node->midChild != NULL && (!node->midChild->fixed) && node->midChild->midChild != NULL && (!node->midChild->midChild->fixed)){
			for(int j=0;j<fixedGroups.size();j++){
				if(fixedGroups[j].find(node->seqID) != fixedGroups[j].end()){
					groupID1 = j;
				}
				if(fixedGroups[j].find(node->midChild->seqID) != fixedGroups[j].end()){
					groupID2 = j;
				}
				if(fixedGroups[j].find(node->midChild->midChild->seqID) != fixedGroups[j].end()){
					groupID3 = j;
				}
			}

			if(groupID1 == groupID2 && groupID1 >=0)
				continue;
			if(groupID1 == groupID3 && groupID1 >=0)
				continue;
			if(groupID2 == groupID3 && groupID2 >=0)
				continue;

			threeBaseFragList.push_back(node);
			threeBaseFragMoveList.push_back(new ThreeBaseMoveLibrary(node->baseType, node->midChild->baseType, node->midChild->midChild->baseType));
		}
	}

}

void BRFoldingTreeBasic::randInit(){
	for(int i=0;i<flexibleConnectionList.size();i++){
		CsMove cm;
		flexibleConnectionList[i]->cm = flexibleConnectionList[i]->mvLib->getRandomMove(cm);
	}

	if(nodes[0]->leftChild != NULL) {
		BRConnectionBasic* ct = nodes[0]->leftChild->upConnection;
		updateConnectionChildCs(ct, ct->cm, false);
		updateConnectionChildTmpCs(ct, ct->cm, false);
	}
	if(nodes[0]->midChild != NULL) {
		BRConnectionBasic* ct = nodes[0]->midChild->upConnection;
		updateConnectionChildCs(ct, ct->cm, false);
		updateConnectionChildTmpCs(ct, ct->cm, false);
	}
	updateEnergies();
}


void BRFoldingTreeBasic::buildFrom2(BRNodeBasic* node) {

	/*
	 * base first, ribose second
	 */


	/*
	 * build left child
	 */
	int seqID = node->seqID;

	if(wcPairPosID[seqID] > 0){
		BRNodeBasic* wc = nodes[wcPairPosID[seqID]];
		if(wc->father == NULL){
			node->leftChild = wc;
			wc->father = node;
			BRConnectionBasic* ct = new BRConnectionBasic(wcLib, node, wc, seqLen);
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
		BRNodeBasic* nwc = nodes[nwcPairPosID[seqID]];
		if(nwc->father == NULL){
			node->leftChild = nwc;
			nwc->father = node;
			BRConnectionBasic* ct = new BRConnectionBasic(nwcLibs[node->baseType*4+nwc->baseType], node, nwc, seqLen);
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
		BRNodeBasic* nb = nodes[seqID+1];
		// helix neighbor pair
		if(wcPairPosID[node->seqID] > 0 && wcPairPosID[seqID+1] > 0 &&
		   wcPairPosID[node->seqID] == wcPairPosID[seqID+1] + 1) {
			node->midChild = nb;
			nb->father = node;
			BRConnectionBasic* ct = new BRConnectionBasic(wcNbLib, node, nb, seqLen);
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
			BRConnectionBasic* ct = new BRConnectionBasic(nwcNbLib, node, nb, seqLen);
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
			BRConnectionBasic* ct = new BRConnectionBasic(loopNbLib,  node, nb, seqLen);
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

	if(seqID < seqLen-1)
		buildFrom2(nodes[seqID+1]);
}

void BRFoldingTreeBasic::buildReverseNodes() {
	for(int i=seqLen-2;i>0;i--){
		if(nodes[i]->father == NULL){
			cout << "build reverse connection: " << i << endl;
			nodes[i]->father = nodes[i+1];
			nodes[i+1]->reverseChild = nodes[i];
			BRConnectionBasic* ct = new BRConnectionBasic(revConnectLib, nodes[i+1], nodes[i], seqLen);
			bool isFixed = false;
			int idA = nodes[i+1]->seqID;
			int idB = nodes[i]->seqID;
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

void BRFoldingTreeBasic::updateEnergies(){
	bool verbose = false;

	int i,j, pi, pj;
	for(i=0;i<seqLen;i++){
		BRNodeBasic* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNodeBasic* nodeB = nodes[j];
			pi = nodeA->seqID*seqLen+nodeB->seqID;
			pj = nodeB->seqID*seqLen+nodeA->seqID;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				allBaseBaseE[pi] = getNeighborBaseBaseEnergyBasic(nodeA, nodeB, et, verbose);
				allBaseBaseE[pj] = allBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				allBaseBaseE[pi] = getNeighborBaseBaseEnergyBasic(nodeB, nodeA, et, verbose);
				allBaseBaseE[pj] = allBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				allBaseBaseE[pi] = getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 2, et, verbose);
				allBaseBaseE[pj] = allBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				allBaseBaseE[pi] = getNonNeighborBaseBaseEnergyBasic(nodeB, nodeA, 2, et, verbose);
				allBaseBaseE[pj] = allBaseBaseE[pi];
			}
			else {
				allBaseBaseE[pi] = getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 3, et, verbose);
				allBaseBaseE[pj] = allBaseBaseE[pi];
			}
		}
	}

	for(int i=0;i<seqLen;i++){
		for(int j=0;j<seqLen;j++){
			tmpBaseBaseE[i*seqLen+j] = allBaseBaseE[i*seqLen+j];
		}
	}
}

void BRFoldingTreeBasic::updateConnectionChildCs(BRConnectionBasic* ct, CsMove& move, bool verbose){
	BRNodeBasic* childNode = ct->childNode;
	BRNodeBasic* fatherNode = ct->fatherNode;
	ct->cm = move;

	if(verbose){
		cout << "update connection child coordinates: " << fatherNode->seqID << "->" << childNode->seqID << endl;
	}

	int i,j,k;
	childNode->cs1 = fatherNode->cs1 + move;


	if(verbose){
		cout << "update node coordinate: " << childNode->seqID << endl;
	}

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
	}


	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];
		CsMove cmv = cct->cm;
		BRNodeBasic* childNode = cct->childNode;
		BRNodeBasic* fatherNode = cct->fatherNode;
		if(verbose){
			cout << "update node coordinate: " << childNode->seqID << endl;
		}
		childNode->cs1 = fatherNode->cs1 + cmv;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
		}

	}


}

void BRFoldingTreeBasic::updateConnectionChildTmpCs(BRConnectionBasic* ct, CsMove& move, bool verbose){
	BRNodeBasic* childNode = ct->childNode;
	BRNodeBasic* fatherNode = ct->fatherNode;

	if(verbose){
		cout << "update connection child tmp coordinates: " << fatherNode->seqID << "->" << childNode->seqID << endl;
	}

	int i,j,k;
	childNode->tmpCs1 = fatherNode->tmpCs1 + move;

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
	}


	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];

		CsMove cmv = cct->cm;
		BRNodeBasic* childNode = cct->childNode;
		BRNodeBasic* fatherNode = cct->fatherNode;
		if(verbose){
			cout << "update node tmp coordinate: " << childNode->seqID << " connectionType: " << cct->connectionType << endl;
		}
		childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
		}
	}

}

void BRFoldingTreeBasic::clearConnectionChildTmpCs(BRConnectionBasic* ct, bool verbose){
	BRNodeBasic* childNode = ct->childNode;
	BRNodeBasic* fatherNode = ct->fatherNode;

	if(verbose){
		cout << "clear connection child tmp coordinates: " << fatherNode->seqID << childNode->seqID << endl;
	}

	int i,j,k, pi;
	childNode->tmpCs1 = childNode->cs1;

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoordsTmp[i] = childNode->baseAtomCoords[i];
	}



	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];
		CsMove cmv = cct->cm;
		BRNodeBasic* childNode = cct->childNode;
		BRNodeBasic* fatherNode = cct->fatherNode;
		if(verbose){
			cout << "clear node coordinate: " << childNode->seqID << endl;
		}
		childNode->tmpCs1 = childNode->cs1;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = childNode->baseAtomCoords[i];
		}

	}

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
		}
	}
}

void BRFoldingTreeBasic::acceptConnectionChildTmpCs(BRConnectionBasic* ct, CsMove& move, bool verbose){
	ct->cm = move;
	BRNodeBasic* childNode = ct->childNode;
	BRNodeBasic* fatherNode = ct->fatherNode;

	if(verbose){
		cout << "accept connection child tmp coordinates: " << fatherNode->seqID << childNode->seqID << endl;
	}

	int i,j,k, pi;
	childNode->cs1 = childNode->tmpCs1;

	for(i=0;i<childNode->baseAtomNum;i++){
		childNode->baseAtomCoords[i] = childNode->baseAtomCoordsTmp[i];
	}

	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];
		BRNodeBasic* childNode = cct->childNode;

		childNode->cs1 = childNode->tmpCs1;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoords[i] = childNode->baseAtomCoordsTmp[i];
		}
	}

	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
		}
	}

}


double BRFoldingTreeBasic::connectionMutEnergy(BRConnectionBasic* selectConnect, bool verbose){
	double tot = 0;
	int i,j, pi, pj;
	int childNum = selectConnect->childList.size();
	int nonChildNum = selectConnect->nonChildList.size();

	/*
	 * base-base energy
	 */
	for(i=0;i<childNum;i++){
		BRNodeBasic* nodeA = nodes[selectConnect->childList[i]];
		for(j=0;j<nonChildNum;j++){
			BRNodeBasic* nodeB = nodes[selectConnect->nonChildList[j]];
			pi = nodeA->seqID*seqLen + nodeB->seqID;
			pj = nodeB->seqID*seqLen + nodeA->seqID;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else {
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 3, et, verbose);
			}
			tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			tot += tmpBaseBaseE[pi] - allBaseBaseE[pi];
		}
	}
	return tot;
}

double BRFoldingTreeBasic::threeBaseMutEnergy(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose){
	double tot = 0;
	int i,j, pi, pj;

	int groupANum = node->threeBaseFragListA.size();
	int groupBNum = node->threeBaseFragListB.size();
	int groupCNum = node->threeBaseFragListC.size();

	/*
	 * base-base energy
	 */
	for(i=0;i<groupANum;i++){
		BRNodeBasic* nodeA = nodes[node->threeBaseFragListA[i]];
		for(j=0;j<groupBNum;j++){
			BRNodeBasic* nodeB = nodes[node->threeBaseFragListB[j]];
			pi = nodeA->seqID*seqLen + nodeB->seqID;
			pj = nodeB->seqID*seqLen + nodeA->seqID;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else {
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 3, et, verbose);
			}
			tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			tot += tmpBaseBaseE[pi] - allBaseBaseE[pi];
		}
	}

	for(i=0;i<groupBNum;i++){
		BRNodeBasic* nodeA = nodes[node->threeBaseFragListB[i]];
		for(j=0;j<groupCNum;j++){
			BRNodeBasic* nodeB = nodes[node->threeBaseFragListC[j]];
			pi = nodeA->seqID*seqLen + nodeB->seqID;
			pj = nodeB->seqID*seqLen + nodeA->seqID;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else {
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 3, et, verbose);
			}
			tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			tot += tmpBaseBaseE[pi] - allBaseBaseE[pi];
		}
	}

	for(i=0;i<groupANum;i++){
		BRNodeBasic* nodeA = nodes[node->threeBaseFragListA[i]];
		for(j=0;j<groupCNum;j++){
			BRNodeBasic* nodeB = nodes[node->threeBaseFragListC[j]];
			pi = nodeA->seqID*seqLen + nodeB->seqID;
			pj = nodeB->seqID*seqLen + nodeA->seqID;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tmpBaseBaseE[pi] = getNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, 2, et, verbose);
				tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			}
			else {
				tmpBaseBaseE[pi] = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 3, et, verbose);
			}
			tmpBaseBaseE[pj] = tmpBaseBaseE[pi];
			tot += tmpBaseBaseE[pi] - allBaseBaseE[pi];
		}
	}
	return tot;
}

double BRFoldingTreeBasic::totalEnergy(bool verbose){
	double tot = 0;
	int i,j;
	for(i=0;i<seqLen;i++){
		BRNodeBasic* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNodeBasic* nodeB = nodes[j];
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tot += getNeighborBaseBaseEnergyBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tot += getNeighborBaseBaseEnergyBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tot += getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 2, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tot += getNonNeighborBaseBaseEnergyBasic(nodeB, nodeA, 2, et, verbose);
			}
			else {
				tot += getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 3, et, verbose);
			}
		}
	}

	return tot;
}

void BRFoldingTreeBasic::updateThreeBaseFragChildCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose){
	BRNodeBasic* child = node->midChild;
	BRNodeBasic* grandChild = child->midChild;

	child->upConnection->cm = move1;
	grandChild->upConnection->cm = move2;

	child->cs1 = node->cs1 + move1;
	grandChild->cs1 = child->cs1+move2;

	int i,j,k;

	for(i=0;i<child->baseAtomNum;i++){
		child->baseAtomCoords[i] = local2global(child->cs1, child->atomCoordLocal[i]);
	}

	for(i=0;i<grandChild->baseAtomNum;i++){
		grandChild->baseAtomCoords[i] = local2global(grandChild->cs1, grandChild->atomCoordLocal[i]);
	}

	BRConnectionBasic* ct = child->upConnection;
	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];
		CsMove cmv = cct->cm;
		BRNodeBasic* childNode = cct->childNode;
		BRNodeBasic* fatherNode = cct->fatherNode;
		if(verbose){
			cout << "update node coordinate: " << childNode->seqID << endl;
		}
		childNode->cs1 = fatherNode->cs1 + cmv;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
		}
	}
	if(child->leftChild != NULL){
		ct = child->leftChild->upConnection;
		child->leftChild->cs1 = child->cs1 + ct->cm;
		for(i=0;i<child->leftChild->baseAtomNum;i++){
			child->leftChild->baseAtomCoords[i] = local2global(child->leftChild->cs1, child->leftChild->atomCoordLocal[i]);
		}

		for(j=0;j<ct->childConnectionList.size();j++){
			BRConnectionBasic* cct = ct->childConnectionList[j];
			CsMove cmv = cct->cm;
			BRNodeBasic* childNode = cct->childNode;
			BRNodeBasic* fatherNode = cct->fatherNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->cs1 = fatherNode->cs1 + cmv;

			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
			}
		}
	}

	if(child->reverseChild != NULL) {
		ct = child->reverseChild->upConnection;
		child->reverseChild->cs1 = child->cs1 + ct->cm;
		for(i=0;i<child->reverseChild->baseAtomNum;i++){
			child->reverseChild->baseAtomCoords[i] = local2global(child->reverseChild->cs1, child->reverseChild->atomCoordLocal[i]);
		}

		for(j=0;j<ct->childConnectionList.size();j++){
			BRConnectionBasic* cct = ct->childConnectionList[j];
			CsMove cmv = cct->cm;
			BRNodeBasic* childNode = cct->childNode;
			BRNodeBasic* fatherNode = cct->fatherNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->cs1 = fatherNode->cs1 + cmv;

			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoords[i] = local2global(childNode->cs1, childNode->atomCoordLocal[i]);
			}
		}
	}

}

void BRFoldingTreeBasic::updateThreeBaseFragChildTmpCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose){
	BRNodeBasic* child = node->midChild;
	BRNodeBasic* grandChild = child->midChild;


	child->tmpCs1 = node->tmpCs1 + move1;
	grandChild->tmpCs1 = child->tmpCs1+move2;

	int i,j,k;

	for(i=0;i<child->baseAtomNum;i++){
		child->baseAtomCoordsTmp[i] = local2global(child->tmpCs1, child->atomCoordLocal[i]);
	}

	for(i=0;i<grandChild->baseAtomNum;i++){
		grandChild->baseAtomCoordsTmp[i] = local2global(grandChild->tmpCs1, grandChild->atomCoordLocal[i]);
	}

	BRConnectionBasic* ct = grandChild->upConnection;
	for(j=0;j<ct->childConnectionList.size();j++){
		BRConnectionBasic* cct = ct->childConnectionList[j];
		CsMove cmv = cct->cm;
		BRNodeBasic* childNode = cct->childNode;
		BRNodeBasic* fatherNode = cct->fatherNode;
		if(verbose){
			cout << "update node coordinate: " << childNode->seqID << endl;
		}
		childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
		}
	}

	if(child->leftChild != NULL){
		ct = child->leftChild->upConnection;
		child->leftChild->tmpCs1 = child->tmpCs1 + ct->cm;
		for(i=0;i<child->leftChild->baseAtomNum;i++){
			child->leftChild->baseAtomCoordsTmp[i] = local2global(child->leftChild->tmpCs1, child->leftChild->atomCoordLocal[i]);
		}

		for(j=0;j<ct->childConnectionList.size();j++){
			BRConnectionBasic* cct = ct->childConnectionList[j];
			CsMove cmv = cct->cm;
			BRNodeBasic* childNode = cct->childNode;
			BRNodeBasic* fatherNode = cct->fatherNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;

			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
			}
		}
	}

	if(child->reverseChild != NULL) {
		ct = child->reverseChild->upConnection;
		child->reverseChild->tmpCs1 = child->tmpCs1 + ct->cm;
		for(i=0;i<child->reverseChild->baseAtomNum;i++){
			child->reverseChild->baseAtomCoordsTmp[i] = local2global(child->reverseChild->tmpCs1, child->reverseChild->atomCoordLocal[i]);
		}

		for(j=0;j<ct->childConnectionList.size();j++){
			BRConnectionBasic* cct = ct->childConnectionList[j];
			CsMove cmv = cct->cm;
			BRNodeBasic* childNode = cct->childNode;
			BRNodeBasic* fatherNode = cct->fatherNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->tmpCs1 = fatherNode->tmpCs1 + cmv;

			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoordsTmp[i] = local2global(childNode->tmpCs1, childNode->atomCoordLocal[i]);
			}
		}
	}
}

void BRFoldingTreeBasic::clearThreeBaseFragChildTmpCs(BRNodeBasic* node, bool verbose){
	BRNodeBasic* child = node->midChild;
	BRNodeBasic* grandChild = child->midChild;

	child->tmpCs1 = child->cs1;
	grandChild->tmpCs1 = grandChild->cs1;

	int i,j,k;

	for(i=0;i<child->baseAtomNum;i++){
		child->baseAtomCoordsTmp[i] = child->baseAtomCoords[i];
	}

	for(i=0;i<grandChild->baseAtomNum;i++){
		grandChild->baseAtomCoordsTmp[i] = grandChild->baseAtomCoords[i];
	}

	BRConnectionBasic* ct = child->upConnection;
	for(j=0;j<ct->childConnectionList.size();j++){
		BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
		if(verbose){
			cout << "update node coordinate: " << childNode->seqID << endl;
		}
		childNode->tmpCs1 = childNode->cs1;
		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoordsTmp[i] = childNode->baseAtomCoords[i];
		}
	}

	if(child->leftChild != NULL){
		ct = child->leftChild->upConnection;
		child->leftChild->tmpCs1 = child->leftChild->cs1;
		for(i=0;i<child->leftChild->baseAtomNum;i++){
			child->leftChild->baseAtomCoordsTmp[i] = child->leftChild->baseAtomCoords[i];
		}
		for(j=0;j<ct->childConnectionList.size();j++){
			BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->tmpCs1 = childNode->cs1;
			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoordsTmp[i] = childNode->baseAtomCoords[i];
			}
		}
	}

	if(child->reverseChild != NULL) {
		ct = child->reverseChild->upConnection;
		child->reverseChild->tmpCs1 = child->reverseChild->cs1;
		for(i=0;i<child->reverseChild->baseAtomNum;i++){
			child->reverseChild->baseAtomCoordsTmp[i] = child->reverseChild->baseAtomCoords[i];
		}
		for(j=0;j<ct->childConnectionList.size();j++){
			BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->tmpCs1 = childNode->cs1;
			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoordsTmp[i] = childNode->baseAtomCoords[i];
			}
		}
	}

	int pi;
	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			tmpBaseBaseE[pi] = allBaseBaseE[pi];
		}
	}
}

void BRFoldingTreeBasic::acceptThreeBaseFragChildTmpCs(BRNodeBasic* node, CsMove& move1, CsMove& move2, bool verbose){
	BRNodeBasic* child = node->midChild;
	BRNodeBasic* grandChild = child->midChild;

	child->upConnection->cm = move1;
	grandChild->upConnection->cm = move2;


	child->cs1 = child->tmpCs1;
	grandChild->cs1 = grandChild->tmpCs1;

	int i,j,k;

	for(i=0;i<child->baseAtomNum;i++){
		child->baseAtomCoords[i] = child->baseAtomCoordsTmp[i];
	}

	for(i=0;i<grandChild->baseAtomNum;i++){
		grandChild->baseAtomCoords[i] = grandChild->baseAtomCoordsTmp[i];
	}

	BRConnectionBasic* ct = child->upConnection;
	for(j=0;j<ct->childConnectionList.size();j++){
		BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
		if(verbose){
			cout << "update node coordinate: " << childNode->seqID << endl;
		}
		childNode->cs1 = childNode->tmpCs1;

		for(i=0;i<childNode->baseAtomNum;i++){
			childNode->baseAtomCoords[i] = childNode->baseAtomCoordsTmp[i];
		}
	}

	if(child->leftChild != NULL){
		ct = child->leftChild->upConnection;
		child->leftChild->cs1 = child->leftChild->tmpCs1;
		for(i=0;i<child->leftChild->baseAtomNum;i++){
			child->leftChild->baseAtomCoords[i] = child->leftChild->baseAtomCoordsTmp[i];
		}
		for(j=0;j<ct->childConnectionList.size();j++){
			BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->cs1 = childNode->tmpCs1;
			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoords[i] = childNode->baseAtomCoordsTmp[i];
			}
		}
	}

	if(child->reverseChild != NULL) {
		ct = child->reverseChild->upConnection;
		child->reverseChild->cs1 = child->reverseChild->tmpCs1;
		for(i=0;i<child->reverseChild->baseAtomNum;i++){
			child->reverseChild->baseAtomCoords[i] = child->reverseChild->baseAtomCoordsTmp[i];
		}
		for(j=0;j<ct->childConnectionList.size();j++){
			BRNodeBasic* childNode = ct->childConnectionList[j]->childNode;
			if(verbose){
				cout << "update node coordinate: " << childNode->seqID << endl;
			}
			childNode->cs1 = childNode->tmpCs1;
			for(i=0;i<childNode->baseAtomNum;i++){
				childNode->baseAtomCoords[i] = childNode->baseAtomCoordsTmp[i];
			}
		}
	}

	int pi;
	for(i=0;i<seqLen;i++){
		for(j=0;j<seqLen;j++){
			pi = i*seqLen+j;
			allBaseBaseE[pi] = tmpBaseBaseE[pi];
		}
	}
}

void BRFoldingTreeBasic::printDetailEnergy(){
	double e = 0;


	bool verbose = false;

	for(int i=0;i<seqLen;i++){
		BRNodeBasic* nodeA = nodes[i];
		for(int j=i+1;j<seqLen;j++){
			BRNodeBasic* nodeB = nodes[j];
			double tot = 0;
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				tot += getNeighborBaseBaseEnergyBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				tot += getNeighborBaseBaseEnergyBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				tot += getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 2, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				tot += getNonNeighborBaseBaseEnergyBasic(nodeB, nodeA, 2, et, verbose);
			}
			else {
				tot += getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 3, et, verbose);
			}



			printf("Pair: %2d-%-2d ene: %9.3f\n",i,j, tot);
			e += tot;
		}
	}

	//printf("total energy: %9.3f\n",  e);
}

void BRFoldingTreeBasic::checkCoordinate(){
	for(int i=0;i<seqLen;i++){
		BRNodeBasic* node = nodes[i];
		if(node->baseAtomCoords[0].distance(node->baseAtomCoordsTmp[0]) > 0.0001){
			cout << "COORD DIFF: " << node->seqID << endl;
		}
	}
}

void BRFoldingTreeBasic::checkTotalEnergy(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e;
	for(i=0;i<seqLen;i++){
		BRNodeBasic* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNodeBasic* nodeB = nodes[j];
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				e = getNeighborBaseBaseEnergyBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				e = getNeighborBaseBaseEnergyBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				e = getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 2, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				e = getNonNeighborBaseBaseEnergyBasic(nodeB, nodeA, 2, et, verbose);
			}
			else {
				e = getNonNeighborBaseBaseEnergyBasic(nodeA, nodeB, 3, et, verbose);
			}
			if(abs(allBaseBaseE[i*seqLen+j] - e) > 0.001){
				cout << "DIFF base-base: " << i << " " << j << " " << allBaseBaseE[i*seqLen+j] << " " << e <<  endl;
			}
		}
	}



}

void BRFoldingTreeBasic::checkTmpTotalEnergy(){
	bool verbose = false;
	double tot = 0;
	int i,j;
	double e;
	for(i=0;i<seqLen;i++){
		BRNodeBasic* nodeA = nodes[i];
		for(j=i+1;j<seqLen;j++){
			BRNodeBasic* nodeB = nodes[j];
			if(nodeA->seqID == nodeB->seqID -1 && connectToDownstream[nodeA->seqID]){
				e = getNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 1 && connectToDownstream[nodeB->seqID]){
				e = getNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID - 2 && connectToDownstream[nodeA->seqID] && connectToDownstream[nodeA->seqID+1]){
				e = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 2, et, verbose);
			}
			else if(nodeA->seqID == nodeB->seqID + 2 && connectToDownstream[nodeB->seqID] && connectToDownstream[nodeB->seqID+1]){
				e = getNonNeighborBaseBaseEnergyTmpBasic(nodeB, nodeA, 2, et, verbose);
			}
			else {
				e = getNonNeighborBaseBaseEnergyTmpBasic(nodeA, nodeB, 3, et, verbose);
			}
			if(abs(tmpBaseBaseE[i*seqLen+j] - e) > 0.001){
				cout << "DIFF base-base: " << i << " " << j << " " << allBaseBaseE[i*seqLen+j] << " " << tmpBaseBaseE[i*seqLen+j] << " " << e <<  endl;
			}
		}
	}

}

BRTreeInfoBasic* BRFoldingTreeBasic::getTreeInfo() {
	BRTreeInfoBasic* info = new BRTreeInfoBasic(seqLen, seq, connectToDownstream, nodes, totalEnergy(false));
	return info;
}


void BRFoldingTreeBasic::printTree(int index) {
	cout << "Node: " << index;
	BRNodeBasic* node = nodes[index];
	if(node->leftChild != NULL) {
		BRConnectionBasic* ct = node->leftChild->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout<< " Left: " << node->leftChild->seqID << " " << ctType;
	}
	if(node->midChild != NULL) {
		BRConnectionBasic* ct = node->midChild->upConnection;
		string ctType = "flex";
		if(ct->fixed)
			ctType = "fixed";
		cout << " Mid: " << node->midChild->seqID << " " << ctType;
	}

	if(node->reverseChild != NULL) {
		BRNodeBasic* rn = node->reverseChild;
		BRConnectionBasic* ct = rn->upConnection;
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
	if(node->reverseChild != NULL){
		printTree(node->reverseChild->seqID);
	}
}

void BRFoldingTreeBasic::printConnections() {
	for(int i=0;i<fixedConnectionList.size();i++) {
		BRConnectionBasic* ct = fixedConnectionList[i];
		cout << ct->fatherNode->seqID << " " << ct->childNode->seqID << " " << ct->connectionType << endl;

		fixedConnectionList[i]->printPartition();
		fixedConnectionList[i]->printChildConnections();
	}
	for(int i=0;i<flexibleConnectionList.size();i++) {
		BRConnectionBasic* ct = flexibleConnectionList[i];
		cout << ct->fatherNode->seqID << " " << ct->childNode->seqID << " " << ct->connectionType << endl;
		flexibleConnectionList[i]->printPartition();
		flexibleConnectionList[i]->printChildConnections();
	}
}

BRFoldingTreeBasic::~BRFoldingTreeBasic() {
	delete [] seq;
	delete [] wcPairPosID;
	delete [] nwcPairPosID;
	delete [] fixed;
	delete [] connectToDownstream;
	if(seqLen > 0) {
		for(int i=0;i<seqLen;i++)
			delete nodes[i];
	}
	delete [] nodes;
}

} /* namespace NSPpred */
