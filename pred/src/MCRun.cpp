/*
 * MCRun.cpp
 *
 *  Created on: Jul 30, 2019
 *      Author: s2982206
 */

#include <pred/MCRun.h>

namespace NSPpred {

void MCRun::optimizeFromInit(const string& keyFile, const string& outFilePrefix, int startID){

	cout << "opt from init: " << endl;
	ifstream f;
	f.open(keyFile, ios::in);
	if(!f.is_open()){
		cout << "fail to open " << keyFile << endl;
		exit(1);
	}

	string line;
	char xx[200];
	int outID = startID;
	int stepNum = this->ft->flexibleConnectionList.size()*5000;
	bool verbose = false;



	vector<BRConnection*> flexConnectList = ft->flexibleConnectionList;
	vector<BRNode*> flexNodes = ft->flexibleNodes;
	int fc = flexConnectList.size();
	int fn = flexNodes.size();
	int tn = ft->seqLen;
	int freeNodeNum = ft->freeNodeIDs.size();
	int fcn = fc + fn;
	int randIndex;

	BRConnection* ct;
	BRNode* node;

	BaseMoveLibrary mvLib;

	CsMove mvMut;
	int randType;

	BaseRotamer* rotMut;
	F2Fragment* frag;
	F3Fragment* f3Frag;


	while(getline(f, line)){

		double connectWT = 0.5;
		double clashWT = 0.5;
		double breakCTWT = 0.1;
		sprintf(xx, "%s-%d.pdb", outFilePrefix.c_str(), outID);
		string outFileName = string(xx);
		outID ++;
		cout << "init from key: " << line << endl;
		this->ft->initFromKey(line);
		double curEne = ft->totalEnergy(1.0, 1.0, 1.0, 0.0, verbose);
		double lastEne = curEne;
		pair<double,double> mutE;

		for(double T=0.2;T>0.02;T=T*0.9){
			int acNum = 0;
			int ac1 = 0;
			int ac2 = 0;
			int ac3 = 0;
			int ac4 = 0;
			int ac5 = 0;
			int ac6 = 0;
			int indexX;
			double phoMutE;
			for(int k=0;k<stepNum;k++){
				randIndex = rand()%fcn;
				if(randIndex < fc){
					ct = flexConnectList[randIndex];
					randType = rand()%10;
					if(randType < 3) {
						if(ct->f2Lib == NULL) continue;
						frag = ct->f2Lib->getRandomFrag();
						ft->updateF2ChildTmpCs(ct, frag, false);
						mutE = ft->f2MutEnergy( ct, breakCTWT, connectWT, clashWT, 0.0, false);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF2ChildTmpCs(ct, verbose);
							acNum ++;
							ac1 ++;
						}
						else {
							ft->clearF2ChildTmpCs(ct, verbose);
						}
					}
					else if(randType < 6) {
						if(ct->f2Lib == NULL) continue;
						frag = ct->f2Lib->getRandomFragLevel2(ct->f2Frag->level1ID);
						ft->updateF2ChildTmpCs(ct, frag, false);
						mutE = ft->f2MutEnergy(ct, breakCTWT,connectWT,clashWT, 0.0, false);


						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF2ChildTmpCs(ct, verbose);
							acNum ++;
							ac2 ++;
						}
						else {
							ft->clearF2ChildTmpCs(ct, verbose);
						}
					}
					else if(randType < 7) {
						mvMut = mvLib.getRandomMove2(ct->cm);
						ft->updateCtChildTmpCs(ct, mvMut, verbose);
						mutE = ft->ctMutEnergy(ct,breakCTWT,connectWT,clashWT,  0.0,verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptCtChildTmpCs(ct, verbose);
							acNum ++;
							ac3 ++;
						}
						else {
							ft->clearCtChildTmpCs(ct, verbose);
						}
					}
					else if(ct->hasThreeBaseFragment){
						f3Frag = ct->f3Lib->getRandomFrag();
						ft->updateF3ChildTmpCs(ct, f3Frag, verbose);

						mutE = ft->f3MutEnergy(ct,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF3ChildTmpCs(ct, verbose);
							acNum ++;
							ac4 ++;
						}
						else {
							ft->clearF3ChildTmpCs(ct, verbose);
						}
					}

				}
				else {
					randType = rand()%10;
					if(randType < 3){
						node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
						mvMut = mvLib.getRandomMove2();
						ft->updateSingleBaseCoordTmp(node, mvMut, false);
						mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptSingleBaseCoordTmp(node, verbose);
							acNum ++;
							ac5 ++;
						}
						else {
							ft->clearSingleBaseCoordTmp(node, verbose);
						}
					}
					else {
						node = ft->nodes[rand()%tn];
						rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
						ft->updateBaseRotamerTmp(node, rotMut, verbose);
						mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptBaseRotamerTmp(node, verbose);
							acNum ++;
							ac6 ++;
							string key = ft->toCtKey();
						}
						else {
							ft->clearBaseRotamerTmp(node, verbose);
						}
					}

				}

			}
			BRTreeInfo* x = ft->getTreeInfo();
			printf("T=%9.5f WT=%6.5f ac1=%5d ac2=%5d ac3=%5d ac4=%5d ac5=%5d ac6=%5d addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,connectWT, ac1,ac2,ac3,ac4,ac5,ac6,curEne,ft->totalEnergy(1.0, 1.0,1.0, 0.0,false), x->rmsd(init));
			delete x;

			if(connectWT < 1.0) connectWT = connectWT*1.1;
			if(clashWT < 1.0) clashWT = clashWT*1.1;
			if(breakCTWT < 1.0) breakCTWT = breakCTWT*1.1;
		}

		BRTreeInfo* x = ft->getTreeInfo();
		x->printPDB(outFileName);
		delete x;
	}
}

void MCRun::generateDecoysRandInit(const string& output){

	ofstream of;
	of.open(output, ios::out);
	if(!of.is_open()){
		cout << "fail to open " << output << endl;
		exit(1);
	}

	double T0 = 3.0;
	double T1 = 0.1;
	int stepNum = ft->flexibleConnectionList.size()*4000;
	double eneCutoff = 0;
	ft->randInit();

	map<string, double> keyToEnergy;
	map<string, string> keyToDetail;

	pair<double,double> mutE;


	vector<BRConnection*> flexConnectList = ft->flexibleConnectionList;
	vector<BRNode*> flexNodes = ft->flexibleNodes;
	int fc = flexConnectList.size();
	int fn = flexNodes.size();
	int tn = ft->seqLen;

	int freeNodeNum = ft->freeNodeIDs.size();

	int fcn = fc + fn;
	int randIndex;
	BRConnection* ct;
	BRNode* node;

	BaseMoveLibrary mvLib;

	CsMove mvMut;
	int randType;

	BaseRotamer* rotMut;
	F2Fragment* frag;
	F3Fragment* f3Frag;

	bool verbose = false;

	cout << "cal total energy: " << endl;
	double curEne = ft->totalEnergy(1.0, 1.0, 1.0, 0.0, verbose);
	cout << "total energy: " << curEne << endl;
	double lastEne = curEne;


	for(int i=0;i<10;i++){
		cout << "round: " << i << endl;
		double clashWT = pow(1.2, rand()%13) * 0.06;
		double connectWT = pow(1.2, rand()%13) * 0.1;
		double breakCTWT = pow(1.2, rand()%13) * 0.002;
		for(double T=T0;T>T1;T=T*0.9){
			for(int k=0;k<stepNum;k++){
				int acNum = 0;
				int ac1 = 0;
				int ac2 = 0;
				int ac3 = 0;
				int ac4 = 0;
				int ac5 = 0;
				int ac6 = 0;
				int indexX;
				randIndex = rand()%fcn;
				if(randIndex < fc){
					ct = flexConnectList[randIndex];
					randType = rand()%10;
					if(randType < 3) {
						if(ct->f2Lib == NULL) continue;
						frag = ct->f2Lib->getRandomFrag();
						ft->updateF2ChildTmpCs(ct, frag, false);
						mutE = ft->f2MutEnergy( ct, breakCTWT,connectWT, clashWT, 0.0,false);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF2ChildTmpCs(ct, verbose);
							string key = ft->toCtKey();
							if(curEne < eneCutoff && keyToEnergy.find(key) != keyToEnergy.end()){
								if(keyToEnergy[key] > curEne){
									keyToEnergy[key] = curEne;
									keyToDetail[key] = ft->toCtDetailString();
								}
							}
							else if(curEne < eneCutoff){
								keyToEnergy[key] = curEne;
								keyToDetail[key] = ft->toCtDetailString();
							}
							acNum ++;
							ac1 ++;
						}
						else {
							ft->clearF2ChildTmpCs(ct, verbose);
						}
					}
					else if(randType < 6) {
						if(ct->f2Lib == NULL) continue;
						frag = ct->f2Lib->getRandomFragLevel2(ct->f2Frag->level1ID);
						ft->updateF2ChildTmpCs(ct, frag, false);
						mutE = ft->f2MutEnergy(ct,breakCTWT, connectWT,clashWT, 0.0,false);


						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF2ChildTmpCs(ct, verbose);
							acNum ++;
							ac2 ++;
							string key = ft->toCtKey();
							if(curEne < eneCutoff && keyToEnergy.find(key) != keyToEnergy.end()){
								if(keyToEnergy[key] > curEne){
									keyToEnergy[key] = curEne;
									keyToDetail[key] = ft->toCtDetailString();
								}
							}
							else if(curEne < eneCutoff){
								keyToEnergy[key] = curEne;
								keyToDetail[key] = ft->toCtDetailString();
							}
						}
						else {
							ft->clearF2ChildTmpCs(ct, verbose);
						}
					}
					else if(randType < 7) {
						if(ct->ctType == "wc" || ct->ctType == "nwc") continue;
						mvMut = mvLib.getRandomMove2(ct->cm);
						ft->updateCtChildTmpCs(ct, mvMut, verbose);
						mutE = ft->ctMutEnergy(ct,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptCtChildTmpCs(ct, verbose);
							acNum ++;
							ac3 ++;
							string key = ft->toCtKey();
							if(curEne < eneCutoff && keyToEnergy.find(key) != keyToEnergy.end()){
								if(keyToEnergy[key] > curEne){
									keyToEnergy[key] = curEne;
									keyToDetail[key] = ft->toCtDetailString();
								}
							}
							else if(curEne < eneCutoff){
								keyToEnergy[key] = curEne;
								keyToDetail[key] = ft->toCtDetailString();
							}
						}
						else {
							ft->clearCtChildTmpCs(ct, verbose);
						}
					}
					else if(ct->hasThreeBaseFragment){
						f3Frag = ct->f3Lib->getRandomFrag();
						ft->updateF3ChildTmpCs(ct, f3Frag, verbose);
						mutE = ft->f3MutEnergy(ct,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptF3ChildTmpCs(ct, verbose);
							acNum ++;
							ac4 ++;
						}
						else {
							ft->clearF3ChildTmpCs(ct, verbose);
						}
					}
				}
				else {
					randType = rand()%10;
					if(randType < 3){
						continue;
						node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
						mvMut = mvLib.getRandomMove2();
						ft->updateSingleBaseCoordTmp(node, mvMut, false);
						mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT, 0.0, verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptSingleBaseCoordTmp(node, verbose);
							acNum ++;
							ac5 ++;
							string key = ft->toCtKey();
							if(curEne < eneCutoff && keyToEnergy.find(key) != keyToEnergy.end()){
								if(keyToEnergy[key] > curEne){
									keyToEnergy[key] = curEne;
									keyToDetail[key] = ft->toCtDetailString();
								}
							}
							else if(curEne < eneCutoff){
								keyToEnergy[key] = curEne;
								keyToDetail[key] = ft->toCtDetailString();
							}
						}
						else {
							ft->clearSingleBaseCoordTmp(node, verbose);
						}
					}
					else {
						node = ft->nodes[rand()%tn];
						rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
						ft->updateBaseRotamerTmp(node, rotMut, verbose);
						mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT,  0.0,verbose);
						if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
							curEne += mutE.first;
							ft->acceptBaseRotamerTmp(node, verbose);
							acNum ++;
							ac6 ++;
							string key = ft->toCtKey();
							if(curEne < eneCutoff && keyToEnergy.find(key) != keyToEnergy.end()){
								if(keyToEnergy[key] > curEne){
									keyToEnergy[key] = curEne;
									keyToDetail[key] = ft->toCtDetailString();
								}
							}
							else if(curEne < eneCutoff){
								keyToEnergy[key] = curEne;
								keyToDetail[key] = ft->toCtDetailString();
							}
						}
						else {
							ft->clearBaseRotamerTmp(node, verbose);
						}
					}
				}

			}

			BRTreeInfo* x = ft->getTreeInfo();
			//ft->printDetailEnergy();
			printf("T=%9.5f WT=%6.5f addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,connectWT,curEne,ft->totalEnergy(1.0, 1.0,1.0, 0.0,false), x->rmsd(init));
			delete x;


			if(clashWT < 1.0)
				clashWT = clashWT*1.07;
			if(connectWT < 1.0)
				connectWT = connectWT*1.07;
			if(breakCTWT < 1.0)
				breakCTWT = breakCTWT*1.1;
		}


	}

	map<string,double>::iterator it;
	for(it = keyToEnergy.begin();it != keyToEnergy.end();it++){
		string dt = keyToDetail[it->first];
		of << it->first << " " << it->second << " " << dt << endl;
	}
	of.close();

}


void MCRun::simpleMC(const string& outPDB, bool printTraj){

	clock_t start = clock();
	ofstream of;
	of.open(outPDB.c_str(), ios::out);

	if(printTraj)
		ft->initTreeInfo->printPDB(of, 0);

	double T0= 2.0;
	double T1 = 0.1;

	int stepNum = 0;

	vector<int> ctChooseIndex;

	for(int i=0;i<ft->flexibleConnectionList.size();i++){
		BRConnection* ct = ft->flexibleConnectionList[i];
		if(ct->ctType == "wc" || ct->ctType == "wcNb")
		{
			stepNum += 400;
			ctChooseIndex.push_back(i);
		}
		else if(ct->ctType == "nwc" || ct->ctType == "nwcNb") {
			stepNum += 2000;
			for(int k=0;k<5;k++){
				ctChooseIndex.push_back(i);
			}
		}
		else {
			stepNum += 4000;
			for(int k=0;k<10;k++){
				ctChooseIndex.push_back(i);
			}
		}
	}

	int ctChooseNum = ctChooseIndex.size();

	if(ctChooseNum == 0) return;

	bool verbose = false;

	ft->randInit();

	double clashWT;
	double connectWT;
	double breakCTWT;

	double shift = ft->et->para.initShift;

	pair<double,double> mutE;


	int fc = ft->flexibleConnectionList.size();
	int fn = ft->flexibleNodes.size();
	int tn = ft->seqLen;
	int rfn = ft->riboFlexibleNodes.size();
	int freeNodeNum = ft->freeNodeIDs.size();

	cout << "flexible connections: " << endl;
	for(int i=0;i<fc;i++){
		cout << ft->flexibleConnectionList[i]->fatherNode->seqID << "->" << ft->flexibleConnectionList[i]->childNode->seqID << " ";
	}
	cout << endl;

	cout << "flexible nodes: " << endl;
	for(int i=0;i<fn;i++){
		cout << ft->flexibleNodes[i]->seqID << " ";
	}
	cout << endl;


	cout << "free nodes: " << endl;
	for(int i=0;i<freeNodeNum;i++){
		cout << ft->freeNodeIDs[i] << " ";
	}
	cout << endl;

	cout << "connection: " << endl;
	ft->printConnections();

	int fcn = fc + fn/2;
	int randIndex;
	BRConnection* ct;
	BRNode* node;

	BaseMoveLibrary mvLib;

	CsMove mvMut;
	int randType;

	BaseRotamer* rotMut;
	F2Fragment* frag;
	F3Fragment* f3Frag;

	int outFreq = ft->et->para.outFreq;
	int modelID = 1;


	ft->updateEnergies(shift);
	double curEne = ft->totalEnergy(1.0, 1.0, 1.0, shift, verbose);
	printf("total energy: %12.3f\n", curEne);

	double lastEne = curEne;


	double ctRate = ft->et->para.connectWTFactor;
	double clashIncreaseRate = ft->et->para.clashWTFactor;

	clashWT = ft->et->para.initClashWT;
	connectWT = ft->et->para.initConnectWT;
	breakCTWT = 0.6*ft->et->para.initConnectWT;

	T0 = ft->et->para.T0;
	double anneal = ft->et->para.anneal;

	for(double T=T0;T>0.1;T=T*anneal){

		shift = shift - ft->et->para.dShift;
		if(shift < 0) shift = 0;
		ft->updateEnergies(shift);
		curEne = ft->totalEnergy(1.0, 1.0, 1.0, shift, verbose);

		if(clashWT < 1.0)
			clashWT = clashWT*clashIncreaseRate;
		if(connectWT < 1.0)
			connectWT = connectWT*ctRate;
		if(breakCTWT < 1.0)
			breakCTWT = breakCTWT*ctRate;
		if(clashWT > 1.0) clashWT = 1.0;
		if(connectWT > 1.0) connectWT = 1.0;
		if(breakCTWT > 1.0) breakCTWT = 1.0;

		int acNum = 0;
		int ac1 = 0;
		int ac2 = 0;
		int ac3 = 0;
		int ac4 = 0;
		int ac5 = 0;
		int ac6 = 0;
		int ac7 = 0;
		int indexX;
		double phoMutE;

		for(int k=0;k<stepNum;k++){
			//cout << "k: " << k << endl;
			if(k % outFreq == 0 && printTraj)
			{
				BRTreeInfo* x = ft->getTreeInfo();
				x->printPDB(of, modelID);
				modelID++;
				delete x;
			}

			randIndex = rand()%fcn;
			if(randIndex < fc){

				ct = ft->flexibleConnectionList[ctChooseIndex[rand()%ctChooseNum]];
				randType = rand()%10;
				if(randType < 4) {
					if(ct->f2Lib == NULL) continue;
					if(ct->ctType == "loopNb" && rand()%3 != 0 && ft->et->para.loopRiboConnectMove){
						CsMove cm = ct->fatherNode->rot->mv12 + ft->fragLib->ribLib->getRandomMove() + ct->childNode->rot->mv31;
						ft->updateCtChildTmpCs(ct, cm, false);
						mutE = ft->ctMutEnergy(ct, breakCTWT, connectWT, clashWT, shift, false);
					}
					else {

						frag = ct->f2Lib->getRandomFrag();
						ft->updateF2ChildTmpCs(ct, frag, false);
						mutE = ft->f2MutEnergy( ct, breakCTWT, connectWT, clashWT,shift,false);
					}

					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF2ChildTmpCs(ct, verbose);
						acNum ++;
						ac1 ++;
					}
					else {
						ft->clearF2ChildTmpCs(ct, verbose);
					}
				}
				else if(randType < 6) {

					if(ct->f2Lib == NULL) continue;
					frag = ct->f2Lib->getRandomFragLevel2(ct->f2Frag->level1ID);
					ft->updateF2ChildTmpCs(ct, frag, false);
					mutE = ft->f2MutEnergy(ct, breakCTWT,connectWT,clashWT,shift,false);


					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF2ChildTmpCs(ct, verbose);
						acNum ++;
						ac2 ++;
					}
					else {
						ft->clearF2ChildTmpCs(ct, verbose);
					}
				}
				else if(randType < 8) {

					if(ct->ctType == "nwc" || ct->ctType == "wc" || ct->ctType == "bulge13" || ct->ctType == "AG" || ct->ctType == "GA") continue;
					if(!ft->et->para.ctRandMove) continue;
					mvMut = mvLib.getRandomMove2(ct->cm);
					ft->updateCtChildTmpCs(ct, mvMut, verbose);
					mutE = ft->ctMutEnergy(ct,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0){
						curEne += mutE.first;
						ft->acceptCtChildTmpCs(ct, verbose);
						acNum ++;
						ac3 ++;
					}
					else {
						ft->clearCtChildTmpCs(ct, verbose);
					}
				}
				else if(ct->hasThreeBaseFragment){
					if(!ft->et->para.f3Move) continue;
					f3Frag = ct->f3Lib->getRandomFrag();
					ft->updateF3ChildTmpCs(ct, f3Frag, verbose);
					mutE = ft->f3MutEnergy(ct,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(0.3*mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF3ChildTmpCs(ct, verbose);
						acNum ++;
						ac4 ++;
					}
					else {
						ft->clearF3ChildTmpCs(ct, verbose);
					}
				}

			}
			else {
				randType = rand()%10;
				if(randType < 1 && freeNodeNum > 0){
					if(!ft->et->para.singleBaseMove) continue;
					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					mvMut = mvLib.getRandomMove2();
					ft->updateSingleBaseCoordTmp(node, mvMut, false);
					mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0){
						curEne += mutE.first;
						ft->acceptSingleBaseCoordTmp(node, verbose);
						acNum ++;
						ac5 ++;
					}
					else {
						ft->clearSingleBaseCoordTmp(node, verbose);
					}
				}
				else if(randType < 3 && freeNodeNum > 0){
					if(!ft->et->para.reverseRotMove) continue;
					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateReverseBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptReverseBaseRotamerTmp(node, verbose);
						acNum ++;
						ac6 ++;
					}
					else {
						ft->clearReverseBaseRotamerTmp(node, verbose);
					}
				}
				else if(rfn > 0){
					node = ft->riboFlexibleNodes[rand()%rfn];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT, shift,verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptBaseRotamerTmp(node, verbose);
						acNum ++;
						ac7 ++;
					}
					else {
						ft->clearBaseRotamerTmp(node, verbose);
					}
				}
			}
		}
		BRTreeInfo* x = ft->getTreeInfo();
		printf("T=%9.5f Ene: %7.3f rmsd: %6.3f\n",T,curEne, x->rmsd(init));
		delete x;
	}

	connectWT = 1.0;
	clashWT = 1.0;
	breakCTWT = 1.0;

	for(double T=0.1;T>0.01;T=T*0.9){
		int acNum = 0;
		int ac1 = 0;
		int ac2 = 0;
		int ac3 = 0;
		int ac4 = 0;
		int ac5 = 0;
		int ac6 = 0;
		int ac7 = 0;
		int indexX;
		double phoMutE;
		for(int k=0;k<stepNum;k++){
			if(k % outFreq == 0 && printTraj) {
				BRTreeInfo* x = ft->getTreeInfo();
				x->printPDB(of, modelID);
				modelID++;
				delete x;
			}
			randIndex = rand()%fcn;
			if(randIndex < fc){
				ct = ft->flexibleConnectionList[ctChooseIndex[rand()%ctChooseNum]];
				mvMut = mvLib.getRandomMove1(ct->cm);
				ft->updateCtChildTmpCs(ct, mvMut, verbose);
				mutE = ft->ctMutEnergy(ct,breakCTWT, connectWT,clashWT, shift, verbose);
				if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
					curEne += mutE.first;
					ft->acceptCtChildTmpCs(ct, verbose);
					acNum ++;
					ac3 ++;
				}
				else {
					ft->clearCtChildTmpCs(ct, verbose);
				}
			}
			else {
				randType = rand()%10;
				if(randType < 2 && freeNodeNum > 0){
					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					mvMut = mvLib.getRandomMove2();
					ft->updateSingleBaseCoordTmp(node, mvMut, false);
					mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT, shift,verbose);
					if(mutE.second < 0){
						curEne += mutE.first;
						ft->acceptSingleBaseCoordTmp(node, verbose);
						acNum ++;
						ac5 ++;
					}
					else {
						ft->clearSingleBaseCoordTmp(node, verbose);
					}
				}
				else if(randType < 3 && freeNodeNum > 0){
					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateReverseBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->singleBaseMutEnergy(node,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptReverseBaseRotamerTmp(node, verbose);
						acNum ++;
						ac6 ++;
					}
					else {
						ft->clearReverseBaseRotamerTmp(node, verbose);
					}
				}
				else {
					node = ft->riboFlexibleNodes[rand()%rfn];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT, shift,verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptBaseRotamerTmp(node, verbose);
						acNum ++;
						ac7++;
					}
					else {
						ft->clearBaseRotamerTmp(node, verbose);
					}
				}
			}
		}
		BRTreeInfo* x = ft->getTreeInfo();
		printf("T=%9.5f Ene: %7.3f rmsd: %6.3f\n",T,curEne, x->rmsd(init));
//		printf("T=%9.5f WT1=%6.5f WT2=%6.5f WT3=%6.5f shift=%5.3f ac1=%4d ac2=%4d ac3=%4d ac4=%4d ac5=%4d ac6=%4d ac7=%4d addEne: %7.3f totEne: %7.3f rmsd: %6.3f\n",T,breakCTWT, connectWT, clashWT,shift, ac1,ac2,ac3,ac4,ac5,ac6, ac7,curEne,ft->totalEnergy(1.0, 1.0,1.0, shift, false), x->rmsd(init));
		delete x;
	}

	clock_t end = clock();
	cout << "time: " << (float)(end-start)/CLOCKS_PER_SEC << "s" << endl;

	if(!printTraj){
		BRTreeInfo* x = ft->getTreeInfo();
		x->printPDB(of, modelID);
		of << "time: " << (float)(end-start)/CLOCKS_PER_SEC << " s" << endl;
		of << "ene: " << x->ene << endl;
		delete x;
	}
	of.close();
}



void MCRun::optimize(double t0, double kStep){

	int stepNum = 0;

	for(int i=0;i<ft->flexibleConnectionList.size();i++){
		BRConnection* ct = ft->flexibleConnectionList[i];
		if(ct->ctType == "wc" || ct->ctType == "wcNb")
		{
			stepNum += 400;
		}
		else if(ct->ctType == "nwc" || ct->ctType == "nwcNb") {
			stepNum += 2000;
		}
		else {
			stepNum += 4000;
		}
	}

	stepNum = (int)(stepNum*kStep);

	bool verbose = false;


	double connectWT = 1.0;
	double clashWT = 1.0;
	double breakCTWT = 1.0;

	double shift = 0.0;

	double curEne = ft->totalEnergy(1.0, 1.0, 1.0, shift,verbose);
	double lastEne = curEne;
	pair<double,double> mutE;


	vector<BRConnection*> flexConnectList = ft->flexibleConnectionList;
	vector<BRNode*> flexNodes = ft->flexibleNodes;
	int fc = flexConnectList.size();
	int fn = flexNodes.size();
	int tn = ft->seqLen;
	int freeNodeNum = ft->freeNodeIDs.size();

	for(int i=0;i<fc;i++){
		printf("flexible connect: %d-%d\n", flexConnectList[i]->fatherNode->seqID, flexConnectList[i]->childNode->seqID);
	}

	for(int i=0;i<fn;i++){
		printf("flexible node: %d\n", flexNodes[i]->seqID);
	}

	for(int i=0;i<freeNodeNum;i++){
		printf("free node: %d\n",  ft->freeNodeIDs[i] );
	}

	int rfn = ft->riboFlexibleNodes.size();
	int fcn = fc + fn;
	int randIndex;
	BRConnection* ct;
	BRNode* node;

	BaseMoveLibrary mvLib;

	CsMove mvMut;
	int randType;

	BaseRotamer* rotMut;
	F2Fragment* frag;
	F3Fragment* f3Frag;

	int outFreq = ft->et->para.outFreq;
	int modelID = 1;


	for(double T=t0;T>0.1;T=T*0.9){
		int acNum = 0;
		int ac1 = 0;
		int ac2 = 0;
		int ac3 = 0;
		int ac4 = 0;
		int ac5 = 0;
		int ac6 = 0;
		int indexX;
		double phoMutE;
		for(int k=0;k<stepNum;k++){

			randIndex = rand()%fcn;
			if(randIndex < fc){
				ct = flexConnectList[randIndex];
				randType = rand()%10;
				if(randType < 3) {
					if(ct->f2Lib == NULL) continue;
					frag = ct->f2Lib->getRandomFrag();
					ft->updateF2ChildTmpCs(ct, frag, false);
					mutE = ft->f2MutEnergy( ct, breakCTWT, connectWT, clashWT,shift,false);

					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF2ChildTmpCs(ct, verbose);
						acNum ++;
						ac1 ++;
					}
					else {
						ft->clearF2ChildTmpCs(ct, verbose);
					}
				}
				else if(randType < 6) {
					if(ct->f2Lib == NULL) continue;
					frag = ct->f2Lib->getRandomFragLevel2(ct->f2Frag->level1ID);
					ft->updateF2ChildTmpCs(ct, frag, false);
					mutE = ft->f2MutEnergy(ct,breakCTWT,  connectWT,clashWT,shift,false);

					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF2ChildTmpCs(ct, verbose);
						acNum ++;
						ac2 ++;
					}
					else {
						ft->clearF2ChildTmpCs(ct, verbose);
					}
				}
				else if(randType < 7) {
					mvMut = mvLib.getRandomMove2(ct->cm);
					ft->updateCtChildTmpCs(ct, mvMut, verbose);
					mutE = ft->ctMutEnergy(ct,breakCTWT, connectWT,clashWT, shift,verbose);

					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptCtChildTmpCs(ct, verbose);
						acNum ++;
						ac3 ++;
					}
					else {
						ft->clearCtChildTmpCs(ct, verbose);
					}
				}

				else if(ct->hasThreeBaseFragment){
					f3Frag = ct->f3Lib->getRandomFrag();
					ft->updateF3ChildTmpCs(ct, f3Frag, verbose);
					mutE = ft->f3MutEnergy(ct,breakCTWT, connectWT,clashWT, shift,verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptF3ChildTmpCs(ct, verbose);
						acNum ++;
						ac4 ++;
					}
					else {
						ft->clearF3ChildTmpCs(ct, verbose);
					}
				}

			}
			else {
				randType = rand()%10;
				if(randType < 3){
					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					mvMut = mvLib.getRandomMove2();
					ft->updateSingleBaseCoordTmp(node, mvMut, false);
					mutE = ft->singleBaseMutEnergy(node,breakCTWT, connectWT,clashWT, shift,verbose);

					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptSingleBaseCoordTmp(node, verbose);
						acNum ++;
						ac5 ++;
					}
					else {
						ft->clearSingleBaseCoordTmp(node, verbose);
					}
				}
				else {
					node = ft->riboFlexibleNodes[rand()%rfn];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptBaseRotamerTmp(node, verbose);
						acNum ++;
					}
					else {
						ft->clearBaseRotamerTmp(node, verbose);
					}
				}
			}
		}

		BRTreeInfo* x = ft->getTreeInfo();
		printf("T=%9.5f Ene: %7.3f rmsd: %6.3f\n",T,curEne, x->rmsd(init));
		//printf("T=%9.5f WT=%6.5f ac1=%5d ac2=%5d ac3=%5d ac4=%5d ac5=%5d ac6=%5d addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,connectWT, ac1,ac2,ac3,ac4,ac5,ac6,curEne,ft->totalEnergy(1.0, 1.0,1.0,shift,false), x->rmsd(init));
		delete x;
	}

	for(double T=0.1;T>0.01;T=T*0.7){
		int acNum = 0;
		int ac1 = 0;
		int ac2 = 0;
		int ac3 = 0;
		int ac4 = 0;
		int ac5 = 0;
		int ac6 = 0;
		int indexX;
		double phoMutE;
		for(int k=0;k<stepNum;k++){
			randIndex = rand()%fcn;
			if(randIndex < fc){
				ct = flexConnectList[randIndex];
				{
					mvMut = mvLib.getRandomMove1(ct->cm);
					ft->updateCtChildTmpCs(ct, mvMut, verbose);
					mutE = ft->ctMutEnergy(ct,breakCTWT, connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptCtChildTmpCs(ct, verbose);
						acNum ++;
						ac3 ++;
					}
					else {
						ft->clearCtChildTmpCs(ct, verbose);
					}
				}
			}
			else {
				randType = rand()%10;
				if(randType < 3){

					node = ft->nodes[ft->freeNodeIDs[rand()%freeNodeNum]];
					mvMut = mvLib.getRandomMove2();
					ft->updateSingleBaseCoordTmp(node, mvMut, false);

					mutE = ft->singleBaseMutEnergy(node,breakCTWT, connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptSingleBaseCoordTmp(node, verbose);
						acNum ++;
						ac5 ++;
					}
					else {
						ft->clearSingleBaseCoordTmp(node, verbose);
					}
				}
				else {
					node = ft->riboFlexibleNodes[rand()%rfn];
					rotMut = ft->rotLib->getRandomRotamerLv1(node->baseType);
					ft->updateBaseRotamerTmp(node, rotMut, verbose);
					mutE = ft->baseRotamerMutEnergy(node,breakCTWT,connectWT,clashWT,shift, verbose);
					if(mutE.second < 0 || rand()*exp(mutE.second/T) < RAND_MAX){
						curEne += mutE.first;
						ft->acceptBaseRotamerTmp(node, verbose);
						acNum ++;
					}
					else {
						ft->clearBaseRotamerTmp(node, verbose);
					}
				}
			}
		}
		BRTreeInfo* x = ft->getTreeInfo();
		printf("T=%9.5f Ene: %7.3f rmsd: %6.3f\n",T,curEne, x->rmsd(init));
		//printf("T=%9.5f WT=%6.5f ac1=%5d ac2=%5d ac3=%5d ac4=%5d ac5=%5d ac6=%5d addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,connectWT, ac1,ac2,ac3,ac4,ac5,ac6,curEne,ft->totalEnergy(1.0, 1.0,1.0,shift,false), x->rmsd(init));
		delete x;
	}

}


/*
void MCRun::optimizeBackbone(const string& output){
	cout << "opt backbone:" << endl;
	double T0= 100.0;
	double T1 = 0.1;
	int stepNum = 20000;

	bool verbose = false;

	double curEne = ft->totalEnergy(verbose);
	double lastEne = curEne;
	double mutE;


	srand(time(0));
	vector<BRConnection*> flexConnectList = ft->flexibleConnectionList;
	vector<BRNode*> flexNodes = ft->flexibleNodes;
	int fc = flexConnectList.size();
	int fn = flexNodes.size();
	int fcn = fc + fn;


	int randIndex;
	BRConnection* ct;
	BRNode* node;

	int randType;
	int randMvIndex1, randMvIndex2;
	int mvType;
	BaseRotamer* rotMut;

	for(double T=T0;T>T1;T=T*0.9){
		//printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f deltaE: %9.3f\n",T,0,curEne,ft->totalEnergy(false), ft->totalEnergy(false)-curEne);

		int acNum = 0;
		for(int k=0;k<stepNum;k++){
			//cout << "step: " << k << " ";
			randIndex = rand()%fn;

			{
			}

			//ft->printDetailEnergy();
			//printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f deltaE: %9.3f mutE: %9.3f index: %d\n",T,0,curEne,ft->totalEnergy(false), ft->totalEnergy(false)-curEne, mutE, randIndex);
		}
		BRTreeInfo* x = ft->getTreeInfo();
		printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,acNum,curEne,ft->totalEnergy(false), x->rmsd(init));
	}

	BRTreeInfo* info = ft->getTreeInfo();
	info->printPDB(output);
	ft->printDetailEnergy();

}
*/

/*
void MCRun::simpleMCBasic(const string& output){

	cout << "run MC basic: " << endl;
	double T0= 100.0;
	double T1 = 0.01;
	int stepNum = 40000;


	bool verbose = false;
	cout << "init energy: " << simpFt->totalEnergy(verbose) << endl;
	//ft->printDetailEnergy();
	simpFt->randInit();
	cout << "rand energy: " << simpFt->totalEnergy(verbose) << endl;
	//ft->printDetailEnergy();

	double curEne = simpFt->totalEnergy(verbose);
	double lastEne = curEne;
	double mutE;


	srand(time(0));
	vector<BRConnectionBasic*> flexConnectList = simpFt->flexibleConnectionList;

	int fc = flexConnectList.size();


	int randIndex;
	BRConnectionBasic* ct;

	int rotNum;
	int randRotIndex;
	CsMove mvMut;

	for(double T=T0;T>T1;T=T*0.9){
		//printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f deltaE: %9.3f\n",T,0,curEne,ft->totalEnergy(false), ft->totalEnergy(false)-curEne);

		int acNum = 0;
		for(int k=0;k<stepNum;k++){
			randIndex = rand()%fc;
			if(randIndex < fc){
				ct = flexConnectList[randIndex];
				rotNum = ct->mvLib->rotNum;
				randRotIndex = rand()%rotNum;
				if(rand()%3 == 0)
					mvMut = ct->mvLib->mvLibLevel1[randRotIndex];
				else
					mvMut = ct->mvLib->mvLibLevel2[ct->mvRotType* rotNum + randRotIndex];


				simpFt->updateConnectionChildTmpCs(ct, mvMut, verbose);
				mutE = simpFt->connectionMutEnergy(ct, verbose);
				if(mutE < 0 || rand()*exp(mutE/T) < RAND_MAX){
					curEne += mutE;
					ct->mvRotType = randRotIndex;
					simpFt->acceptConnectionChildTmpCs(ct, mvMut, verbose);
					acNum ++;
				}
				else {
					simpFt->clearConnectionChildTmpCs(ct, verbose);
				}
			}


			//ft->printDetailEnergy();
			//printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f deltaE: %9.3f mutE: %9.3f index: %d\n",T,0,curEne,ft->totalEnergy(false), ft->totalEnergy(false)-curEne, mutE, randIndex);

		}
		BRTreeInfoBasic* b = simpFt->getTreeInfo();
		printf("T=%9.5f accept=%5d addEne: %9.3f totEne: %9.3f rmsd: %6.3f\n",T,acNum,curEne,simpFt->totalEnergy(false), b->rmsd(initBasic));
	}

	BRTreeInfoBasic* info = simpFt->getTreeInfo();
	info->printPDB(output);
	simpFt->printDetailEnergy();

}
*/

void MCRun::debug(){

	vector<BRConnection*> flexConnectList = ft->flexibleConnectionList;
	vector<BRNode*> flexNodes = ft->flexibleNodes;

	cout << "print connections: " << endl;

	ft->printConnections();
	//ft->printNodesPartition();

	double connectWT = 0.1;
	double clashWT = 0.5;
	double breakCTWT = 0.3;
	clock_t start = clock();

	double shift = 0.0;

	double initEnergy = ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false);
	double totE = initEnergy;

	cout << "check total energy: " << endl;
	ft->checkTotalEnergy(shift);
	cout << "check tmp energy: " << endl;
	ft->checkTmpTotalEnergy(shift);
	cout << "print detail energy: " << endl;
	ft->printDetailEnergy();

	cout << "check connection: " << endl;
	ft->checkConnection();
	ft->checkNode();
	ft->checkRibose();
	ft->checkPho();

	BaseMoveLibrary* moveLib = new BaseMoveLibrary();
	// debug single base move
	for(int i=0;i<10000;i++) {

		{
		cout << "step: " << i << endl;
		cout << "single base move" << endl;

		BRNode* node = flexNodes[rand()%flexNodes.size()];
		cout << "select node: " << node->seqID << endl;
		CsMove cm = moveLib->getRandomMove1();

		cout << "update node child coordinate" << endl;
		ft->updateSingleBaseCoordTmp(node, cm, true);

		cout << "track coordinate change" << endl;
		ft->trackCoordinateChangeSingleBase(node);

		cout << "track finished" << endl;

		pair<double,double> mutE = ft->singleBaseMutEnergy(node, breakCTWT, connectWT,clashWT,shift,true);

		cout << "mutE: ";
		cout << mutE.first << endl;
		cout << mutE.second << endl;

		if(i%2==0){
			cout << "accept connection node: "  << endl;
			ft->acceptSingleBaseCoordTmp(node, true);
			totE += mutE.second;
		}
		else {
			cout << "reject connection node: " << endl;
			ft->clearSingleBaseCoordTmp(node, false);
		}

		cout << "check energy: " << endl;
		ft->checkTotalEnergy(shift);
		cout << "check tmp energy: " << endl;
		ft->checkTmpTotalEnergy(shift);

		cout << "check connection: " << endl;
		ft->checkConnection();
		ft->checkNode();
		ft->checkRibose();
		ft->checkPho();

		ft->printDetailEnergy();

		printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));
		}

		{

			BRConnection* ct = flexConnectList[rand()%flexConnectList.size()];

			cout << "ct move" << endl;

			cout << "select connection: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << endl;

			CsMove cm = moveLib->getRandomMove1(ct->cm);

			cout << "update connection child coordinate" << endl;
			ft->updateCtChildTmpCs(ct, cm, true);

			cout << "track coordinate change" << endl;
			ft->trackCoordinateChangeCt(ct);

			cout << "track finished" << endl;
			pair<double,double> mutE = ft->ctMutEnergy(ct,breakCTWT, connectWT, clashWT,shift,false);

			cout << "mutE: ";
			cout << mutE.first << endl;
			cout << mutE.second << endl;

			if(i%2==0){
				cout << "accept connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
				ft->acceptCtChildTmpCs(ct, false);
				totE += mutE.second;
			}
			else {
				cout << "reject connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
				ft->clearCtChildTmpCs(ct, false);
			}

			cout << "check energy: " << endl;
			ft->checkTotalEnergy(shift);
			cout << "check tmp energy: " << endl;
			ft->checkTmpTotalEnergy(shift);

			cout << "check connection: " << endl;
			ft->checkConnection();
			ft->checkNode();
			ft->checkRibose();
			ft->checkPho();
			ft->printDetailEnergy();
			printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));

		}

	}

	// debug connection move
	for(int i=0;i<10000;i++){

		cout << endl << "step: " << i << endl;
		cout << "conection move" << endl;

		BRConnection* ct = flexConnectList[rand()%flexConnectList.size()];
		cout << "select connection: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << endl;

		CsMove cm = moveLib->getRandomMove1(ct->cm);

		cout << "update connection child coordinate" << endl;
		ft->updateCtChildTmpCs(ct, cm, true);

		cout << "track coordinate change" << endl;
		ft->trackCoordinateChangeCt(ct);

		cout << "track finished" << endl;
		pair<double,double> mutE = ft->ctMutEnergy(ct, breakCTWT, connectWT,clashWT,shift,false);

		cout << "mutE: ";
		cout << mutE.first << endl;
		cout << mutE.second << endl;

		if(i%2==0){
			cout << "accept connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
			ft->acceptCtChildTmpCs(ct, false);
			totE += mutE.second;
		}
		else {
			cout << "reject connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
			ft->clearCtChildTmpCs(ct, false);
		}

		cout << "check energy: " << endl;
		ft->checkTotalEnergy(shift);
		cout << "check tmp energy: " << endl;
		ft->checkTmpTotalEnergy(shift);

		cout << "check connection: " << endl;
		ft->checkConnection();
		ft->checkNode();
		ft->checkRibose();
		ft->checkPho();

		cout << "print detail energy: " << endl;
		ft->printDetailEnergy();

		printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));

	}

	// debug f2 fragment change
	for(int i=0;i<10000;i++){

		cout << endl << "step: " << i << endl;
		cout << "f2 fragment move" << endl;

		BRConnection* ct = flexConnectList[rand()%flexConnectList.size()];
		cout << "select connection: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID << endl;

		F2Fragment* frag = ct->f2Lib->getRandomFrag();

		cout << "update connection child coordinate" << endl;
		ft->updateF2ChildTmpCs(ct, frag, true);

		cout << "track coordinate change" << endl;
		ft->trackCoordinateChangeF2(ct);

		cout << "mutE: ";

		pair<double,double> mutE = ft->f2MutEnergy(ct, breakCTWT, connectWT,clashWT,shift,false);
		cout << mutE.first << endl;

		CsMove cm1 = ft->nodes[4]->cs1 - ft->nodes[3]->cs1;
		CsMove cm2 = ft->nodes[4]->tmpCs1 - ft->nodes[3]->tmpCs1;

		if(i%2==0){
			cout << "accept connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
			ft->acceptF2ChildTmpCs(ct, false);
			totE += mutE.second;
		}
		else {
			cout << "reject connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << endl;
			ft->clearF2ChildTmpCs(ct, false);
		}

		cout << "check energy: " << endl;
		ft->checkTotalEnergy(shift);
		cout << "check tmp energy: " << endl;
		ft->checkTmpTotalEnergy(shift);

		cout << "check connection: " << endl;
		ft->checkConnection();
		ft->checkNode();
		ft->checkRibose();
		ft->checkPho();

		cout << "print detail energy: " << endl;
		ft->printDetailEnergy();

		printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));

	}


	//debug f3 fragment change
	for(int i=0;i<10000;i++){

		cout << endl << "step: " << i << endl;
		cout << "f3 fragment move" << endl;
		BRConnection* ct = flexConnectList[rand()%flexConnectList.size()];
		if(!ct->hasThreeBaseFragment) {
			i--;
			continue;
		}

		cout << "select connection: " << ct->fatherNode->seqID << "->" << ct->childNode->seqID  << "->" << ct->childNode->midChild->seqID << endl;

		F3Fragment* frag = ct->f3Lib->getRandomFrag();

		cout << "update child coords" << endl;
		ft->updateF3ChildTmpCs(ct, frag, true);

		cout << "track coord change" << endl;
		ft->trackCoordinateChangeF3(ct);


		cout << "mutE: ";
		pair<double,double> mutE = ft->f3MutEnergy(ct,breakCTWT, connectWT,clashWT, shift,true);
		cout << mutE.first << endl;

		if(i%2==0){
			cout << "accept connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << "->" << ct->childNode->midChild->seqID << endl;
			ft->acceptF3ChildTmpCs(ct, false);
			totE += mutE.second;
		}
		else {
			cout << "reject connection cm: " << ct->fatherNode->seqID << " " << ct->childNode->seqID << "->" << ct->childNode->midChild->seqID << endl;
			ft->clearF3ChildTmpCs(ct, false);
		}

		cout << "check energy: " << endl;
		ft->checkTotalEnergy(shift);
		cout << "check tmp energy: " << endl;
		ft->checkTmpTotalEnergy(shift);

		cout << "check connection: " << endl;
		ft->checkConnection();
		ft->checkNode();
		ft->checkRibose();
		ft->checkPho();

		//cout << "print detail energy: " << endl;
		ft->printDetailEnergy();

		printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));

	}


	// debug base rotamer change
	for(int i=0;i<10000;i++) {
		cout << "step: " << i << endl;
		cout << "rotamer move" << endl;
		BRNode* node = flexNodes[rand()%flexNodes.size()];
		cout << "select node: " << node->seqID << endl;
		CsMove cm = moveLib->getRandomMove1();

		cout << "update node child coordinate" << endl;
		BaseRotamer* rot = ft->rotLib->getRandomRotamerLv1(node->baseType);
		ft->updateBaseRotamerTmp(node, rot, true);


		cout << "track coordinate change" << endl;
		ft->trackCoordinateChangeRotamer(node);


		pair<double,double> mutE = ft->baseRotamerMutEnergy(node,breakCTWT, connectWT,clashWT, shift,true);

		cout << "mutE: ";
		cout << mutE.first << endl;

		CsMove cm1 = ft->nodes[4]->cs1 - ft->nodes[3]->cs1;
		CsMove cm2 = ft->nodes[4]->tmpCs1 - ft->nodes[3]->tmpCs1;

		if(i%2==0){
			cout << "accept connection node: "  << endl;
			ft->acceptBaseRotamerTmp(node, false);
			totE += mutE.second;
		}
		else {
			cout << "reject connection node: " << endl;
			ft->clearBaseRotamerTmp(node, false);
		}

		cout << "check energy: " << endl;
		ft->checkTotalEnergy(shift);
		cout << "check tmp energy: " << endl;
		ft->checkTmpTotalEnergy(shift);

		cout << "check connection: " << endl;
		ft->checkConnection();
		ft->checkNode();
		ft->checkRibose();
		ft->checkPho();

		//cout << "print detail energy: " << endl;
		ft->printDetailEnergy();
		printf("initE: %10.3f addE: %10.3f finalE: %10.3f\n", initEnergy, totE, ft->totalEnergy(breakCTWT, connectWT,clashWT,shift,false));
	}

	delete moveLib;
}


MCRun::~MCRun() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpred */
