/*
 * AssignRNASS.cpp
 *
 *  Created on: 2020Äê11ÔÂ17ÈÕ
 *      Author: pengx
 */

#include "model/AssignRNASS.h"

namespace NSPmodel {

AssignRNASS::AssignRNASS(RNAPDB* pdb, RnaAtomLib* atLib) {
	// TODO Auto-generated constructor stub

	this->baseList = pdb->getValidBaseList(atLib);

	this->len = this->baseList.size();

	//cout << "length: " << len << endl;
	char cseq[this->len+1];

	cseq[this->len] = '\0';

	for(int i=0;i<len;i++){
		cseq[i] = this->baseList[i]->baseType;
	}

	this->seq = string(cseq);
	//cout << "seq : " << seq << endl;

	if(len == 0){
		this->seq = "";
		this->ssSeq = "";
		this->nwcSeq = "";
		return;
	}

	this->pairIndex = new int[len];
	this->pairDistanceToWC = new double[len];
	this->nwcPairIndex = new int[len];
	this->nwcPairHbondNum = new int[len];
	this->connectWithNeighbor = new bool[len];

	for(int i=0;i<len;i++){
		this->pairIndex[i] = -1;
		this->nwcPairIndex[i] = -1;
		this->nwcPairHbondNum[i] = 0;
		this->pairDistanceToWC[i] = 9.9;
	}
	for(int i=0;i<len;i++){
		this->connectWithNeighbor[i] = false;
	}

	for(int i=0;i<baseList.size();i++){
		RNABase* baseA = baseList[i];
		LocalFrame csA = baseA->getCoordSystem();
		for(int j=i+1;j<baseList.size();j++){
			RNABase* baseB = baseList[j];
			LocalFrame csB = baseB->getCoordSystem();
			if(j==i+1 && baseA->connectToNeighbor(*baseB))
				this->connectWithNeighbor[i] = true;
			BasePair bp(baseA, baseB, atLib);
			if(bp.isWCPair()){
				if(pairIndex[i] == -1 && pairIndex[j] == -1){
					pairIndex[i] = j;
					pairIndex[j] = i;
					pairDistanceToWC[i] = bp.distanceToWCPair();
					pairDistanceToWC[j] = bp.distanceToWCPair();
				}
			}
			else if(pairIndex[i]>=0 && pairIndex[j] == -1){
				double d = bp.distanceToWCPair();
				if(d < pairDistanceToWC[i]) {
					pairIndex[pairIndex[i]] = -1;
					pairDistanceToWC[pairIndex[i]] = 9.9;
					pairIndex[i] = j;
					pairIndex[j] = i;
					pairDistanceToWC[i] = d;
					pairDistanceToWC[j] = pairDistanceToWC[i];
				}
			}
			else if(pairIndex[j]>=0 && pairIndex[i] == -1){
				double d = bp.distanceToWCPair();
				if(d < pairDistanceToWC[j]) {
					pairIndex[pairIndex[j]] = -1;
					pairDistanceToWC[pairIndex[j]] = 9.9;
					pairIndex[i] = j;
					pairIndex[j] = i;
					pairDistanceToWC[i] = d;
					pairDistanceToWC[j] = pairDistanceToWC[i];
				}
			}
		}
	}


	for(int i=0;i<len-1;i++){
		if(!connectWithNeighbor[i]) breakList.push_back(i);
	}


	for(int i=0;i<baseList.size();i++){
		RNABase* baseA = baseList[i];
		LocalFrame csA = baseA->getCoordSystem();
		for(int j=i+1;j<baseList.size();j++){
			RNABase* baseB = baseList[j];
			LocalFrame csB = baseB->getCoordSystem();
			if(j==i+1 && baseA->connectToNeighbor(*baseB))
				this->connectWithNeighbor[i] = true;
			BasePair bp(baseA, baseB, atLib);
			if(bp.isWCPair()) continue;

			if(bp.isHbondedPair()){
				if(pairIndex[i] >=0 || pairIndex[j] >=0) continue;
				if(bp.hbNum <= nwcPairHbondNum[i] || bp.hbNum <= nwcPairHbondNum[j]) continue;
				if(nwcPairIndex[i] >= 0 && nwcPairIndex[j] <0) {
					nwcPairIndex[nwcPairIndex[i]] = -1;
				}
				else if(nwcPairIndex[i] <0 && nwcPairIndex[j] >=0) {
					nwcPairIndex[nwcPairIndex[j]] = -1;
				}
				else if(nwcPairIndex[i] >=0 && nwcPairIndex[j] >=0) {
					nwcPairIndex[nwcPairIndex[i]] = -1;
					nwcPairIndex[nwcPairIndex[j]] = -1;
				}

				nwcPairHbondNum[i] = bp.hbNum;
				nwcPairHbondNum[j] = bp.hbNum;

				nwcPairIndex[i] = j;
				nwcPairIndex[j] = i;
			}

		}
	}


	this->ssSeq = indexToBractString(this->pairIndex);
	this->nwcSeq = indexToBractString(this->nwcPairIndex);
}

string AssignRNASS::indexToBractString(int* index){

	char ss[this->len+1];
	for(int i=0;i<len;i++){
		ss[i] = '.';
	}

	ss[len] = '\0';
	for(int i=0;i<len;i++){
		int pi = index[i];
		if(pi < i) continue;
		if(pi == -1) continue;
		int idA = i;
		int idB = pi;
		int a = 0;
		int b = 0;
		int c = 0;
		int d = 0;
		int e = 0;
		int f = 0;
		int g = 0;
		for(int k=idA+1;k<idB;k++){
			char cc = ss[k];
			if(cc == '(') a++;
			if(cc == ')') a--;
			if(cc == '[') b++;
			if(cc == ']') b--;
			if(cc == '{') c++;
			if(cc == '}') c--;
			if(cc == '<') d++;
			if(cc == '>') d--;
			if(cc == 'A') e++;
			if(cc == 'a') e--;
			if(cc == 'B') f++;
			if(cc == 'b') f--;
			if(cc == 'C') g++;
			if(cc == 'c') g--;
		}

		if(a==0) {
			ss[idA] = '(';
			ss[idB] = ')';
		}
		else if(b == 0){
			ss[idA] = '[';
			ss[idB] = ']';
		}
		else if(c == 0){
			ss[idA] = '{';
			ss[idB] = '}';
		}
		else if(d == 0){
			ss[idA] = '<';
			ss[idB] = '>';
		}
		else if(e == 0){
			ss[idA] = 'A';
			ss[idB] = 'a';
		}
		else if(f == 0){
			ss[idA] = 'B';
			ss[idB] = 'b';
		}
		else if(g == 0){
			ss[idA] = 'C';
			ss[idB] = 'c';
		}
	}
	return string(ss);
}

void AssignRNASS::printInfo(const string& outFile){
	ofstream of;
	of.open(outFile.c_str(), ios::out);
	if(!of.is_open()){
		cout << "fail to open file: " << outFile << endl;
		exit(1);
	}
	of << "seq " << seq << endl;
	of << "sec " << ssSeq << endl;
	of << "nwc " << nwcSeq << endl;
	if(breakList.size() >0){
		of <<"break";
		for(int i=0;i<breakList.size();i++){
			of << " " << breakList[i];
		}
		of << endl;
	}
	of.close();
}

AssignRNASS::~AssignRNASS() {
	// TODO Auto-generated destructor stub
	delete [] connectWithNeighbor;
	delete [] pairDistanceToWC;
	delete [] pairIndex;
	delete [] nwcPairHbondNum;
	delete [] nwcPairIndex;
}

} /* namespace NSPforcefield */
