/*
 * MoveMutator.cpp
 *
 *  Created on: Feb 7, 2019
 *      Author: s2982206
 */

#include <pred/MoveMutator.h>

namespace NSPpred {

MoveMutator::MoveMutator(const string& type, int typeA, int typeB, double theta){
	int pairType = typeA*4+typeB;
	srand(time(NULL));
	string path = NSPdataio::datapath();
	char s[3];
	char sA[3];
	char sB[3];

	sprintf(sA, "%d", typeA);
	sprintf(sB, "%d", typeB);

	sprintf(s, "%d", pairType);
	string fileName;
	if(type == "helixNb")
		fileName = path + "/rnaDM/helixNbMove/helixNb-0.5.cs-" + string(s);
	else if(type == "loopNb" || type == "jump")
		fileName = path + "/rnaDM/loopNbMove/loopNb-0.3.cs-" + string(s);
	else if(type == "wc")
		fileName = path + "/rnaDM/wcMove/wc-0.5.cs-" + string(s);
	else if(type == "nwc")
		fileName = path + "/rnaDM/nwcMove/nwc-0.5.cs-" + string(s);
	else if(type == "nwcNb")
		fileName = path + "/rnaDM/nwcNbMove/nwcNb-0.5.cs-" + string(s);
	else if(type == "base2Ribo")
		fileName = path + "/rnaDM/baseRiboMove/B2R.cs-" + string(sA);
	else if(type == "ribo2Base")
		fileName = path + "/rnaDM/riboBaseMove/R2B.cs-" + string(sB);
	else if(type == "ribo2Ribo")
		fileName = path + "/rnaDM/riboRiboMove/R2R-0.3.cs-" + string(s);

	ifstream f;
	f.open(fileName, ios::in);
	if(f.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}
	string line;
	while(getline(f,line)){
		this->moveList.push_back(new CsMove(line));
	}
	f.close();

	for(int i=0;i<100000;i++) {
		TransMatrix tm;
		tm = tm.randomRot(theta);
		double rx = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		double ry = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		double rz = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		XYZ ori(rx,ry,rz);
		this->disturbanceList.push_back(new CsMove(ori,tm));
	}

	this->moveNum = moveList.size();
	this->dbNum = disturbanceList.size();
	this->proportion = 1.0;
}

MoveMutator::MoveMutator(const string& type, CsMove* wildTypeMove, int typeA, int typeB, double  cutoff){
	int pairType = typeA*4+typeB;
	srand(time(NULL));
	string path = NSPdataio::datapath();
	char s[3];

	sprintf(s, "%d", pairType);
	string fileName;
	if(type == "helixNb")
		fileName = path + "/rnaDM/helixNbMove/helixNb-0.5.cs-" + string(s);
	else if(type == "loopNb" || type == "jump")
		fileName = path + "/rnaDM/loopNbMove/loopNb-0.3.cs-" + string(s);
	else if(type == "wc")
		fileName = path + "/rnaDM/wcMove/wc-0.5.cs-" + string(s);
	else if(type == "nwc")
		fileName = path + "/rnaDM/nwcMove/nwc-0.5.cs-" + string(s);

	ifstream f;
	f.open(fileName, ios::in);
	if(f.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}
	int totNum = 0;
	int accept = 0;
	string line;

	BaseDistanceMatrix wtDM(*wildTypeMove);

	while(getline(f,line)){
		CsMove* mv = new CsMove(line);
		BaseDistanceMatrix dm(*mv);
		totNum++;
		if(wtDM.distanceTo(dm) < cutoff) {
			this->moveList.push_back(mv);
			accept ++;
		}
		else {
			delete mv;
		}
	}
	f.close();

	printf("%s %d %d %5.3f\n",type.c_str(), accept, totNum,1.0*accept/totNum );

	double theta = 2.0;
	for(int i=0;i<100000;i++) {
		TransMatrix tm;
		tm = tm.randomRot(theta);
		double rx = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		double ry = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		double rz = 0.1*theta*(1.0*rand()/RAND_MAX - 0.5);
		XYZ ori(rx,ry,rz);
		this->disturbanceList.push_back(new CsMove(ori,tm));
	}
	this->moveNum = moveList.size();
	this->dbNum = disturbanceList.size();
	this->proportion = 1.0*accept/totNum;
}

CsMove* MoveMutator::randomMove(){
	int randIndex = rand()%this->moveNum;
	return moveList[randIndex];
}

CsMove* MoveMutator::randomDisturbance(){
	int randIndex = rand()%this->dbNum;
	return disturbanceList[randIndex];
}

MoveMutator::~MoveMutator() {
	for(int i=0;i<moveList.size();i++){
		delete moveList[i];
	}
	for(int i=0;i<disturbanceList.size();i++) {
		delete disturbanceList[i];
	}
}


} /* namespace NSPpred */
