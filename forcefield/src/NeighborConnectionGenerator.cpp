/*
 * NeighborConnectionGenerator.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: s2982206
 */

#include "forcefield/NeighborConnectionGenerator.h"

namespace NSPforcefield {

NeighborMove::NeighborMove(const string& line){
	vector<string> spt;
	NSPtools::splitString(line, " ", &spt);

	this->move = CsMove(line);
	double x1 = atof(spt[12].c_str());
	double y1 = atof(spt[13].c_str());
	double z1 = atof(spt[14].c_str());
	double x2 = atof(spt[15].c_str());
	double y2 = atof(spt[16].c_str());
	double z2 = atof(spt[17].c_str());
	double ene = atof(spt[18].c_str());
	this->localP1 = XYZ(x1,y1,z1);
	this->localP2 = XYZ(x2,y2,z2);
	this->energy = ene;
}

NeighborMove::~NeighborMove(){

}

NeighborConnectionGenerator::NeighborConnectionGenerator(int type) {
	// TODO Auto-generated constructor stub

	srand(time(NULL));

	string path = NSPdataio::datapath();

	if(type == 1 || type == 2 || type == 3) {
		for(int i=0;i<16;i++) {
			vector<NeighborMove*> list;
			char s[5];
			char t[5];
			sprintf(s, "%d", i);
			sprintf(t, "%d", type);
			string fileName = path + "/rnaDM/nbMove" + t + "/nb-move-" + string(s);
			ifstream f;
			f.open(fileName, ios::in);
			if(f.fail()){
				cout << "can't open file: " << fileName << endl;
				exit(0);
			}
			string line;
			while(getline(f,line)){
				list.push_back(new NeighborMove(line));
			}
			f.close();
			moveListList.push_back(list);
		}
	}

	string s = "1.5580  3.7013  -3.9455 0.8518  0.5154  -0.0943 -0.5010 0.8538  0.1412  0.1533  -0.0730 0.9855   -0.673  -0.901  -4.538   0.315  -5.037  -1.024 -6.500";
	this->move0 = new NeighborMove(s);
}

NeighborMove* NeighborConnectionGenerator::getRandomMove(int pairType){
	int n = this->moveListList[pairType].size();
	int randIndex = rand()%n;
	return this->moveListList[pairType][randIndex];
}

NeighborConnectionGenerator::~NeighborConnectionGenerator() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<16;i++) {
		for(int j=0;j<this->moveListList[i].size();j++) {
			delete this->moveListList[16][j];
		}
	}
}

} /* namespace NSPforceField */
