/*
 * BaseStackingEnergyTable.cpp
 *
 *  Created on: Aug 16, 2019
 *      Author: s2982206
 */

#include <forcefield/BaseStackingEnergyTable.h>

namespace NSPforcefield {

BaseStackingEnergyTable::BaseStackingEnergyTable() {
	string path = NSPdataio::datapath() + "stacking/";
	vector<string> spt;
	int idX, idY, idZ;
	for(int i=0;i<24;i++){
		for(int j=0;j<125000;j++)
			this->etList[i].push_back(0);

		char ss[10];
		sprintf(ss, "%d", i);
		string energyFile = path + "stacking.et-" + string(ss);
		ifstream file;
		file.open(energyFile);
		if(file.fail()){
			cout << "can't open file: " << energyFile << endl;
			exit(0);
		}
		string line;
		while(getline(file, line)){
			NSPtools::splitString(line, " ", &spt);
			idX = atoi(spt[0].c_str());
			idY = atoi(spt[1].c_str());
			idZ = atoi(spt[2].c_str());
			this->etList[i][idX*2500+idY*5+idZ] = atof(spt[3].c_str());
		}
		file.close();
	}
	this->wt = 1.0;
}

BaseStackingEnergyTable::~BaseStackingEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
