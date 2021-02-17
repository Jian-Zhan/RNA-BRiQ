/*
 * BasePairEnergyTable.cpp
 *
 *  Created on: Jul 12, 2019
 *      Author: s2982206
 */

#include "forcefield/BasePairEnergyTable.h"

namespace NSPforcefield {

BasePairEnergyTable::BasePairEnergyTable() {
	string path = NSPdataio::datapath();
	string augc = "AUGC";

	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			string ss = augc.substr(i,1)+augc.substr(j,1);
			ifstream file;
			string energyFile1 = path+"/keyEnergy/nbPairs/bw1.2/nbPair-"+ss+".keys.energy";
			string energyFile2 = path+"/keyEnergy/nnbPairs/bw1.2/nnbPair-"+ss+".keys.energy";
			string pFile = path+"/keyEnergy/nbPairs/bw1.2/nbPair-"+ss+".keys.P";

			file.open(energyFile1);
			if(file.fail()) {
				cout << "can't open file: " << energyFile1 << endl;
				exit(0);
			}
			string line;
			while(getline(file, line)){
				nbMapList[i*4+j][line.substr(0,9)] = atof(line.substr(10, line.length()-10).c_str());
			}
			file.close();

			file.open(energyFile2);
			if(file.fail()) {
				cout << "can't open file: " << energyFile2 << endl;
				exit(0);
			}
			while(getline(file, line)){
				nnbMapList[i*4+j][line.substr(0,9)] = atof(line.substr(10, line.length()-10).c_str());
			}
			file.close();


			file.open(pFile);
			if(file.fail()){
				cout << "can't open file: " << pFile << endl;
				exit(0);
			}

			string key;
			double x,y,z;
			while(file >> key >> x >> y >> z){
				this->nbKeyToP[i*4+j][key] = XYZ(x,y,z);
			}
			file.close();
		}
	}

}

void BasePairEnergyTable::setWeight(double wtNb, double wtNnb){
	for(int i=0;i<16;i++){
		for(it=nbMapList[i].begin();it!=nbMapList[i].end();++it){
			nbMapList[i][it->first] = it->second*wtNb;
		}
		for(it=nnbMapList[i].begin();it!=nnbMapList[i].end();++it){
			nnbMapList[i][it->first] = it->second*wtNnb;
		}
	}
}

BasePairEnergyTable::~BasePairEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
