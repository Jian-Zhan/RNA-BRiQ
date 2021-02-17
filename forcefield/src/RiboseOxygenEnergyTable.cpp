/*
 * RiboseOxygenEnergyTable.cpp
 *
 *  Created on: Jul 26, 2019
 *      Author: s2982206
 */

#include <forcefield/RiboseOxygenEnergyTable.h>

namespace NSPforcefield {

RiboseOxygenEnergyTable::RiboseOxygenEnergyTable() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath() + "riboseOxygen/";
	vector<string> oxygenNames;
	oxygenNames.push_back("O2");
	oxygenNames.push_back("O3");
	oxygenNames.push_back("O4");
	oxygenNames.push_back("O5");
	oxygenNames.push_back("OP");
	string augc = "AUGC";

	vector<string> spt;
	double ene;
	int directIndex;
	for(int i=0;i<4;i++){
		for(int j=0;j<5;j++){
			string energyFile = path + augc.substr(i,1) + "-" + oxygenNames[j]+".ene";
			ifstream file;
			file.open(energyFile, ios::in);
			if(file.fail()) {
				cout << "can't open file: " << energyFile << endl;
				exit(0);
			}
			string line;

			while(getline(file, line)){
				NSPtools::splitString(line, " ", &spt);
				ene = atof(spt[0].c_str());
				ene = energyRescale(ene);
				this->etList[i*5+j].push_back(ene);
				if(spt.size() < 2) {
					this->directDep[i*5+j].push_back(false);
					XYZ t;
					this->hbAtom[i*5+j].push_back(t);
					this->ang1List[i*5+j].push_back(0.0);
					this->ang2List[i*5+j].push_back(0.0);
					this->ang3List[i*5+j].push_back(0.0);
				}
				else {
					this->directDep[i*5+j].push_back(true);
					XYZ t(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str()));
					double ang1 = atof(spt[4].c_str());
					double ang2 = atof(spt[5].c_str());
					double ang3 = atof(spt[6].c_str());
					this->hbAtom[i*5+j].push_back(t);
					this->ang1List[i*5+j].push_back(ang1);
					this->ang2List[i*5+j].push_back(ang2);
					this->ang3List[i*5+j].push_back(ang3);
				}
			}
			file.close();
		}
	}

	/*
	 * sep = 1
	 */
	for(int i=0;i<4;i++){
		string energyFile = path + "sep1/" + augc.substr(i,1) + "-OP.ene";
		ifstream file;
		file.open(energyFile, ios::in);
		if(file.fail()) {
			cout << "can't open file: " << energyFile << endl;
			exit(0);
		}
		string line;

		while(getline(file, line)){
			NSPtools::splitString(line, " ", &spt);
			ene = atof(spt[0].c_str());
			ene = energyRescale(ene);
			this->etListP1[i].push_back(ene);
			if(spt.size() < 2) {
				this->directDepP1[i].push_back(false);
				XYZ t;
				this->hbAtomP1[i].push_back(t);
				this->ang1ListP1[i].push_back(0.0);
				this->ang2ListP1[i].push_back(0.0);
				this->ang3ListP1[i].push_back(0.0);
			}
			else {
				this->directDepP1[i].push_back(true);
				XYZ t(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str()));
				double ang1 = atof(spt[4].c_str());
				double ang2 = atof(spt[5].c_str());
				double ang3 = atof(spt[6].c_str());
				this->hbAtomP1[i].push_back(t);
				this->ang1ListP1[i].push_back(ang1);
				this->ang2ListP1[i].push_back(ang2);
				this->ang3ListP1[i].push_back(ang3);
			}
		}
		file.close();
	}

	/*
	 * sep = -1
	 */

	vector<string> oxygenNamesM1;
	oxygenNamesM1.push_back("O4");
	oxygenNamesM1.push_back("OP");
	for(int i=0;i<4;i++){
		for(int j=0;j<2;j++){
			string energyFile = path + "sep-1/" + augc.substr(i,1) + "-" + oxygenNamesM1[j]+".ene";
			ifstream file;
			file.open(energyFile, ios::in);
			if(file.fail()) {
				cout << "can't open file: " << energyFile << endl;
				exit(0);
			}
			string line;

			while(getline(file, line)){
				NSPtools::splitString(line, " ", &spt);
				ene = atof(spt[0].c_str());
				ene = energyRescale(ene);
				this->etListM1[i*2+j].push_back(ene);
				if(spt.size() < 2) {
					this->directDepM1[i*2+j].push_back(false);
					XYZ t;
					this->hbAtomM1[i*2+j].push_back(t);
					this->ang1ListM1[i*2+j].push_back(0.0);
					this->ang2ListM1[i*2+j].push_back(0.0);
					this->ang3ListM1[i*2+j].push_back(0.0);
				}
				else {
					this->directDepM1[i*2+j].push_back(true);
					XYZ t(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str()));
					double ang1 = atof(spt[4].c_str());
					double ang2 = atof(spt[5].c_str());
					double ang3 = atof(spt[6].c_str());
					this->hbAtomM1[i*2+j].push_back(t);
					this->ang1ListM1[i*2+j].push_back(ang1);
					this->ang2ListM1[i*2+j].push_back(ang2);
					this->ang3ListM1[i*2+j].push_back(ang3);
				}
			}
			file.close();
		}
	}

}



RiboseOxygenEnergyTable::~RiboseOxygenEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
