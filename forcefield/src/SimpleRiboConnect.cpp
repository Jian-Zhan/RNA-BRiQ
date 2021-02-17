/*
 * SimpleRiboConnect.cpp
 *
 *  Created on: Aug 7, 2019
 *      Author: s2982206
 */

#include <forcefield/SimpleRiboConnect.h>

namespace NSPforcefield {

SimpleRiboConnect::SimpleRiboConnect() {

	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	string fileName;
	this->wt = 1.0;

	fileName = path + "simpleRiboConnect.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	double ene;
	while(getline(file,s)){
		ene = atof(s.c_str());
		this->eneList.push_back(ene);
	}
	file.close();

	fileName = path + "phoLib/pho1.rotLib1";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	getline(file, s);
	double x,y,z;
	for(int i=0;i<10;i++){
		file >> x >> y >> z;
		XYZ t(x,y,z);
		phoLib1.push_back(t);
		phoLib1Ene.push_back(et1.getEnergy(t));
	}
	file.close();

	for(int id=0;id<10;id++){
		char ss[10];
		sprintf(ss, "%d", id);
		fileName = path + "phoLib/pho1.rotLib2-" + string(ss);
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()){
			cout << "fail to open file: " << fileName << endl;
			exit(1);
		}
		getline(file, s);
		vector<XYZ> tList;
		vector<double> eList;
		for(int i=0;i<10;i++){
			file >> x >> y >> z;
			XYZ t(x,y,z);
			tList.push_back(t);
			eList.push_back(et1.getEnergy(t));
		}
		phoLib2.push_back(tList);
		phoLib2Ene.push_back(eList);
		file.close();
	}
	string augc = "AUGC";
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			string ss = augc.substr(i,1)+augc.substr(j,1);
			ifstream file;
			string pFile = path+"/split6D/n1Keys/"+ss+".P.ene";
			file.open(pFile);
			if(file.fail()){
				cout << "can't open file: " << pFile << endl;
				exit(0);
			}

			int bpIndex;
			double x,y,z, ene;
			while(file >> bpIndex >> x >> y >> z >> ene){
				XYZ p = XYZ(x,y,z);
				this->nbKeyToP[i*4+j][bpIndex] = PhoBasicLocal(p, ene);
			}
			file.close();
		}
	}


}

SimpleRiboConnect::~SimpleRiboConnect() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
