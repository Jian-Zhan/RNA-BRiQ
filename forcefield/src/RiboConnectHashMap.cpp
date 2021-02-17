/*
 * RiboConnectHashMap.cpp
 *
 *  Created on: Jul 23, 2019
 *      Author: s2982206
 */

#include <forcefield/RiboConnectHashMap.h>

namespace NSPforcefield {

RiboConnectHashMap::RiboConnectHashMap() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	string fileName = path+"keyEnergy/riboConnect/key.po4.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open ribo-connect key file " << fileName << endl;
		exit(1);
	}
	int index = 0;
	vector<string> spt;
	while(getline(file,s)){
		NSPtools::splitString(s," ",&spt);
		string key = spt[0];
		keyToIndex[key] = index;
		PList.push_back(XYZ(atof(spt[1].c_str()), atof(spt[2].c_str()), atof(spt[3].c_str())));
		O5List.push_back(XYZ(atof(spt[4].c_str()), atof(spt[5].c_str()), atof(spt[6].c_str())));
		OP1List.push_back(XYZ(atof(spt[7].c_str()), atof(spt[8].c_str()), atof(spt[9].c_str())));
		OP2List.push_back(XYZ(atof(spt[10].c_str()), atof(spt[11].c_str()), atof(spt[12].c_str())));
		eneList.push_back(atof(spt[13].c_str()));
		index++;
	}

	file.close();

	fileName = path+"keyEnergy/riboConnect/rep.keys";
	file.open(fileName.c_str(), ios::in);
	while(getline(file, s)){
		string key = s.substr(0,9);
		repKeys.push_back(key);

		index = keyToIndex[key];
		repEnes.push_back(eneList[index]);
		repPList.push_back(PList[index]);
		repO5List.push_back(O5List[index]);
		repOP1List.push_back(OP1List[index]);
		repOP2List.push_back(OP2List[index]);
	}
	file.close();

	repNum = repKeys.size();
}

void RiboConnectHashMap::setWeight(double wt){
	for(int i=0;i<eneList.size();i++){
		eneList[i] = eneList[i]*wt;
	}
	for(int i=0;i<repEnes.size();i++){
		repEnes[i] = repEnes[i]*wt;
	}
}

RiboConnectHashMap::~RiboConnectHashMap() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
