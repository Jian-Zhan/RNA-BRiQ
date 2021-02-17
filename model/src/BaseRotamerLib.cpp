/*
 * BaseRotamerLib.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: s2982206
 */

#include <model/BaseRotamerLib.h>

namespace NSPmodel {

BaseRotamerLib::BaseRotamerLib() {
	// TODO Auto-generated constructor stub
	string augc = "AUGC";
	string path = NSPdataio::datapath();
	ifstream file;
	string s;

	Parameter para;
	for(int i=0;i<4;i++){
		string baseType = augc.substr(i,1);
		string fileName = path+"baseRotamer/"+baseType+".rot";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()) {
			cout << "fail to open rotamer library file " << fileName << endl;
			exit(1);
		}
		int index = 0;
		while(getline(file,s)){
			BaseRotamer* rot = new BaseRotamer(s);

			rot->resType = i;
			rot->rotTypeLv1 = 0;
			rot->rotType = index;
			rotLib[i][index] = rot;
			rot->energy = rot->energy*para.wtRotamer;
			index++;
		}
		file.close();
	}
}


BaseRotamerLib::~BaseRotamerLib() {
	// TODO Auto-generated destructor stub
	for(int i=0;i<4;i++){
		for(int j=0;j<1500;j++){
			delete rotLib[i][j];
		}
	}
}

} /* namespace NSPforcefield */
