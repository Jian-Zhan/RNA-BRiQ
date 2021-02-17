/*
 * RnaAtomLib.cpp
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#include "model/RnaAtomLib.h"

namespace NSPmodel {

RnaAtomLib::RnaAtomLib() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath();
	string libFile = path+"rnaAtoms";
	ifstream file;
	this->rn = RNABaseName();

	file.open(libFile.c_str(), ios::in);
	if(! file.is_open()) {
		cout << "fail to open rna atom library file " << libFile << endl;
		exit(1);
	}

	for(int i=0;i<4;i++) {
		vector<int>* ids = new vector<int>();
		vector<string>* atomNames = new vector<string>();
		vector<string>* scAtomNames = new vector<string>();
		this->baseUniqueIDs.push_back(ids);
		this->baseAtomNames.push_back(atomNames);
		this->baseScAtomNames.push_back(scAtomNames);
	}

	string s,prefix,uniqueName,atomName;
	int uniqueID = 0;
	while(getline(file,s)) {
		prefix = s.substr(0,4);
		if(prefix == "BASE") continue;
		else {
			AtomProperty* ap = new AtomProperty(s);
			int baseType = rn.sinToInt(s[0]);
			this->baseUniqueIDs[baseType]->push_back(uniqueID);
			uniqueName = ap->atomUniqueName;
			atomName = ap->atomName;
			this->uniqueAtomNames.push_back(uniqueName);
			this->properties.push_back(ap);
			this->uniqueNamesToUniqueID[uniqueName] = uniqueID;
			this->baseAtomNames[baseType]->push_back(atomName);
			if(atomName == "P" || atomName == "OP1" || atomName == "OP2" ||
			   (atomName.length() == 3 && atomName[2] == '\''))
				continue;
			this->baseScAtomNames[baseType]->push_back(atomName);
			uniqueID ++;
		}
	}
	file.close();
}

AtomProperty* RnaAtomLib::getAtomProperty(string& uniqueName) const{

	map<string,int>::const_iterator it;
	it = uniqueNamesToUniqueID.find(uniqueName);
	if(it != uniqueNamesToUniqueID.end())
	{
		int id = it->second;
		return properties.at(id);
	}
	else {
		cerr << "undefined atom: " << uniqueName << endl;
		return NULL;
	}
}

vector<string>* RnaAtomLib::getAtomNames(int intType) const{
	if(intType > 3) {
		cerr << "invalid type: " << intType << endl;
		return NULL;
	}
	return this->baseAtomNames[intType];
}

vector<string>* RnaAtomLib::getSidechainAtoms(int intType) const{
	if(intType > 3) {
		cerr << "invalid type: " << intType << endl;
		return NULL;
	}
	return this->baseScAtomNames[intType];
}

string RnaAtomLib::getUniqueName(int id){
	if(id < 0 || id >= uniqueAtomNames.size()){
		cerr << "invalid index: " << id << endl;
		return "X-X";
	}
	return uniqueAtomNames[id];
}

int RnaAtomLib::getUniqueID(string& uniqueName){
	if(uniqueNamesToUniqueID.find(uniqueName) != uniqueNamesToUniqueID.end()){
		return uniqueNamesToUniqueID[uniqueName];
	}
	else {
		cerr << "undefined atom: " << uniqueName << endl;
		return -1;
	}
}


RnaAtomLib::~RnaAtomLib() {
	for(unsigned int i=0;i<this->properties.size();i++){
		delete this->properties[i];
	}

	for(unsigned int i=0;i<this->baseAtomNames.size();i++) {
		delete this->baseAtomNames[i];
		delete this->baseUniqueIDs[i];
		delete this->baseScAtomNames[i];
	}
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
