/*
 * AtomLib.cpp
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#include "model/AtomLib.h"

namespace NSPmodel {

AtomLib::AtomLib() {
	// TODO Auto-generated constructor stub
	string path = NSPdataio::datapath();
	string libFile = path+"resLib.dat";
	ifstream file;
	this-> rn = ResName();

	file.open(libFile.c_str(), ios::in);
	if(! file.is_open())
	{
		 cout << "fail to open file " << libFile << endl;
		 exit (1);
	}


	for(int i=0;i<20;i++)
	{
		vector<int>* ids = new vector<int>();
		vector<string>* atomNames = new vector<string>();
		vector<string>* scAtomNames = new vector<string>();
		this->aaUniqueIDs.push_back(ids);
		this->aaAtomNames.push_back(atomNames);
		this->aaScAtomNames.push_back(scAtomNames);
	}

	string s,prefix,uniqueName,atomName;
	int uniqueID=0;
	while(getline(file,s))
	{
		prefix = s.substr(0,3);
		if(prefix == "RES") continue;
		else
		{
			AtomProperty* ap = new AtomProperty(s);
			int aaType = rn.triToInt(prefix);
			this->aaUniqueIDs.at(aaType)->push_back(uniqueID);

			uniqueName = ap->atomUniqueName;
			atomName = ap->atomName;
			this->uniqueAtomNames.push_back(uniqueName);
			this->properties.push_back(ap);
			this->uniqueNameToUniqueID[uniqueName] = uniqueID;
			this->aaAtomNames.at(aaType)->push_back(atomName);
			if(atomName!= "N" && atomName!="CA" && atomName!="C" && atomName!="O")
				this->aaScAtomNames.at(aaType)->push_back(atomName);
			uniqueID++;
		}
	}

	file.close();
}

bool AtomLib::atomDefined(const string& uniqueName) const
{
	map<string,int>::const_iterator it;
	it = uniqueNameToUniqueID.find(uniqueName);
	if(it != uniqueNameToUniqueID.end())
		return true;
	else
		return false;
}

int AtomLib::uniqueNameToID(const string& uniqueName) const
{
	map<string,int>::const_iterator it;
	it = uniqueNameToUniqueID.find(uniqueName);
	if(it != uniqueNameToUniqueID.end())
		return it->second;
	else
	{
		cerr << "Warning: unknown uniqueName " << uniqueName << endl;
		return -1;
	}
}

string AtomLib::uniqueIDToName(int uniqueID) const
{
	if(uniqueID < 0 || uniqueID >= (int)this->uniqueAtomNames.size())
	{
		cerr << "invalid uniqueID: " << uniqueID << endl;
		exit(0);
	}
	return this->uniqueAtomNames.at(uniqueID);
}

vector<int>* AtomLib::getAminoAcidAtomIDs(const string& triName)
{
	int aaType = rn.triToInt(triName);
	if(aaType < 0 || aaType >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaType << endl;
		return NULL;
	}
	return this->aaUniqueIDs.at(aaType);
}

vector<string>* AtomLib::getAminoAcidAtomNames(const string& triName)
{
	int aaType = rn.triToInt(triName);
	if(aaType < 0 || aaType >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaType << endl;
		return NULL;
	}
	return this->aaAtomNames.at(aaType);
}

vector<string>* AtomLib::getAminoAcidSidechainAtomNames(const string& triName)
{
	int aaType = rn.triToInt(triName);
	if(aaType < 0 || aaType >= (int)this->aaUniqueIDs.size())
	{
		cerr << "invalid atom type: " << aaType << endl;
		return NULL;
	}
	return this->aaScAtomNames.at(aaType);
}

AtomProperty* AtomLib::getAtomProperty(const string& uniqueName)
{
	map<string,int>::const_iterator it;
	it = uniqueNameToUniqueID.find(uniqueName);
	if(it != uniqueNameToUniqueID.end())
	{
		int id = it->second;
		return properties.at(id);
	}

	else
	{
		cerr << "invalid atom name: " << uniqueName << endl;
		return NULL;
	}

}

AtomLib::~AtomLib() {

	for(unsigned int i=0;i<this->properties.size();i++)
	{
		delete this->properties.at(i);
	}
	for(int i=0;i<20;i++)
	{
		delete this->aaAtomNames.at(i);
		delete this->aaUniqueIDs.at(i);
		delete this->aaScAtomNames.at(i);
	}

}

} /* namespace NSPdesignseq */
