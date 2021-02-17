/*
 * ProteinRep.cpp
 *
 *  Created on: 2018/10/30
 *      Author: notxp
 */

#include "model/ProteinRep.h"

namespace NSPmodel {

Atom::Atom(){
	this->name = "X";
	this->type = "X";
	this->resType = "";
	this->alt = ' ';
	this->coord = XYZ();
}

Atom::Atom(const string& line){
	string s = line.substr(12,4);
	this->name = trimString(s);
	if(this->name == "OXT" || this->name == "OXT1")
		this->name = "O";
	this->alt = line.at(16);
	this->type = "";
	float x = atof(line.substr(30,8).c_str());
	float y = atof(line.substr(38,8).c_str());
	float z = atof(line.substr(46,8).c_str());
	this->coord = XYZ(x,y,z);
	guessAtomType();
	s = line.substr(17,3);
	this->resType = trimString(s);
}

Atom::Atom(const string& line, int fileType){
	if(fileType == 0){
		//pdb file
		string s = line.substr(12,4);
		this->name = trimString(s);
		if(this->name == "OXT" || this->name == "OXT1")
			this->name = "O";
		this->alt = line.at(16);
		this->type = "";
		float x = atof(line.substr(30,8).c_str());
		float y = atof(line.substr(38,8).c_str());
		float z = atof(line.substr(46,8).c_str());
		this->coord = XYZ(x,y,z);
		guessAtomType();
		s = line.substr(17,3);
		this->resType = trimString(s);
	}
	else {
		//cif file
		vector<string> spt;
		splitString(line, " ", &spt);
		this->name = spt[3];
		if(this->name[0] == '"')
			this->name = name.substr(1,this->name.length()-2);

		if(this->name == "OXT" || this->name == "OXT1")
			this->name = "O";
		this->alt = spt[9][0];
		if(alt == '?') this->alt = ' ';
		this->type = spt[2];
		float x = atof(spt[10].c_str());
		float y = atof(spt[11].c_str());
		float z = atof(spt[12].c_str());
		this->coord = XYZ(x,y,z);
		guessAtomType();
		this->resType = spt[5];
	}
}

Atom::Atom(string name, const XYZ& coord){
	this->name = name;
	this->coord = coord;
	this->type = "X";
	this->resType = "UNK";
	this->alt = ' ';
}


Atom& Atom::operator=(const Atom& other)
{
	if(this == &other)
		return *this;
	this->name = other.name;
	this->type = other.type;
	this->coord = other.coord;
	this->resType = other.resType;
	return *this;
}

bool Atom::isBackboneAtom() const
{
	return name=="N" || name=="CA" || name=="C" || name=="O" || name=="OXT" || name=="OT1" || name=="OT2";
}

void Atom::guessAtomType()
{
	char c = name.at(0);
	if(name.length() == 1)
		this->type = name;
	else if(name == "FE" || name == "ZN" || name=="MG" || name=="MN" || name=="CL" || name == "SE" || name=="CU" || name=="NA")
		this->type = name;
	else if(c>='A' && c <= 'Z')
	{
	    this->type = string(1,c);

	}
	else
		this->type = string(1,name.at(1));
}

bool Atom::isIon() const
{
	if(type == "FE" || type == "ZN" || type == "MG" || type == "MN" || type == "CA" || type == "CU" || type == "NI")
		return true;
	return false;
}

string& Atom::getName()
{
	return this->name;
}

string& Atom::getType()
{
	return this->type;
}

string& Atom::getResType()
{
	return this->resType;
}



XYZ& Atom::getCoord()
{
	return this->coord;
}

void Atom::setCoord(const XYZ& coord)
{
	this->coord = coord;
}

float Atom::distance(const Atom& other) const
{
	return this->coord.distance(other.coord);
}

void Atom::setAtomType(string type)
{
	this->type = type;
}

void Atom::setResType(string resType)
{
	this->resType = resType;
}

string Atom::nameString() const
{
	char s[10];
	if(this->type.length() > 1)
		sprintf(s,"%-5s",this->name.c_str());
	else
		sprintf(s," %-4s",this->name.c_str());
	string ss = string(s);
	return ss;
}

Atom::~Atom(){

}

Residue::Residue() {
	this->resID = "1";
	this->resSeqID = 0;
	this->atomNum = 0;
	this->triName = "UNK";
	this->chainID = '-';
	this->hasLocalFrame = false;
	this->coordSys;
	this->altLoc = ' ';
	this->hasAltConf = false;

}

Residue::Residue(string resID, char chainID, string triName)
{
	this->resID = resID;
	this->chainID = chainID;
	this->triName = triName;
	this->resSeqID = 0;
	this->atomNum = 0;
	this->hasLocalFrame = false;
	this->altLoc = ' ';
	this->hasAltConf = false;
}

void Residue::addAtom(Atom* a)
{
	this->atomList.push_back(a);
	if(a->isBackboneAtom())
		this->backboneAtoms.push_back(a);
	else
		this->sidechainAtoms.push_back(a);
	this->atomMap[a->getName()] = a;
	this->atomNum ++;
}

void Residue::setResSeqID(int id)
{
    this->resSeqID = id;
}

bool Residue::hasThreeCoreAtoms() const
{
	if(atomMap.find("N") == atomMap.end())
		return false;
	if(atomMap.find("CA") == atomMap.end())
		return false;
	if(atomMap.find("C") == atomMap.end())
		return false;
	return true;
}

void Residue::updateCoordSystem()
{
	if(!this->hasThreeCoreAtoms())
		return;
	Atom* N = this->atomMap.at("N");
	Atom* C = this->atomMap.at("C");
	Atom* CA = this->atomMap.at("CA");

	XYZ n = N->getCoord();
	XYZ c = C->getCoord();
	XYZ ca = CA->getCoord();

	XYZ can = ~(N->getCoord() - CA->getCoord());
	XYZ cac = ~(C->getCoord() - CA->getCoord());

	XYZ z = ~(can^cac);
	XYZ x = ~(can+cac);

	XYZ y = ~(z^x);

	this->coordSys = LocalFrame(ca, x, y, z);

}

bool Residue::hasAtom(const string& atomName) const
{
    if(atomMap.find(atomName) == atomMap.end())
		return false;
	else
        return true;
}

Atom* Residue::getAtom(const string& atomName)
{
    map<string,Atom*>::const_iterator it = atomMap.find(atomName);
    if(it == atomMap.end())
        return NULL;
    else
        return it->second;
}

vector<Atom*>* Residue::getAtomList()
{
    return &this->atomList;
}

vector<Atom*>* Residue::getBackboneAtoms()
{
    return &this->backboneAtoms;
}

vector<Atom*>* Residue::getSidechainAtoms()
{
    return &this->sidechainAtoms;
}

bool Residue::sidechainComplete(AtomLib* atomLib) const
{
	vector<string>* scAtoms = atomLib->getAminoAcidSidechainAtomNames(this->triName);
	if(scAtoms == NULL) return false;
	vector<string>::iterator it;
	for(it=scAtoms->begin();it<scAtoms->end();it++)
	{
		string s = *it;
		if(!hasAtom(s))
			return false;
	}
	return true;
}

XYZ Residue::getCbCoord()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	XYZ localCB = XYZ(-0.942, 0.009, 1.208);
	return this->coordSys.local2globalcrd(localCB);
}

LocalFrame& Residue::getCoordSystem()
{
	if(!hasThreeCoreAtoms())
	{
		cerr << "lack backbone atom information" << endl;
		exit(1);
	}
	if(!this->hasLocalFrame)
		updateCoordSystem();
	return this->coordSys;
}

char Residue::getChainID() const
{
	return this->chainID;
}

int Residue::getResSeqID() const
{
	return this->resSeqID;
}

string Residue::getResID() const
{
	return this->resID;
}

string Residue::getType() const
{
	return this->triName;
}

int Residue::printPDBFormat(ofstream& out, int startAtomID) const
{
    char c = this->resID.at(resID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
    if(c >= '0' && c <= '9')
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    else
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s%3s %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->triName.c_str(),this->chainID,this->resID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    return atomID;
}

Residue::~Residue() {
}

RNABase::RNABase() {
	this->baseID = "1";
	this->chainID = 'A';
	this->baseType = 'N';
	this->baseTypeInt = 4;
	this->atomNum = 0;
	this->baseSeqID = -1;
	this->hasLocalFrame = false;
	this->hasAltConf = false;
	this->altLoc = '-';
}

RNABase::RNABase(const string& baseID, const string& chainID, char baseType){
	this->baseID = baseID;
	this->chainID = chainID;
	this->baseType = baseType;
	if(baseType == 'A' || baseType == 'a')
		this->baseTypeInt = 0;
	else if(baseType == 'U' || baseType == 'u')
		this->baseTypeInt = 1;
	else if(baseType == 'G' || baseType == 'g')
		this->baseTypeInt = 2;
	else if(baseType == 'C' || baseType == 'c')
		this->baseTypeInt = 3;
	else
		this->baseTypeInt = 4;

	this->atomNum = 0;
	this->baseSeqID = -1;
	this->hasLocalFrame = false;
	this->hasAltConf = false;
	this->altLoc = '-';
}

void RNABase::addAtom(Atom* a) {
	this->atomList.push_back(a);
	string atomName = a->name;
	if(atomName == "P" || atomName == "OP1" || atomName == "OP2") {
		this->backboneAtoms.push_back(a);
	}
	else if(atomName.length() == 3 && atomName[2] == '\'')
		this->backboneAtoms.push_back(a);
	else
		this->sidechainAtoms.push_back(a);
	atomMap[a->name] = a;
}

void RNABase::updateCoordSystem() {
	if(hasLocalFrame) return;
	if(this->baseType == 'A' || this->baseType == 'G') {
		if(hasAtom("C1'") && hasAtom("N9") && hasAtom("C4")) {
			XYZ c1 = this->atomMap["C1'"]->coord;
			XYZ n9 = this->atomMap["N9"]->coord;
			XYZ c4 = this->atomMap["C4"]->coord;

			XYZ x = ~(n9-c1);
			XYZ z = ~((n9-c4)^x);
			XYZ y = z^x;
			this->coordSys = LocalFrame(c1, x, y, z);
			this->hasLocalFrame = true;
		}
	}
	else if(this->baseType == 'U' || this->baseType == 'C') {
		if(hasAtom("C1'") && hasAtom("N1") && hasAtom("C2")) {
			XYZ c1 = this->atomMap["C1'"]->coord;
			XYZ n1 = this->atomMap["N1"]->coord;
			XYZ c2 = this->atomMap["C2"]->coord;

			XYZ x = ~(n1-c1);
			XYZ z = ~((n1-c2)^x);
			XYZ y = z^x;
			this->coordSys = LocalFrame(c1, x, y, z);
			this->hasLocalFrame = true;
		}
	}
}

bool RNABase::sidechainComplete(RnaAtomLib* atLib) const{
	vector<string>* atomNames = atLib->getSidechainAtoms(this->baseTypeInt);
	for(unsigned int i=0;i<atomNames->size();i++) {
		if(this->atomMap.find(atomNames->at(i)) == this->atomMap.end())
			return false;
	}
	return true;
}

int RNABase::printPDBFormat(ofstream& out, int startAtomID) const{
    char c = this->baseID.at(baseID.length()-1);
    char s[100];
    vector<Atom*>::const_iterator it;
    int atomID = startAtomID;
    if(c >= '0' && c <= '9')
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s  %c %c%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->baseType,this->chainID[0],this->baseID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    else
    {
        for(it=atomList.begin();it!=atomList.end();it++)
        {
        	(*it)->guessAtomType();
            XYZ& coord = (*it)->getCoord();
            sprintf(s,"ATOM%7d %5s  %c %c%5s   %8.3f%8.3f%8.3f  1.00  0.00          %2s",atomID,(*it)->nameString().c_str(),this->baseType,this->chainID[0],this->baseID.c_str(),coord[0],coord[1],coord[2],(*it)->type.c_str());
            out << s << endl;
            atomID++;
        }
    }
    return atomID;
}

RNABase::~RNABase() {
}

PolarAtom::PolarAtom(Residue* res, string atomName) {

	this->uniqueName = res->triName + "-" + atomName;
	this->isDonor = false;
	this->isAcceptor = false;

	if(uniqueName == "ASP-OD1" || uniqueName == "ASP-OD2" || uniqueName == "ASN-OD1")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "ASN-ND2")
	{
		this->support = res->getAtom("CG")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "GLU-OE1" || uniqueName == "GLU-OE2" || uniqueName == "GLN-OE1")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isAcceptor = true;
	}
	else if(uniqueName == "GLN-NE2")
	{
		this->support = res->getAtom("CD")->getCoord();
		this->core = res->getAtom(atomName)->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "HIS-ND1")
	{
		XYZ n =  res->getAtom("ND1")->getCoord();
		XYZ sup1 = res->getAtom("CG")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();

		this->core = n;
		this->support = sup1 + sup2 -n;

		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;

	}
	else if(uniqueName == "HIS-NE2")
	{
		XYZ n = res->getAtom("NE2")->getCoord();

		XYZ sup1 = res->getAtom("CD2")->getCoord();
		XYZ sup2 = res->getAtom("CE1")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isAcceptor = true;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "LYS-NZ")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CE")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "ARG-NE")
	{
		XYZ n = res->getAtom("NE")->getCoord();
		XYZ sup1 = res->getAtom("CD")->getCoord();
		XYZ sup2 = res->getAtom("CZ")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "ARG-NH1" || uniqueName == "ARG-NH2")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isDonor = true;
	}
	else if(uniqueName == "SER-OG" || uniqueName == "THR-OG1")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CB")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(uniqueName == "TRP-NE1")
	{
		XYZ n = res->getAtom("NE1")->getCoord();
		XYZ sup1 = res->getAtom("CD1")->getCoord();
		XYZ sup2 = res->getAtom("CE2")->getCoord();
		this->core = n;
		this->support = sup1 + sup2 -n;
		this->isDonor = true;
	//	this->HCore = true;
	}
	else if(uniqueName == "TYR-OH")
	{
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("CZ")->getCoord();
		this->isAcceptor = true;
		this->isDonor = true;
	}
	else if(atomName == "O" || atomName == "OXT" || atomName == "OT1" || atomName == "OT2")
	{
		this->uniqueName = res->triName + "-O";
		this->core = res->getAtom(atomName)->getCoord();
		this->support = res->getAtom("C")->getCoord();
		this->isAcceptor = true;
	}
	else
	{
		cerr << "not polar atom " << uniqueName << endl;
		exit(0);
	}

	if(atomName.at(0) == 'O')
		this->vdwRadius = 1.6;
	else
		this->vdwRadius = 1.7;

}

PolarAtom::PolarAtom(RNABase* base, string atomName){

	char xx[20];
	this->isAcceptor = false;
	this->isDonor = false;
	this->core = XYZ(0, 0, 0);
	this->support = XYZ(0, 0, 0);

	sprintf(xx, "%c-%s", base->baseType, atomName.c_str());
	this->uniqueName = string(xx);
	if(base->baseTypeInt == 0){
		if(atomName == "N1"){
			XYZ c = base->getAtom("N1")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N3"){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N6"){
			XYZ c = base->getAtom("N6")->coord;
			XYZ sup = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup;
			this->isDonor = true;
		}
		else if(atomName == "N7"){
			XYZ c = base->getAtom("N7")->coord;
			XYZ sup1 = base->getAtom("C5")->coord;
			XYZ sup2 = base->getAtom("C8")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
	}
	else if(base->baseTypeInt == 1){
		if(atomName == "O2"){
			XYZ c = base->getAtom("O2")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup1;
			this->isAcceptor = true;
		}
		else if(atomName == "N3"){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "O4"){
			XYZ c = base->getAtom("O4")->coord;
			XYZ sup = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
	}
	else if(base->baseTypeInt == 2){
		if(atomName == "N1"){
			XYZ c = base->getAtom("N1")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N2"){
			XYZ c = base->getAtom("N2")->coord;
			XYZ sup = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup;
			this->isDonor = true;
		}
		else if(atomName == "N3"){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "O6"){
			XYZ c = base->getAtom("O6")->coord;
			XYZ sup = base->getAtom("C6")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
		else if(atomName == "N7"){
			XYZ c = base->getAtom("N7")->coord;
			XYZ sup1 = base->getAtom("C5")->coord;
			XYZ sup2 = base->getAtom("C8")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
	}
	else if(base->baseTypeInt == 3){
		if(atomName == "O2"){
			XYZ c = base->getAtom("O2")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			this->core = c;
			this->support = sup1;
			this->isAcceptor = true;
		}
		else if(atomName == "N3"){
			XYZ c = base->getAtom("N3")->coord;
			XYZ sup1 = base->getAtom("C2")->coord;
			XYZ sup2 = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup1 + sup2 - c;
			this->isDonor = true;
		}
		else if(atomName == "N4"){
			XYZ c = base->getAtom("N4")->coord;
			XYZ sup = base->getAtom("C4")->coord;
			this->core = c;
			this->support = sup;
			this->isAcceptor = true;
		}
	}
	if(atomName[0] == 'N')
		this->vdwRadius = 1.7;
	else
		this->vdwRadius = 1.6;

}

PolarAtom::PolarAtom(Residue* res, string atomName, Atom* preC)
{
	if(atomName != "N" || res->triName == "PRO")
	{
		cerr << "constructor for only backbone polar N" << res->triName << "-" <<  atomName << endl;
		exit(0);
	}
	this->uniqueName = res->triName + "-N";
	XYZ n = res->getAtom("N")->getCoord();
	XYZ sup1 = res->getAtom("CA")->getCoord();
	XYZ sup2 = preC->getCoord();
	this->core = n;
	this->support = sup1 + sup2 -n;
	this->isDonor = true;
	this->isAcceptor = false;
	this->vdwRadius = 1.7;
//	this->HCore = true;
}

PolarAtom::~PolarAtom() {
}

ProteinChain::ProteinChain() {
	this->pdbID = "xxxx";
	this->chainID = '-';
	this->chainLen = 0;
}

ProteinChain::ProteinChain(char chainID) {
	this->pdbID = "xxxx";
	this->chainID = chainID;
	this->chainLen = 0;
}

ProteinChain::ProteinChain(string pdbID, char chainID){
	this->pdbID = pdbID;
	this->chainID = chainID;
	this->chainLen = 0;
}

void ProteinChain::setPDBID(string pdbID){
	this->pdbID = pdbID;
}

void ProteinChain::setChainID(char c){
	this->chainID = c;
}

string ProteinChain::getPDBID() const{
	return this->pdbID;
}

char ProteinChain::getChainID() const{
	return this->chainID;
}

int ProteinChain::getChainLength() const{
	return this->chainLen;
}

vector<Residue*>& ProteinChain::getResList(){
	return this->resList;
}

Residue* ProteinChain::getResidue(const string& resID) {
	map<string,Residue*>::const_iterator it = resMap.find(resID);
	if(it != resMap.end())
		return it->second;
	else
		return NULL;
}

void ProteinChain::addResidue(Residue* res){
	this->resList.push_back(res);
	this->resMap[res->resID] = res;
	res->setResSeqID(this->resList.size()-1);
	this->chainLen ++;
}

string ProteinChain::getSequence() const{
	string s = "";
	ResName rn = ResName();
	for(int i=0;i<chainLen;i++){
		char c = rn.triToSin(this->resList.at(i)->triName);
		s = s + c;
	}
	return s;
}

int ProteinChain::printPDBFormat(ofstream& out, int startAtomID) const{

	for(int i=0;i<resList.size();i++){
		startAtomID = resList.at(i)->printPDBFormat(out,startAtomID);
	}
	return startAtomID;
}

int ProteinChain::printPDBFormatNoHydrogen(ofstream& out, int startAtomID) const{

	return 0;
}


ProteinChain::~ProteinChain() {
}

RNAChain::RNAChain() {
	this->pdbID = "pdbx";
	this->chainID = "A";
	this->chainLen = 0;
}


RNAChain::RNAChain(const string& chainID) {
	this->pdbID = "pdbx";
	this->chainID = chainID;
	this->chainLen = 0;
}

RNAChain::RNAChain(const string& pdbID, const string& chainID) {
	this->pdbID = pdbID;
	this->chainID = chainID;
	this->chainLen = 0;
}

int RNAChain::printPDBFormat(ofstream& out, int startAtomID) const{
	for(int i=0;i<baseList.size();i++) {
		startAtomID = baseList[i]->printPDBFormat(out, startAtomID);
	}
	return startAtomID;
}

PDB::PDB() {
	this->pdbID = "pdbx";
}

PDB::PDB(const string& pdbFile, const string& pdbID)
{
	this->pdbID = pdbID;
	ifstream input;
	input.open(pdbFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << pdbFile << endl;
        exit (1);
    }
    string s;
    int models = 0;
    int len;
    char lastChainID = '@';
    string lastResID = "XXX";

    char curChainID;
    string curResID;
    char altLoc;
    string resName;

    ProteinChain* curChain;
    Residue* curResidue;

    ResName rn = ResName();
	while(getline(input,s))
    {
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        resName = s.substr(17,3);
        if(resName == "HOH") continue;

        if(!rn.isStandardAminoAcid(resName) && resName != "MSE") continue;


        curChainID = s.at(21);
        curResID = trimString(s.substr(22,5));
        altLoc = s.at(16);


        Atom* a = new Atom(s);
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new ProteinChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

                  lastChainID = curChainID;
                  lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new Residue(curResID,curChainID,resName);
            curChain->addResidue(curResidue);
            residues.push_back(curResidue);
            lastResID = curResID;
        }
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') {
        	curResidue->hasAltConf = true;
        	continue;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }

    input.close();
}

PDB& PDB::operator=(const PDB& other){
	cout << "operator '=' is inhibited in class PDB" << endl;
	abort();
	return *this;
}

vector<ProteinChain*>& PDB::getChains()
{
    return chains;
}

ProteinChain* PDB::getFirstChain()
{
    if(this->chains.size() > 0)
        return this->chains.at(0);
    return NULL;
}

ProteinChain* PDB::getChain(char c)
{
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        if(this->chains.at(i)->getChainID() == c)
            return chains.at(i);
    }
    return NULL;
}

string PDB::getFirstSeq()
{
    return getFirstChain()->getSequence();
}

vector<Residue*>& PDB::getResList()
{
    return residues;
}

void PDB::printPDBFormat(ofstream& out) const
{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
        ProteinChain* pc = this->chains.at(i);
        for(int j=0;j<pc->getChainLength();j++)
        {
            Residue* res = pc->getResList().at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}


PDB::~PDB() {
	unsigned int i,j;
	ProteinChain* p;
	for(i = 0;i<chains.size();i++)
	{
		p = chains.at(i);
		delete p;

	}

	Residue* p2;
	Atom* p3;
	for(i=0;i<residues.size();i++)
	{

		p2 = residues.at(i);

		for(j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;

	}
}

RNAPDB::RNAPDB() {
	this->pdbID = "pdbx";
}

RNAPDB::RNAPDB(const string& pdbFile){
	this->pdbID = "pdbx";
	if(pdbFile.substr(pdbFile.length()-3, 3) == "pdb"){
		readPDB(pdbFile);
	}
	else if(pdbFile.substr(pdbFile.length()-3, 3) == "cif"){
		readCIF(pdbFile);
	}
}

RNAPDB::RNAPDB(const string& pdbFile, const string& pdbID) {

	this->pdbID = pdbID;
	if(pdbFile.substr(pdbFile.length()-3, 3) == "pdb"){
		readPDB(pdbFile);
	}
	else if(pdbFile.substr(pdbFile.length()-3, 3) == "cif"){
		readCIF(pdbFile);
	}
}

void RNAPDB::readPDB(const string& pdbFile){
	ifstream input;
	input.open(pdbFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << pdbFile << endl;
        exit (1);
    }
    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID;
    string curResID;
    int curResSeqID = 0;
    char altLoc;
    string rawResName;
    string resName;

    RNAChain* curChain;
    RNABase* curResidue;
    ResName rn;

	while(getline(input,s))
    {
		//cout << s << endl;
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;
        rawResName = s.substr(17,3);
        if(rawResName == "HOH") continue;

        rawResName = trimString(rawResName);
        if(!rn.isRNABase(rawResName)) continue;
        resName = rn.toStandardBase(rawResName);


        curChainID = s.substr(21,1);
        curResID = trimString(s.substr(22,5));
        altLoc = s.at(16);


        Atom* a = new Atom(s);
        if(rawResName == "4SU" && a->name == "S4"){
        	a->name = "O4";
        	a->type = "O";
        }
        if(rawResName == "QUO" && a->name == "C7"){
        	a->name = "N7";
        	a->type = "N";
        }

        if(a->type == "H") continue;
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new RNAChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

             lastChainID = curChainID;
             lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new RNABase(curResID, curChainID, resName[0]);
        	curResidue->setResSeqID(curResSeqID);
        	curResSeqID++;
        	curChain->addBase(curResidue);
            baseList.push_back(curResidue);
            lastResID = curResID;
        }
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') {
        	curResidue->hasAltConf = true;
        	continue;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }

    input.close();
}

void RNAPDB::readCIF(const string& cifFile){
	ifstream input;
	input.open(cifFile.c_str(),ios::in);

	if (! input.is_open())
    {
        cout << "fail to open file " << cifFile << endl;
        exit (1);
    }

    string s;
    int models = 0;
    int len;
    string lastChainID = "@";
    string lastResID = "XXX";

    string curChainID;
    string curResID;
    int curResSeqID = 0;
    char altLoc;
    string rawResName;
    string resName;

    RNAChain* curChain;
    RNABase* curResidue;
    ResName rn;

    vector<string> spt;

	while(getline(input,s))
    {
		//cout << s << endl;
        len = s.length();
        if(len < 6) continue;
        string prefix = s.substr(0,6);
        if(prefix == "MODEL ")
        {
            if(models > 0)
                break;
            else
                models = 1;
        }
        if(len < 54) continue;
        if(prefix != "ATOM  " && prefix != "HETATM") continue;

        Atom* a = new Atom(s, 1);
        splitString(s, " ", &spt);
        if(spt.size() < 19) continue;

        rawResName = spt[5];
        if(rawResName == "HOH") continue;

        rawResName = trimString(rawResName);
        if(!rn.isRNABase(rawResName)) continue;
        resName = rn.toStandardBase(rawResName);


        curChainID = spt[18];

        curResID = spt[16];
        altLoc = spt[9][0];

        if(altLoc == '?')
        	altLoc = ' ';

        if(rawResName == "4SU" && a->name == "S4"){
        	a->name = "O4";
        	a->type = "O";
        }
        if(rawResName == "QUO" && a->name == "C7"){
        	a->name = "N7";
        	a->type = "N";
        }

        if(a->type == "H") continue;
        if(curChainID != lastChainID)
        {
             if(getChain(curChainID) == NULL)
             {
                  curChain = new RNAChain(curChainID);
                  this->chains.push_back(curChain);
             }
             else
            	  curChain = getChain(curChainID);

             lastChainID = curChainID;
             lastResID = "XXX";
        }

        if(curResID != lastResID)
        {
        	curResidue = new RNABase(curResID, curChainID, resName[0]);
        	curResidue->setResSeqID(curResSeqID);
        	curResSeqID++;
        	curChain->addBase(curResidue);
            baseList.push_back(curResidue);
            lastResID = curResID;
        }
        if(altLoc != ' ' && altLoc != 'A' && altLoc != '1') {
        	curResidue->hasAltConf = true;
        	continue;
        }

        curResidue->addAtom(a);
        curResidue->setAltLoc(altLoc);
    }

    input.close();
}



void RNAPDB::printPDBFormat(ofstream& out) const{
    int startID = 1;
    for(unsigned int i=0;i<this->chains.size();i++)
    {
    	RNAChain* pc = this->chains.at(i);
        for(int j=0;j<pc->getChainLength();j++)
        {
            RNABase* res = pc->getBaseList().at(j);
            startID = res->printPDBFormat(out,startID);
        }
    }
}

RNAPDB::~RNAPDB() {
	unsigned int i,j;
	RNAChain* p;
	for(i = 0;i<chains.size();i++)
	{
		p = chains.at(i);
		delete p;

	}

	RNABase* p2;
	Atom* p3;
	for(i=0;i<baseList.size();i++)
	{

		p2 = baseList.at(i);

		for(j=0;j<p2->getAtomList()->size();j++)
		{
			p3 = p2->getAtomList()->at(j);
			delete p3;

		}
		delete p2;
	}
}


float Phipsi::distance(const Phipsi& other) const{
	float delX = this->phi - other.phi;
	float delY = this->psi - other.psi;
	if(delX > 180)
		delX = 360 - delX;
	else if(delX < -180)
		delX = 360 + delX;

	if(delY > 180)
		delY = 360 - delY;
	else if(delY < -180)
		delY = 360 + delY;
	return sqrt(delX*delX+delY*delY);
}

char Phipsi::regionAB() const{
	if(phi < 0 && psi > -100 && psi < 60)
	{
		return 'A';
	}
	return 'B';
}

Phipsi::~Phipsi() {
}

PhipsiLib::PhipsiLib(){
	this->pointNum = 200;
	this->ppList.push_back(new Phipsi( -62.71 ,  -41.57));
	this->ppList.push_back(new Phipsi( -60.44 ,  -36.89));
	this->ppList.push_back(new Phipsi( -65.79 ,  -36.65));
	this->ppList.push_back(new Phipsi( -67.17 ,  -41.65));
	this->ppList.push_back(new Phipsi( -59.91 ,  -45.30));
	this->ppList.push_back(new Phipsi( -56.70 ,  -41.11));
	this->ppList.push_back(new Phipsi( -64.61 ,  -46.65));
	this->ppList.push_back(new Phipsi( -59.68 ,  -50.44));
	this->ppList.push_back(new Phipsi( -54.61 ,  -46.94));
	this->ppList.push_back(new Phipsi( -71.59 ,  -46.40));
	this->ppList.push_back(new Phipsi( -72.03 ,  -37.93));
	this->ppList.push_back(new Phipsi( -62.47 ,  -30.49));
	this->ppList.push_back(new Phipsi( -54.86 ,  -32.60));
	this->ppList.push_back(new Phipsi( -69.62 ,  -31.26));
	this->ppList.push_back(new Phipsi( -48.60 ,  -39.19));
	this->ppList.push_back(new Phipsi( -65.64 ,  -55.36));
	this->ppList.push_back(new Phipsi( -53.82 ,  -55.75));
	this->ppList.push_back(new Phipsi( -46.67 ,  -48.77));
	this->ppList.push_back(new Phipsi( -57.50 ,  -24.42));
	this->ppList.push_back(new Phipsi( -66.32 ,  -23.90));
	this->ppList.push_back(new Phipsi( -78.56 ,  -32.35));
	this->ppList.push_back(new Phipsi( -81.49 ,  -42.96));
	this->ppList.push_back(new Phipsi( -75.23 ,  -23.39));
	this->ppList.push_back(new Phipsi( -62.73 ,  -17.29));
	this->ppList.push_back(new Phipsi( -70.93 ,  -15.55));
	this->ppList.push_back(new Phipsi( -80.44 ,  -56.71));
	this->ppList.push_back(new Phipsi( -86.95 ,  -23.82));
	this->ppList.push_back(new Phipsi( -92.42 ,  -35.65));
	this->ppList.push_back(new Phipsi( -81.19 ,  -14.33));
	this->ppList.push_back(new Phipsi( -66.67 ,   -7.74));
	this->ppList.push_back(new Phipsi( -76.53 ,   -6.37));
	this->ppList.push_back(new Phipsi( -96.72 ,  -50.52));
	this->ppList.push_back(new Phipsi( -91.95 ,  -12.06));
	this->ppList.push_back(new Phipsi(-100.84 ,  -23.37));
	this->ppList.push_back(new Phipsi( -85.93 ,   -3.37));
	this->ppList.push_back(new Phipsi(-108.71 ,  -38.01));
	this->ppList.push_back(new Phipsi( -76.47 ,    6.23));
	this->ppList.push_back(new Phipsi(-104.47 ,  -10.15));
	this->ppList.push_back(new Phipsi( -96.36 ,   -0.46));
	this->ppList.push_back(new Phipsi( -89.35 ,    6.81));
	this->ppList.push_back(new Phipsi(-116.16 ,  -21.63));
	this->ppList.push_back(new Phipsi( -38.10 ,  -59.62));
	this->ppList.push_back(new Phipsi(-107.80 ,    3.09));
	this->ppList.push_back(new Phipsi(-100.12 ,   11.45));
	this->ppList.push_back(new Phipsi(-117.52 ,   -6.73));
	this->ppList.push_back(new Phipsi(-117.12 ,  -57.58));
	this->ppList.push_back(new Phipsi( -95.66 ,  -77.23));
	this->ppList.push_back(new Phipsi( -91.15 ,   20.75));
	this->ppList.push_back(new Phipsi(-129.66 ,  -36.11));
	this->ppList.push_back(new Phipsi(-112.26 ,   14.81));
	this->ppList.push_back(new Phipsi(-133.22 ,  -11.58));
	this->ppList.push_back(new Phipsi(-124.53 ,    6.65));
	this->ppList.push_back(new Phipsi(-105.87 ,   25.91));
	this->ppList.push_back(new Phipsi( -58.02 ,  -83.01));
	this->ppList.push_back(new Phipsi(-124.36 ,   23.11));
	this->ppList.push_back(new Phipsi( -78.31 ,   41.41));
	this->ppList.push_back(new Phipsi(-142.95 ,  -64.95));
	this->ppList.push_back(new Phipsi(-143.81 ,   13.38));
	this->ppList.push_back(new Phipsi(-119.01 ,   38.30));
	this->ppList.push_back(new Phipsi(-119.93 ,  -94.85));
	this->ppList.push_back(new Phipsi(-138.92 ,   38.56));
	this->ppList.push_back(new Phipsi(-172.29 ,  -27.81));
	this->ppList.push_back(new Phipsi(  -5.37 ,  -81.80));
	this->ppList.push_back(new Phipsi(-109.17 ,  132.79));
	this->ppList.push_back(new Phipsi(-108.81 ,  123.49));
	this->ppList.push_back(new Phipsi(-117.57 ,  127.44));
	this->ppList.push_back(new Phipsi(-118.05 ,  136.97));
	this->ppList.push_back(new Phipsi(-118.00 ,  118.36));
	this->ppList.push_back(new Phipsi(-125.28 ,  131.48));
	this->ppList.push_back(new Phipsi(-127.56 ,  122.79));
	this->ppList.push_back(new Phipsi(-108.62 ,  114.24));
	this->ppList.push_back(new Phipsi(-100.17 ,  128.43));
	this->ppList.push_back(new Phipsi( -99.45 ,  118.02));
	this->ppList.push_back(new Phipsi(-110.11 ,  143.24));
	this->ppList.push_back(new Phipsi(-100.40 ,  139.45));
	this->ppList.push_back(new Phipsi(-128.65 ,  139.96));
	this->ppList.push_back(new Phipsi(-121.40 ,  146.83));
	this->ppList.push_back(new Phipsi(-134.95 ,  131.25));
	this->ppList.push_back(new Phipsi(-127.40 ,  111.21));
	this->ppList.push_back(new Phipsi(-115.32 ,  105.65));
	this->ppList.push_back(new Phipsi(-138.60 ,  119.63));
	this->ppList.push_back(new Phipsi(-138.74 ,  141.26));
	this->ppList.push_back(new Phipsi(-131.24 ,  150.02));
	this->ppList.push_back(new Phipsi(-113.22 ,  153.33));
	this->ppList.push_back(new Phipsi( -91.01 ,  134.29));
	this->ppList.push_back(new Phipsi( -90.81 ,  123.32));
	this->ppList.push_back(new Phipsi(-101.02 ,  152.94));
	this->ppList.push_back(new Phipsi(-123.72 ,  157.97));
	this->ppList.push_back(new Phipsi( -90.39 ,  145.82));
	this->ppList.push_back(new Phipsi(-100.97 ,  105.75));
	this->ppList.push_back(new Phipsi( -90.37 ,  110.74));
	this->ppList.push_back(new Phipsi( -81.99 ,  129.04));
	this->ppList.push_back(new Phipsi( -80.58 ,  139.50));
	this->ppList.push_back(new Phipsi( -81.23 ,  117.30));
	this->ppList.push_back(new Phipsi(-110.58 ,  164.44));
	this->ppList.push_back(new Phipsi( -88.14 ,  157.53));
	this->ppList.push_back(new Phipsi( -78.65 ,  149.95));
	this->ppList.push_back(new Phipsi( -73.22 ,  125.73));
	this->ppList.push_back(new Phipsi( -96.67 ,  167.33));
	this->ppList.push_back(new Phipsi( -71.22 ,  135.23));
	this->ppList.push_back(new Phipsi( -70.37 ,  144.53));
	this->ppList.push_back(new Phipsi(-134.04 ,  160.70));
	this->ppList.push_back(new Phipsi(-140.74 ,  152.22));
	this->ppList.push_back(new Phipsi(-123.92 ,  169.98));
	this->ppList.push_back(new Phipsi( -76.83 ,  160.13));
	this->ppList.push_back(new Phipsi(-146.93 ,  130.97));
	this->ppList.push_back(new Phipsi( -81.98 ,  169.43));
	this->ppList.push_back(new Phipsi( -68.82 ,  154.34));
	this->ppList.push_back(new Phipsi(-108.56 , -179.84));
	this->ppList.push_back(new Phipsi(-149.54 ,  142.95));
	this->ppList.push_back(new Phipsi(-136.72 ,  172.17));
	this->ppList.push_back(new Phipsi(-145.19 ,  163.00));
	this->ppList.push_back(new Phipsi(-151.06 ,  154.55));
	this->ppList.push_back(new Phipsi(-108.83 ,   93.20));
	this->ppList.push_back(new Phipsi( -68.43 ,  116.31));
	this->ppList.push_back(new Phipsi( -78.79 ,  103.37));
	this->ppList.push_back(new Phipsi( -63.16 ,  128.77));
	this->ppList.push_back(new Phipsi( -90.94 ,   95.27));
	this->ppList.push_back(new Phipsi( -61.74 ,  138.40));
	this->ppList.push_back(new Phipsi( -61.79 ,  147.62));
	this->ppList.push_back(new Phipsi( -68.50 ,  165.94));
	this->ppList.push_back(new Phipsi(-126.07 ,   94.90));
	this->ppList.push_back(new Phipsi( -90.20 , -176.79));
	this->ppList.push_back(new Phipsi(-141.05 ,  103.07));
	this->ppList.push_back(new Phipsi( -58.84 ,  157.36));
	this->ppList.push_back(new Phipsi( -53.54 ,  132.85));
	this->ppList.push_back(new Phipsi( -52.94 ,  143.20));
	this->ppList.push_back(new Phipsi( -74.20 ,  179.78));
	this->ppList.push_back(new Phipsi( -54.33 ,  121.37));
	this->ppList.push_back(new Phipsi(-126.43 , -170.97));
	this->ppList.push_back(new Phipsi(-155.86 ,  116.09));
	this->ppList.push_back(new Phipsi(-160.77 ,  148.35));
	this->ppList.push_back(new Phipsi(-162.21 ,  134.10));
	this->ppList.push_back(new Phipsi(-153.57 ,  171.62));
	this->ppList.push_back(new Phipsi(-159.20 ,  161.49));
	this->ppList.push_back(new Phipsi(-145.25 , -175.88));
	this->ppList.push_back(new Phipsi( -80.19 ,   82.18));
	this->ppList.push_back(new Phipsi(-122.33 ,   77.75));
	this->ppList.push_back(new Phipsi( -93.64 ,   75.65));
	this->ppList.push_back(new Phipsi(-138.24 ,   79.85));
	this->ppList.push_back(new Phipsi( -54.25 ,   96.28));
	this->ppList.push_back(new Phipsi( -40.61 ,  127.84));
	this->ppList.push_back(new Phipsi(-104.40 , -156.29));
	this->ppList.push_back(new Phipsi( -81.68 , -159.28));
	this->ppList.push_back(new Phipsi(-173.28 ,  156.71));
	this->ppList.push_back(new Phipsi(-167.83 ,  170.51));
	this->ppList.push_back(new Phipsi(-159.27 ,   91.56));
	this->ppList.push_back(new Phipsi(-162.77 , -174.54));
	this->ppList.push_back(new Phipsi(-148.38 , -155.90));
	this->ppList.push_back(new Phipsi( -81.65 ,   64.08));
	this->ppList.push_back(new Phipsi(-130.51 ,   61.59));
	this->ppList.push_back(new Phipsi(-152.79 ,   64.32));
	this->ppList.push_back(new Phipsi(-100.02 ,   51.12));
	this->ppList.push_back(new Phipsi(-122.55 , -135.41));
	this->ppList.push_back(new Phipsi( -92.38 , -124.40));
	this->ppList.push_back(new Phipsi(-157.80 , -114.40));
	this->ppList.push_back(new Phipsi(  79.97 ,   11.05));
	this->ppList.push_back(new Phipsi(  67.21 ,   13.11));
	this->ppList.push_back(new Phipsi(  75.20 ,   24.73));
	this->ppList.push_back(new Phipsi(  77.13 ,   -1.28));
	this->ppList.push_back(new Phipsi(  92.68 ,   16.58));
	this->ppList.push_back(new Phipsi(  91.44 ,    1.17));
	this->ppList.push_back(new Phipsi(  60.91 ,   25.89));
	this->ppList.push_back(new Phipsi(  67.07 ,   37.97));
	this->ppList.push_back(new Phipsi(  55.06 ,   37.75));
	this->ppList.push_back(new Phipsi(  87.84 ,  -12.94));
	this->ppList.push_back(new Phipsi( 111.79 ,    6.45));
	this->ppList.push_back(new Phipsi( 102.18 ,  -12.68));
	this->ppList.push_back(new Phipsi(  97.17 ,   46.11));
	this->ppList.push_back(new Phipsi(  58.70 ,   52.03));
	this->ppList.push_back(new Phipsi(  48.74 ,   46.96));
	this->ppList.push_back(new Phipsi(  43.41 ,   60.83));
	this->ppList.push_back(new Phipsi(  62.67 ,   78.25));
	this->ppList.push_back(new Phipsi( 104.48 ,  -30.61));
	this->ppList.push_back(new Phipsi(  73.14 ,  -44.39));
	this->ppList.push_back(new Phipsi( 131.11 ,  -17.09));
	this->ppList.push_back(new Phipsi( 160.05 ,   42.52));
	this->ppList.push_back(new Phipsi( 142.52 ,  -65.36));
	this->ppList.push_back(new Phipsi(   4.13 ,   98.81));
	this->ppList.push_back(new Phipsi( 123.61 , -160.45));
	this->ppList.push_back(new Phipsi(  92.55 , -154.26));
	this->ppList.push_back(new Phipsi( 100.14 ,  178.37));
	this->ppList.push_back(new Phipsi( 128.04 ,  168.58));
	this->ppList.push_back(new Phipsi(  76.51 , -174.33));
	this->ppList.push_back(new Phipsi(  84.21 ,  162.79));
	this->ppList.push_back(new Phipsi( 103.76 ,  146.02));
	this->ppList.push_back(new Phipsi(  67.39 , -154.11));
	this->ppList.push_back(new Phipsi( 152.15 , -168.57));
	this->ppList.push_back(new Phipsi( 114.53 , -122.60));
	this->ppList.push_back(new Phipsi( 146.81 , -140.10));
	this->ppList.push_back(new Phipsi(  77.05 , -123.32));
	this->ppList.push_back(new Phipsi(  56.48 , -134.67));
	this->ppList.push_back(new Phipsi( 159.19 ,  165.30));
	this->ppList.push_back(new Phipsi( 175.66 , -154.48));
	this->ppList.push_back(new Phipsi( 176.13 , -177.66));
	this->ppList.push_back(new Phipsi( 147.35 ,  128.99));
	this->ppList.push_back(new Phipsi(  52.69 , -117.81));
	this->ppList.push_back(new Phipsi( 106.36 ,  108.41));
	this->ppList.push_back(new Phipsi(  68.96 ,  118.77));
	this->ppList.push_back(new Phipsi(  81.03 ,  -73.88));
	creatIndexTable();
}
int PhipsiLib::findNearestPointsWithoutIndex(Phipsi* pp)
{
	float minDis = 10000.0;
	int minIndex = -1;
	float d;
	for(int i=0;i<pointNum;i++)
	{
		d = pp->distance(*(ppList.at(i)));
		//cout << d << endl;
		if(d < minDis)
		{
			minDis = d;
			minIndex = i;
		}
	}
	return minIndex;
}

void PhipsiLib::creatIndexTable()
{
	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			for(int k=0;k<20;k++)
			{
				this->indexTable[i][j][k] = -1;
			}
		}
	}

	set<unsigned int> possibleNeighbor;
	Phipsi *pp, *pp1, *pp2, *pp3, *pp4;
	set<unsigned int>::const_iterator it;

	for(int i=0;i<36;i++)
	{
		for(int j=0;j<36;j++)
		{
			possibleNeighbor.clear();
			float phi0 = i*10-180;
			float psi0 = j*10-180;
			for(unsigned int a=0;a<ppList.size();a++)
			{
				pp = ppList.at(a);
				if((pp->phi) > phi0 && (pp->phi) < phi0+10 && (pp->psi) > psi0 && (pp->psi) < psi0+10)
					possibleNeighbor.insert(a);
			}

			for(int a=0;a<11;a++)
			{
				pp1 = new Phipsi(phi0+a,psi0);
				pp2 = new Phipsi(phi0+a, psi0+10);
				pp3 = new Phipsi(phi0, psi0+a);
				pp4 = new Phipsi(phi0+10, psi0+a);
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp1));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp2));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp3));
				possibleNeighbor.insert(findNearestPointsWithoutIndex(pp4));
				delete pp1;
				delete pp2;
				delete pp3;
				delete pp4;
			}

			int n=0;
			for(it = possibleNeighbor.begin();it!=possibleNeighbor.end();it++)
			{
				int k = *it;
				this->indexTable[i][j][n] = k;
				n++;
			}
		}
	}
}

Phipsi* PhipsiLib::indexToPhipsi(int id) const
{
	return this->ppList.at(id);
}

int PhipsiLib::phipsiToIndex(const Phipsi* pp) const
{
	int i = (int)((pp->phi+180)/10);
	int j = (int)((pp->psi+180)/10);
	double minD = 10000.0;
	int minIndex = -1;
	for(int k=0;k<20;k++)
	{
		int id = this->indexTable[i][j][k];
		if(id < 0)
			break;
		double dist = pp->distance(* ppList.at(id));
		if(dist < minD)
		{
			minD = dist;
			minIndex = id;
		}
	}
	return minIndex;
}

vector<pair<int,double>> PhipsiLib::neighborPhipsiIndexList(const Phipsi* pp) const
{
	int n = 3;
	int indexList[n];
	double distList[n];
	for(int i=0;i<n;i++){
		indexList[i] = -1;
		distList[i] = 999999.9;
	}

	double d;
	int k;
	for(int i=0;i<pointNum;i++){
		d = pp->distance(*ppList[i]);
		//cout << "index: " << i << " distance: " << d << endl;
		if(d == 0)
			d = 0.001;
        if(d < distList[n-1]) {
            distList[n-1] = d;
            indexList[n-1] = i;
        }
        else
            continue;

        for(int j=n-2;j>=0;j--) {
            if(distList[j+1] < distList[j]) {
                d = distList[j];
                distList[j] = distList[j+1];
                distList[j+1] = d;
                k = indexList[j];
                indexList[j] = indexList[j+1];
                indexList[j+1] = k;
            }
            else
                break;
        }
	}

	vector<pair<int,double>> result;
	for(int i=0;i<n;i++){
		pair<int,double> p(indexList[i],distList[i]);
		result.push_back(p);
	}
	return result;
}

PhipsiLib::~PhipsiLib() {
	for(int i=0;i<this->pointNum;i++)
	{
		delete this->ppList.at(i);
	}
}


ResPairOrientation::ResPairOrientation(){
	for(int i=0;i<10;i++){
		points.push_back(XYZ(0,0,0));
	}
}

ResPairOrientation::ResPairOrientation(LocalFrame& csA, LocalFrame& csB){

	CsMove cs = csA.getMove(csB);
	LocalFrame newCsB = LocalFrame(cs.oriMove, cs.tm);
	addAtoms(newCsB);
}

ResPairOrientation::ResPairOrientation(LocalFrame& csB){
	addAtoms(csB);
}


ResPairOrientation::ResPairOrientation(string& s){
	vector<string> spt;
	splitString(s," ",&spt);
	if(spt.size() != 15){
		cout << "orientation string error: " + s << endl;
		exit(1);
	}

	points.push_back(XYZ(0.826, -1.204, 0));
	points.push_back(XYZ(0,0,0));
	points.push_back(XYZ(0.862, 1.257, 0));
	points.push_back(XYZ(-1.040, 0.006, 1.345));
	points.push_back(XYZ(-2.332, -0.663, 0.981));
	points.push_back(XYZ(atof(spt[0].c_str()), atof(spt[1].c_str()), atof(spt[2].c_str())));
	points.push_back(XYZ(atof(spt[3].c_str()), atof(spt[4].c_str()), atof(spt[5].c_str())));
	points.push_back(XYZ(atof(spt[6].c_str()), atof(spt[7].c_str()), atof(spt[8].c_str())));
	points.push_back(XYZ(atof(spt[9].c_str()), atof(spt[10].c_str()), atof(spt[11].c_str())));
	points.push_back(XYZ(atof(spt[12].c_str()), atof(spt[13].c_str()), atof(spt[14].c_str())));

}


void ResPairOrientation::addAtoms(LocalFrame& csB){

	points.push_back(XYZ(0.826, -1.204, 0));
	points.push_back(XYZ(0,0,0));
	points.push_back(XYZ(0.862, 1.257, 0));
	points.push_back(XYZ(-1.040, 0.006, 1.345));
	points.push_back(XYZ(-2.332, -0.663, 0.981));

	points.push_back(csB.local2globalcrd(points[0]));
	points.push_back(csB.local2globalcrd(points[1]));
	points.push_back(csB.local2globalcrd(points[2]));
	points.push_back(csB.local2globalcrd(points[3]));
	points.push_back(csB.local2globalcrd(points[4]));

}


XYZ ResPairOrientation::getTern(int index){
	if(index >=0 && index < 10)
		return points[index];
	else {
		cout << "index out of range" << endl;
		exit(1);
	}
}

string ResPairOrientation::toString() const{
	string rpoString = "";
	char s[30];
	for(int i=0;i<10;i++){
		sprintf(s, " %7.3f %7.3f %7.3f", points[i][0], points[i][1], points[i][2]);
		string ss(s);
		rpoString = rpoString + ss;
	}
	return rpoString;
}




} /* namespace NSPmodel */
