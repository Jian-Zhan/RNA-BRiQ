/*
 * BRNodeBasic.cpp
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#include <pred/BRNodeBasic.h>

namespace NSPpred {

BRNodeBasic& BRNodeBasic::operator =(const BRNodeBasic& other){
	this->baseType = other.baseType;
	this->seqID = other.seqID;
	this->cs1 = other.cs1;
	this->tmpCs1 = other.tmpCs1;

	this->fixed = other.fixed;
	this->father = other.father;
	this->leftChild = other.leftChild;
	this->midChild = other.midChild;
	this->reverseChild = other.reverseChild;
	this->upConnection = other.upConnection;
	this->baseAtomNum = other.baseAtomNum;
	for(int i=0;i<baseAtomNum;i++){
		this->atomCoordLocal[i] = other.atomCoordLocal[i];
		this->baseAtomCoords[i] = other.baseAtomCoords[i];
		this->baseAtomCoordsTmp[i] = other.baseAtomCoordsTmp[i];
	}


	return *this;
}

void BRNodeBasic::copyValueFrom(const BRNodeBasic& other) {
	this->baseType = other.baseType;
	this->seqID = other.seqID;
	this->cs1 = other.cs1;
	this->tmpCs1 = other.tmpCs1;
	this->fixed = other.fixed;
	this->father = other.father;
	this->leftChild = other.leftChild;
	this->midChild = other.midChild;
	this->reverseChild = other.reverseChild;
	this->upConnection = other.upConnection;
	this->baseAtomNum = other.baseAtomNum;
	for(int i=0;i<baseAtomNum;i++){
		this->atomCoordLocal[i] = other.atomCoordLocal[i];
		this->baseAtomCoords[i] = other.baseAtomCoords[i];
		this->baseAtomCoordsTmp[i] = other.baseAtomCoordsTmp[i];
	}

}

vector<Atom*> BRNodeBasic::toAtomList(RnaAtomLib& atLib) {

    vector<Atom*> list;
    vector<string>* names = atLib.getSidechainAtoms(this->baseType);
    vector<XYZ> tList;
    if(baseType == 0){

            XYZ a = XYZ(1.468,   0.000,   0.000); //A-N9
            XYZ b = XYZ(2.306,  -1.084,   0.000); //A-C8
            XYZ c = XYZ(3.577,  -0.770,   0.000); //A-N7
            XYZ d = XYZ(3.576,   0.615,   0.000); //A-C5
            XYZ e = XYZ(4.614,   1.556,   0.000); //A-C6
            XYZ f = XYZ(5.904,   1.230,   0.000); //A-N6
            XYZ g = XYZ(4.276,   2.861,   0.000); //A-N1
            XYZ h = XYZ(2.980,   3.184,   0.000); //A-C2
            XYZ i = XYZ(1.914,   2.391,   0.000); //A-N3
            XYZ j = XYZ(2.285,   1.103,   0.000); //A-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
    }

    else if(baseType == 1){
            XYZ a = XYZ(1.478,   0.000,   0.000); //U-N1
            XYZ b = XYZ(2.122,   1.221,   0.000); //U-C2
            XYZ c = XYZ(1.528,   2.282,   0.000); //U-O2
            XYZ d = XYZ(3.491,   1.159,   0.000); //U-N3
            XYZ e = XYZ(4.265,   0.020,   0.000); //U-C4
            XYZ f = XYZ(5.490,   0.123,   0.000); //U-O4
            XYZ g = XYZ(3.526,  -1.204,   0.000); //U-C5
            XYZ h = XYZ(2.191,  -1.173,   0.000); //U-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }
    else if(baseType == 2) {
            XYZ a = XYZ(1.468,   0.000,   0.000); //G-N9
            XYZ b = XYZ(2.295,  -1.094,   0.000); //G-C8
            XYZ c = XYZ(3.560,  -0.779,   0.000); //G-N7
            XYZ d = XYZ(3.570,   0.606,   0.000); //G-C5
            XYZ e = XYZ(4.655,   1.514,   0.000); //G-C6
            XYZ f = XYZ(5.864,   1.265,   0.000); //G-O6
            XYZ g = XYZ(4.221,   2.832,   0.000); //G-N1
            XYZ h = XYZ(2.909,   3.225,   0.000); //G-C2
            XYZ i = XYZ(2.690,   4.543,   0.000); //G-N2
            XYZ j = XYZ(1.886,   2.389,   0.000); //G-N3
            XYZ k = XYZ(2.287,   1.103,   0.000); //G-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
            tList.push_back(k);
    }
    else if(baseType == 3) {
            XYZ a = XYZ(1.478,   0.000,   0.000); //C-N1
            XYZ b = XYZ(2.151,   1.224,   0.000); //C-C2
            XYZ c = XYZ(1.490,   2.271,   0.000); //C-O2
            XYZ d = XYZ(3.503,   1.239,   0.000); //C-N3
            XYZ e = XYZ(4.178,   0.091,   0.000); //C-C4
            XYZ f = XYZ(5.508,   0.150,   0.000); //C-N4
            XYZ g = XYZ(3.519,  -1.170,   0.000); //C-C5
            XYZ h = XYZ(2.181,  -1.170,   0.000); //C-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }

    vector<Atom*> atomList;
    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), local2global(cs1, tList[i])));
    }
    atomList.push_back(new Atom("C1'", cs1.origin_));

    return atomList;
}

vector<Atom*> BRNodeBasic::toTmpAtomList(RnaAtomLib& atLib) {

    vector<Atom*> list;
    vector<string>* names = atLib.getSidechainAtoms(this->baseType);
    vector<XYZ> tList;
    if(baseType == 0){

            XYZ a = XYZ(1.468,   0.000,   0.000); //A-N9
            XYZ b = XYZ(2.306,  -1.084,   0.000); //A-C8
            XYZ c = XYZ(3.577,  -0.770,   0.000); //A-N7
            XYZ d = XYZ(3.576,   0.615,   0.000); //A-C5
            XYZ e = XYZ(4.614,   1.556,   0.000); //A-C6
            XYZ f = XYZ(5.904,   1.230,   0.000); //A-N6
            XYZ g = XYZ(4.276,   2.861,   0.000); //A-N1
            XYZ h = XYZ(2.980,   3.184,   0.000); //A-C2
            XYZ i = XYZ(1.914,   2.391,   0.000); //A-N3
            XYZ j = XYZ(2.285,   1.103,   0.000); //A-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
    }
    else if(baseType == 1){
            XYZ a = XYZ(1.478,   0.000,   0.000); //U-N1
            XYZ b = XYZ(2.122,   1.221,   0.000); //U-C2
            XYZ c = XYZ(1.528,   2.282,   0.000); //U-O2
            XYZ d = XYZ(3.491,   1.159,   0.000); //U-N3
            XYZ e = XYZ(4.265,   0.020,   0.000); //U-C4
            XYZ f = XYZ(5.490,   0.123,   0.000); //U-O4
            XYZ g = XYZ(3.526,  -1.204,   0.000); //U-C5
            XYZ h = XYZ(2.191,  -1.173,   0.000); //U-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }
    else if(baseType == 2) {
            XYZ a = XYZ(1.468,   0.000,   0.000); //G-N9
            XYZ b = XYZ(2.295,  -1.094,   0.000); //G-C8
            XYZ c = XYZ(3.560,  -0.779,   0.000); //G-N7
            XYZ d = XYZ(3.570,   0.606,   0.000); //G-C5
            XYZ e = XYZ(4.655,   1.514,   0.000); //G-C6
            XYZ f = XYZ(5.864,   1.265,   0.000); //G-O6
            XYZ g = XYZ(4.221,   2.832,   0.000); //G-N1
            XYZ h = XYZ(2.909,   3.225,   0.000); //G-C2
            XYZ i = XYZ(2.690,   4.543,   0.000); //G-N2
            XYZ j = XYZ(1.886,   2.389,   0.000); //G-N3
            XYZ k = XYZ(2.287,   1.103,   0.000); //G-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
            tList.push_back(k);
    }
    else if(baseType == 3) {
            XYZ a = XYZ(1.478,   0.000,   0.000); //C-N1
            XYZ b = XYZ(2.151,   1.224,   0.000); //C-C2
            XYZ c = XYZ(1.490,   2.271,   0.000); //C-O2
            XYZ d = XYZ(3.503,   1.239,   0.000); //C-N3
            XYZ e = XYZ(4.178,   0.091,   0.000); //C-C4
            XYZ f = XYZ(5.508,   0.150,   0.000); //C-N4
            XYZ g = XYZ(3.519,  -1.170,   0.000); //C-C5
            XYZ h = XYZ(2.181,  -1.170,   0.000); //C-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }


    vector<Atom*> atomList;
    for(int i=0;i<tList.size();i++){
    	atomList.push_back(new Atom(names->at(i), local2global(tmpCs1, tList[i])));
    }
    atomList.push_back(new Atom("C1'", tmpCs1.origin_));
    return atomList;
}

BRNodeBasic::~BRNodeBasic() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
