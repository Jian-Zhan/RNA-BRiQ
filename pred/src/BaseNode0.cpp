/*
 * BaseNode.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#include <pred/BaseNode0.h>

namespace NSPpred {

BaseNode::BaseNode() {
	this->baseType = 0;
	this->father = NULL;
	this->seqID = -1;
	this->leftChild = NULL;
	this->rightChild = NULL;
	this->fixed = false;
	// TODO Auto-generated constructor stub
}

vector<Atom*> BaseNode::toAtomList(RnaAtomLib& atLib){
	vector<Atom*> list;
	vector<string>* names = atLib.getSidechainAtoms(this->baseType);
	vector<XYZ> tList;

	if(baseType == 0){

		XYZ a = XYZ(1.468,   0.000,   0.000);
		XYZ b = XYZ(2.306,  -1.084,   0.000);
		XYZ c = XYZ(3.577,  -0.770,   0.000);
		XYZ d = XYZ(3.576,   0.615,   0.000);
		XYZ e = XYZ(4.614,   1.556,   0.000);
		XYZ f = XYZ(5.904,   1.230,   0.000);
		XYZ g = XYZ(4.276,   2.861,   0.000);
		XYZ h = XYZ(2.980,   3.184,   0.000);
		XYZ i = XYZ(1.914,   2.391,   0.000);
		XYZ j = XYZ(2.285,   1.103,   0.000);
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
		XYZ a = XYZ(1.478,   0.000,   0.000);
		XYZ b = XYZ(2.122,   1.221,   0.000);
		XYZ c = XYZ(1.528,   2.282,   0.000);
		XYZ d = XYZ(3.491,   1.159,   0.000);
		XYZ e = XYZ(4.265,   0.020,   0.000);
		XYZ f = XYZ(5.490,   0.123,   0.000);
		XYZ g = XYZ(3.526,  -1.204,   0.000);
		XYZ h = XYZ(2.191,  -1.173,   0.000);
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
		XYZ a = XYZ(1.468,   0.000,   0.000);
		XYZ b = XYZ(2.295,  -1.094,   0.000);
		XYZ c = XYZ(3.560,  -0.779,   0.000);
		XYZ d = XYZ(3.570,   0.606,   0.000);
		XYZ e = XYZ(4.655,   1.514,   0.000);
		XYZ f = XYZ(5.864,   1.265,   0.000);
		XYZ g = XYZ(4.221,   2.832,   0.000);
		XYZ h = XYZ(2.909,   3.225,   0.000);
		XYZ i = XYZ(2.690,   4.543,   0.000);
		XYZ j = XYZ(1.886,   2.389,   0.000);
		XYZ k = XYZ(2.287,   1.103,   0.000);
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
		XYZ a = XYZ(1.478,   0.000,   0.000);
		XYZ b = XYZ(2.151,   1.224,   0.000);
		XYZ c = XYZ(1.490,   2.271,   0.000);
		XYZ d = XYZ(3.503,   1.239,   0.000);
		XYZ e = XYZ(4.178,   0.091,   0.000);
		XYZ f = XYZ(5.508,   0.150,   0.000);
		XYZ g = XYZ(3.519,  -1.170,   0.000);
		XYZ h = XYZ(2.181,  -1.170,   0.000);
		tList.push_back(a);
		tList.push_back(b);
		tList.push_back(c);
		tList.push_back(d);
		tList.push_back(e);
		tList.push_back(f);
		tList.push_back(g);
		tList.push_back(h);
	}
	if(tList.size() != names->size()){
		cout << "atom name error: " << tList.size() << " " << names->size() << endl;
	}
	list.push_back(new Atom("C1'", this->cs.origin_));
	list.push_back(new Atom("P", this->tPho));
	for(int i=0;i<names->size();i++){
		list.push_back(new Atom(names->at(i), this->cs.local2globalcrd(tList[i])));
	}
	//list.push_back(new Atom("P", this->tPho));
	return list;
}

vector<Atom*> BaseNode::toTmpAtomList(RnaAtomLib& atLib){
	vector<Atom*> list;
	vector<string>* names = atLib.getSidechainAtoms(this->baseType);
	vector<XYZ> tList;

	if(baseType == 0){

		XYZ a = XYZ(1.468,   0.000,   0.000);
		XYZ b = XYZ(2.306,  -1.084,   0.000);
		XYZ c = XYZ(3.577,  -0.770,   0.000);
		XYZ d = XYZ(3.576,   0.615,   0.000);
		XYZ e = XYZ(4.614,   1.556,   0.000);
		XYZ f = XYZ(5.904,   1.230,   0.000);
		XYZ g = XYZ(4.276,   2.861,   0.000);
		XYZ h = XYZ(2.980,   3.184,   0.000);
		XYZ i = XYZ(1.914,   2.391,   0.000);
		XYZ j = XYZ(2.285,   1.103,   0.000);
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
		XYZ a = XYZ(1.478,   0.000,   0.000);
		XYZ b = XYZ(2.122,   1.221,   0.000);
		XYZ c = XYZ(1.528,   2.282,   0.000);
		XYZ d = XYZ(3.491,   1.159,   0.000);
		XYZ e = XYZ(4.265,   0.020,   0.000);
		XYZ f = XYZ(5.490,   0.123,   0.000);
		XYZ g = XYZ(3.526,  -1.204,   0.000);
		XYZ h = XYZ(2.191,  -1.173,   0.000);
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
		XYZ a = XYZ(1.468,   0.000,   0.000);
		XYZ b = XYZ(2.295,  -1.094,   0.000);
		XYZ c = XYZ(3.560,  -0.779,   0.000);
		XYZ d = XYZ(3.570,   0.606,   0.000);
		XYZ e = XYZ(4.655,   1.514,   0.000);
		XYZ f = XYZ(5.864,   1.265,   0.000);
		XYZ g = XYZ(4.221,   2.832,   0.000);
		XYZ h = XYZ(2.909,   3.225,   0.000);
		XYZ i = XYZ(2.690,   4.543,   0.000);
		XYZ j = XYZ(1.886,   2.389,   0.000);
		XYZ k = XYZ(2.287,   1.103,   0.000);
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
		XYZ a = XYZ(1.478,   0.000,   0.000);
		XYZ b = XYZ(2.151,   1.224,   0.000);
		XYZ c = XYZ(1.490,   2.271,   0.000);
		XYZ d = XYZ(3.503,   1.239,   0.000);
		XYZ e = XYZ(4.178,   0.091,   0.000);
		XYZ f = XYZ(5.508,   0.150,   0.000);
		XYZ g = XYZ(3.519,  -1.170,   0.000);
		XYZ h = XYZ(2.181,  -1.170,   0.000);
		tList.push_back(a);
		tList.push_back(b);
		tList.push_back(c);
		tList.push_back(d);
		tList.push_back(e);
		tList.push_back(f);
		tList.push_back(g);
		tList.push_back(h);
	}
	if(tList.size() != names->size()){
		cout << "atom name error: " << tList.size() << " " << names->size() << endl;
	}
	list.push_back(new Atom("C1'", this->tmpCs.origin_));
	for(int i=0;i<names->size();i++){
		list.push_back(new Atom(names->at(i), this->tmpCs.local2globalcrd(tList[i])));
	}
	//list.push_back(new Atom("P", this->tPho));
	return list;
}

BaseNode::~BaseNode() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPpred */
