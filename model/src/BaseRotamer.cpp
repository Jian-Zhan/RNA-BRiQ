/*
 * BaseRotamer.cpp
 *
 *  Created on: Jul 17, 2019
 *      Author: s2982206
 */

#include <model/BaseRotamer.h>

namespace NSPmodel {

BaseRotamer::BaseRotamer(RNABase* base) {
	vector<string> nameList;
	nameList.push_back("C1'");
	nameList.push_back("C2'");
	nameList.push_back("C3'");
	nameList.push_back("C4'");
	nameList.push_back("O4'");
	nameList.push_back("O2'");
	nameList.push_back("O3'");
	nameList.push_back("C5'");

	LocalFrame cs = base->getCoordSystem();
	for(int i=0;i<8;i++){
		Atom* a = base->getAtom(nameList[i]);
		if(a == NULL) {
			cout << "base backbone incomplete: " << base->baseID << endl;
			exit(1);
		}
		tList1[i] = cs.global2localcrd(a->coord);
		//cout << "base rotamer: " <<  tList1[i].toString() << endl;
	}

	double x=0,y=0,z=0;
	for(int i=0;i<8;i++){
		x += tList1[i].x_;
		y += tList1[i].y_;
		z += tList1[i].z_;
	}
	center = XYZ(x*0.125, y*0.125, z*0.125);
	this->energy = 0.0;

	LocalFrame cs1;
	LocalFrame cs2 = LocalFrame(tList1[1], tList1[2], tList1[6]);
	LocalFrame cs3 = LocalFrame(tList1[4], tList1[3], tList1[7]);
	mv12 = cs2 - cs1;
	mv13 = cs3 - cs1;
	mv31 = cs1 - cs3;
	mv21 = cs1 - cs2;
	this->rotTypeLv1 = 0;
	this->rotType = 0;
	this->resType = base->baseTypeInt;
	double ang = dihedral(tList1[0], tList1[1], tList1[3], tList1[2]);
	if(ang > 0)
		ang = 180 - ang;
	else
		ang = -ang - 180;
	this->improper = ang;
	XYZ a  = tList1[4];
	XYZ b = tList1[0];
	XYZ c(1.48, 0, 0);
	XYZ d(2.318, -1.073, 0);
	this->chi = dihedral(a,b,c,d);

	XYZ tmp = tList1[1] - tList1[5];
	tList1[8] = tList1[5] + tmp * (1.0/tmp.length());
	tmp = tList1[2] - tList1[6];
	tList1[9] = tList1[6] + tmp * (1.0/tmp.length());
	tmp = tList1[0] - tList1[4] + tList1[3] - tList1[4];
	tList1[10] = tList1[4] + tmp * (1.0/tmp.length());
}

BaseRotamer::BaseRotamer(const string& line) {

	map<char,int> typeToInt;
	typeToInt['A'] = 0;
	typeToInt['U'] = 1;
	typeToInt['G'] = 2;
	typeToInt['C'] = 3;
	this->resType = typeToInt[line[0]];
	vector<string> spt;
	splitString(line," ",&spt);

	tList1[0] = XYZ();

	for(int i=1;i<8;i++){
		tList1[i] = XYZ(atof(spt[i*3-2].c_str()), atof(spt[i*3-1].c_str()), atof(spt[i*3].c_str()));
	}

	double x=0,y=0,z=0;
	for(int i=0;i<8;i++){
		x += tList1[i].x_;
		y += tList1[i].y_;
		z += tList1[i].z_;
	}
	center = XYZ(x*0.125, y*0.125, z*0.125);
	this->energy = atof(spt[22].c_str());

	LocalFrame cs1;
	LocalFrame cs2 = LocalFrame(tList1[1], tList1[2], tList1[6]);
	LocalFrame cs3 = LocalFrame(tList1[4], tList1[3], tList1[7]);
	mv12 = cs2 - cs1;
	mv13 = cs3 - cs1;
	mv31 = cs1 - cs3;
	mv21 = cs1 - cs2;
	this->rotTypeLv1 = 0;
	this->rotType = 0;
	double ang = dihedral(tList1[0], tList1[1], tList1[3], tList1[2]);
	if(ang > 0)
		ang = 180 - ang;
	else
		ang = -ang - 180;
	this->improper = ang;
	XYZ a  = tList1[4];
	XYZ b = tList1[0];
	XYZ c(1.48, 0, 0);
	XYZ d(2.318, -1.073, 0);
	this->chi = dihedral(a,b,c,d);

	XYZ tmp = tList1[1] - tList1[5];
	tList1[8] = tList1[5] + tmp * (1.0/tmp.length());
	tmp = tList1[2] - tList1[6];
	tList1[9] = tList1[6] + tmp * (1.0/tmp.length());
	tmp = tList1[0] - tList1[4] + tList1[3] - tList1[4];
	tList1[10] = tList1[4] + tmp * (1.0/tmp.length());

}

BaseRotamer::~BaseRotamer() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
