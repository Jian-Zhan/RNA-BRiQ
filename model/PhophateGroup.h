/*
 * PhophateGroup.h
 *
 *  Created on: Jul 23, 2019
 *      Author: s2982206
 */

#ifndef MODEL_PHOPHATEGROUP_H_
#define MODEL_PHOPHATEGROUP_H_
#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "geometry/localframe.h"
#include "tools/StringTool.h"

namespace NSPmodel {
using namespace NSPgeometry;
using namespace NSPtools;

class PhophateGroupLocal {
public:
	XYZ p;
	XYZ o5;
	XYZ op1;
	XYZ op2;

	XYZ o5Support;
	XYZ op1Support;
	XYZ op2Support;
	double e;

	PhophateGroupLocal(){
		p = XYZ(0.805 ,  0.299 , -1.356);
		o5 = XYZ(0.473 , -0.961 , -2.271);
		op1 = XYZ(2.245 ,  0.278 , -1.032);
		op2 = XYZ(0.218 ,  1.498 , -1.991);

		XYZ tmp = p - o5;
		o5Support = o5 + tmp*(1.0/tmp.length());
		tmp = p - op1;
		op1Support = op1 + tmp*(1.0/tmp.length());
		tmp = p - op2;
		op2Support = op2 + tmp*(1.0/tmp.length());

		e = 0;
	}

	PhophateGroupLocal(const XYZ& p, const XYZ& o5, const XYZ& op1, const XYZ& op2, double e){
		this->p = p;
		this->o5 = o5;
		this->op1 = op1;
		this->op2 = op2;

		XYZ tmp = p - o5;
		o5Support = o5 + tmp*(1.0/tmp.length());
		tmp = p - op1;
		op1Support = op1 + tmp*(1.0/tmp.length());
		tmp = p - op2;
		op2Support = op2 + tmp*(1.0/tmp.length());

		this->e = e;
	}

	PhophateGroupLocal(const string& s){
		vector<string> spt;
		splitString(s, " ", &spt);
		this->p = XYZ(atof(spt[0].c_str()), atof(spt[1].c_str()), atof(spt[2].c_str()));
		this->o5 = XYZ(atof(spt[3].c_str()), atof(spt[4].c_str()), atof(spt[5].c_str()));
		this->op1 = XYZ(atof(spt[6].c_str()), atof(spt[7].c_str()), atof(spt[8].c_str()));
		this->op2 = XYZ(atof(spt[9].c_str()), atof(spt[10].c_str()), atof(spt[11].c_str()));
		XYZ tmp = p - o5;
		o5Support = o5 + tmp*(1.0/tmp.length());
		tmp = p - op1;
		op1Support = op1 + tmp*(1.0/tmp.length());
		tmp = p - op2;
		op2Support = op2 + tmp*(1.0/tmp.length());
		this->e = atof(spt[12].c_str());
	}

	PhophateGroupLocal& operator=(const PhophateGroupLocal& other){
		this->p = other.p;
		this->o5 = other.o5;
		this->op1 = other.op1;
		this->op2 = other.op2;
		this->o5Support = other.o5Support;
		this->op1Support = other.op1Support;
		this->op2Support = other.op2Support;
		this->e = other.e;
		return *this;
	}

	bool equalTo(PhophateGroupLocal& other){
		if(squareDistance(p, other.p) > 0.000001) return false;
		if(squareDistance(o5, other.o5) >  0.000001) return false;
		if(squareDistance(op1, other.op1) >  0.000001) return false;
		if(squareDistance(op2, other.op2) >  0.000001) return false;
		if(squareDistance(o5Support, other.o5Support) >  0.000001) return false;
		if(squareDistance(op1Support, other.op1Support) >  0.000001) return false;
		if(squareDistance(op2Support, other.op2Support) >  0.000001) return false;

		return true;
	}

	string toString(){
		char xx[200];
		sprintf(xx, "%s%s%s%s %6.3f", p.toString().c_str(), o5.toString().c_str(), op1.toString().c_str(), op2.toString().c_str(), e);
		return string(xx);
	}

	virtual ~PhophateGroupLocal();
};

class PhophateGroup {
public:
	/*
	 * atom order
	 * P       0
	 * O5'     1
	 * OP1     2
	 * OP2     3
	 * O5Supp  4
	 * OP1Supp 5
	 * OP2Supp 6
	 */
	XYZ tList[7];

	double e;

	PhophateGroup() {
		tList[0] = XYZ(0.805 ,  0.299 , -1.356);
		tList[1] = XYZ(0.473 , -0.961 , -2.271);
		tList[2] = XYZ(2.245 ,  0.278 , -1.032);
		tList[3] = XYZ(0.218 ,  1.498 , -1.991);
		XYZ tmp = tList[0] - tList[1];
		tList[4] = tList[1] + tmp*(1.0/tmp.length());
		tmp = tList[0] - tList[2];
		tList[5] = tList[2] + tmp*(1.0/tmp.length());
		tmp = tList[0] - tList[3];
		tList[6] = tList[3] + tmp*(1.0/tmp.length());
		e = 0;
	}

	PhophateGroup(const LocalFrame& cs){
		PhophateGroupLocal local;
		tList[0] = local2global(cs, local.p);
		tList[1] = local2global(cs, local.o5);
		tList[2] = local2global(cs, local.op1);
		tList[3] = local2global(cs, local.op2);
		tList[4] = local2global(cs, local.o5Support);
		tList[5] = local2global(cs, local.op1Support);
		tList[6] = local2global(cs, local.op2Support);
		this->e = 0;
	}

	PhophateGroup(const PhophateGroupLocal& local, const LocalFrame& cs){
		tList[0] = local2global(cs, local.p);
		tList[1] = local2global(cs, local.o5);
		tList[2] = local2global(cs, local.op1);
		tList[3] = local2global(cs, local.op2);
		tList[4] = local2global(cs, local.o5Support);
		tList[5] = local2global(cs, local.op1Support);
		tList[6] = local2global(cs, local.op2Support);
		this->e = local.e;
	}

	PhophateGroup& operator=(const PhophateGroup& other){
		tList[0] = other.tList[0];
		tList[1] = other.tList[1];
		tList[2] = other.tList[2];
		tList[3] = other.tList[3];
		tList[4] = other.tList[4];
		tList[5] = other.tList[5];
		tList[6] = other.tList[6];
		this->e = other.e;
		return *this;
	}

	bool equalTo(PhophateGroup& other){
		for(int i=0;i<7;i++){
			if(squareDistance(this->tList[i], other.tList[i]) > 0.0000001) return false;
		}
		return true;
	}

	virtual ~PhophateGroup();
};

} /* namespace NSPmodel */

#endif /* MODEL_PHOPHATEGROUP_H_ */
