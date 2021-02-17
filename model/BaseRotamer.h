/*
 * BaseRotamer.h
 *
 *  Created on: Jul 17, 2019
 *      Author: s2982206
 */

#ifndef MODEL_BASEROTAMER_H_
#define MODEL_BASEROTAMER_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "model/ProteinRep.h"
#include <vector>
#include "tools/StringTool.h"


namespace NSPmodel {


using namespace NSPgeometry;

class BaseRotamer {

	/*
	 * Atom Order:
	 * 0:  C1'
	 * 1:  C2'
	 * 2:  C3'
	 * 3:  C4'
	 * 4:  O4'
	 * 5:  O2'
	 * 6:  O3'
	 * 7:  C5'
	 * 8:  O2' support
	 * 9:  O3' support
	 * 10: O4' support
	 */
public:

	int resType;
	int rotTypeLv1;
	int rotType;

	XYZ tList1[11]; //local coordinate in cs1
	XYZ center;
	CsMove mv12;
	CsMove mv13;
	CsMove mv31;
	CsMove mv21;
	double chi;
	double improper;
	double energy;

	BaseRotamer(){
		this->resType = 0;
		this->rotTypeLv1 = 0;
		this->rotType = 0;
		this->energy = 0;
		this->improper = 0;
		this->chi = 0;
	}

	BaseRotamer(RNABase* base);
	BaseRotamer(const string& line);


	BaseRotamer& operator=(const BaseRotamer& other){
		this->resType = other.resType;
		this->rotTypeLv1 = other.rotTypeLv1;
		this->rotType = other.rotType;
		for(int i=0;i<8;i++){
			tList1[i] = other.tList1[i];
		}
		center = other.center;
		mv12 = other.mv12;
		mv13 = other.mv13;
		mv31 = other.mv31;
		mv21 = other.mv21;
		energy = other.energy;
		return *this;
	}

	bool equalTo(const BaseRotamer& other){
		if(this->resType != other.resType) return false;
		if(this->rotTypeLv1 != other.rotTypeLv1) return false;
		if(this->rotType != other.rotType) return false;
		for(int i=0;i<11;i++){
			if(squareDistance(tList1[i], other.tList1[i]) > 0.0001) return false;
		}
		if(squareDistance(center, other.center) > 0.0001) return false;
		return true;
	}

	double distanceTo(BaseRotamer* other){
		double d = 0;
		for(int i=0;i<8;i++){
			d += tList1[i].squaredDistance(other->tList1[i]);
		}
		return sqrt(d*0.125);
	}

	void checkRotamer(){
		LocalFrame cs1;
		LocalFrame cs2 = LocalFrame(tList1[1], tList1[2], tList1[6]);
		LocalFrame cs3 = LocalFrame(tList1[4], tList1[3], tList1[7]);

		LocalFrame csX2 = cs1 + mv12;
		if(!cs2.equalTo(csX2)){
			cout << "mv12 error" << endl;
			cs2.print();
			cout << endl;
			csX2.print();
			cout << endl;
		}
		LocalFrame csX3 = cs1 + mv13;
		if(!cs3.equalTo(csX3)){
			cout << "mv13 error" << endl;
			cs3.print();
			cout << endl;
			csX3.print();
			cout << endl;
		}
		LocalFrame csX1 = cs3 + mv31;
		if(!cs1.equalTo(csX1)) {
			cout << "mv31 error" << endl;
			cs1.print();
			cout << endl;
			csX1.print();
			cout << endl;
		}

		LocalFrame csY1 = cs2 + mv21;
		if(!cs1.equalTo(csY1)) {
			cout << "mv21 error" << endl;
			cs1.print();
			cout << endl;
			csY1.print();
			cout << endl;
		}

		XYZ tmp = tList1[1] - tList1[5];
		XYZ t8 = tList1[5] + tmp * (1.0/tmp.length());
		tmp = tList1[2] - tList1[6];
		XYZ t9 = tList1[6] + tmp * (1.0/tmp.length());
		tmp = tList1[0] - tList1[4] + tList1[3] - tList1[4];
		XYZ t10 = tList1[4] + tmp * (1.0/tmp.length());

		if(t8.distance(tList1[8]) > 0.01)
			cout << "O2' support error" << endl;
		if(t9.distance(tList1[9]) > 0.01)
			cout << "O3' support error" << endl;
		if(t10.distance(tList1[10]) > 0.01)
			cout << "O4' support error" << endl;


	}
	virtual ~BaseRotamer();
};

} /* namespace NSPforcefield */

#endif /* MODEL_BASEROTAMER_H_ */
