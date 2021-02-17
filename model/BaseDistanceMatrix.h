/*
 * BaseDistanceMatrix.h
 *
 *  Created on: Nov 5, 2018
 *      Author: s2982206
 */

#ifndef MODEL_BASEDISTANCEMATRIX_H_
#define MODEL_BASEDISTANCEMATRIX_H_

#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/ProteinRep.h"
#include "tools/StringTool.h"
#include "geometry/TransMatrix.h"
#include <string>


namespace NSPmodel {

using namespace NSPgeometry;
using namespace std;

class BaseDistanceMatrix {
public:
	double dm[16];
	string hashKey;

	BaseDistanceMatrix();
	BaseDistanceMatrix(const string& line);
	BaseDistanceMatrix(LocalFrame& csA, LocalFrame& csB);
	BaseDistanceMatrix(CsMove& move);
	BaseDistanceMatrix(RNABase& baseA, RNABase& baseB);
	BaseDistanceMatrix(vector<double>& dList);
	BaseDistanceMatrix(double dm[]);
	BaseDistanceMatrix& operator=(const BaseDistanceMatrix& other){
		for(int i=0;i<16;i++){
			this->dm[i] = other.dm[i];
		}
		return *this;
	}
	double distanceTo(BaseDistanceMatrix& other);
	double squareDistance(BaseDistanceMatrix& other);
	string riboConnectHashKey(){
		char ss[17];
		ss[16] = '\0';
		double d;
		d = dm[0];
		if(d < 5)
			ss[0] = 'A' + (int)(d*5.0);
		else
			ss[0] = 'Z';

		for(int i=1;i<16;i++){
			d = dm[i];
			if(d < 10.0)
				ss[i] = 'A' + (int)(d*2.5);
			else
				ss[i] = 'Z';
		}
		return string(ss);
	}

	string toString() {
		char a[100];
		sprintf(a, "%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f", dm[0],dm[1],dm[2],dm[3],dm[4],dm[5],dm[6],dm[7],dm[8],dm[9],dm[10],dm[11],dm[12],dm[13],dm[14],dm[15]);
		return string(a);
	}

	virtual ~BaseDistanceMatrix();
};

} /* namespace NSPtest */

#endif /* MODEL_BASEDISTANCEMATRIX_H_ */
