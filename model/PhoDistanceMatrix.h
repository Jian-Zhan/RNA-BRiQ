/*
 * PhoDistanceMatrix.h
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#ifndef MODEL_PHODISTANCEMATRIX_H_
#define MODEL_PHODISTANCEMATRIX_H_
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "model/ProteinRep.h"
#include <string>

namespace NSPmodel {

class PhoDistanceMatrix {
public:
	double dm[3];
	PhoDistanceMatrix();
	PhoDistanceMatrix(vector<double> dList){
		for(int i=0;i<3;i++){
			dm[i] = dList[i];
		}
	}

	PhoDistanceMatrix(double dList[]){
		for(int i=0;i<3;i++){
			dm[i] = dList[i];
		}
	}

	PhoDistanceMatrix(const XYZ& p){
		dm[0] = p.x_;
		dm[1] = p.y_;
		dm[2] = p.z_;
	}

	double distanceTo(const PhoDistanceMatrix& other){
		double a = this->dm[0]- other.dm[0];
		double b = this->dm[1]- other.dm[1];
		double c = this->dm[2]- other.dm[2];
		return sqrt(a*a+b*b+c*c);
	}



	virtual ~PhoDistanceMatrix();
};

} /* namespace NSPtest */

#endif /* MODEL_PHODISTANCEMATRIX_H_ */
