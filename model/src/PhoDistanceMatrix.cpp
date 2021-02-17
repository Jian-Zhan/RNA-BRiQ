/*
 * PhoDistanceMatrix.cpp
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#include "model/PhoDistanceMatrix.h"

namespace NSPmodel {

PhoDistanceMatrix::PhoDistanceMatrix() {
	for(int i=0;i<3;i++){
		dm[i] = 0;
	}
}


PhoDistanceMatrix::~PhoDistanceMatrix() {
}

} /* namespace NSPtest */
