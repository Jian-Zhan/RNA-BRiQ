/*
 * ThreeBaseMoveLibrary.h
 *
 *  Created on: Sep 13, 2019
 *      Author: s2982206
 */

#ifndef PRED_THREEBASEMOVELIBRARY_H_
#define PRED_THREEBASEMOVELIBRARY_H_

#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "dataio/datapaths.h"
#include <vector>
#include <fstream>
#include <time.h>

namespace NSPpred {

using namespace NSPgeometry;

class ThreeBaseMoveLibrary {
public:

	int rotNum1; //baseA contact with baseC
	int rotNum2; //all f3 fragments

	vector<CsMove> mvLib1AB;
	vector<CsMove> mvLib1BC;

	vector<CsMove> mvLib2AB;
	vector<CsMove> mvLib2BC;

	ThreeBaseMoveLibrary(int typeA, int typeB, int typeC);

	pair<CsMove,CsMove> getRandomMoveLib1(){
		int randIndex = rand()%rotNum1;
		pair<CsMove, CsMove> p(mvLib1AB[randIndex], mvLib1BC[randIndex]);
		return p;
	}

	pair<CsMove,CsMove> getRandomMoveLib2(){
		int randIndex = rand()%rotNum2;
		pair<CsMove, CsMove> p(mvLib2AB[randIndex], mvLib2BC[randIndex]);
		return p;
	}


	virtual ~ThreeBaseMoveLibrary();
};

} /* namespace NSPforcefield */

#endif /* PRED_THREEBASEMOVELIBRARY_H_ */
