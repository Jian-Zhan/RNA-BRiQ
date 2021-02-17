/*
 * BaseMoveLibrary.h
 *
 *  Created on: Jul 24, 2019
 *      Author: s2982206
 */

#ifndef PRED_BASEMOVELIBRARY_H_
#define PRED_BASEMOVELIBRARY_H_

#include "geometry/CsMove.h"
#include "geometry/TransMatrix.h"
#include "dataio/datapaths.h"
#include <vector>
#include <fstream>
#include <time.h>

namespace NSPpred {

using namespace NSPgeometry;

class BaseMoveLibrary {
public:

	vector<CsMove> randMoveList1;
	vector<CsMove> randMoveList2;
	int randMoveNum;

	BaseMoveLibrary();

	CsMove getRandomMove1(){
		return randMoveList1[rand()%randMoveNum];
	}

	CsMove getRandomMove1(const CsMove& oldMove){
		return oldMove + randMoveList1[rand()%randMoveNum];
	}

	CsMove getRandomMove2(){
		return randMoveList2[rand()%randMoveNum];
	}

	CsMove getRandomMove2(const CsMove& oldMove){
		return oldMove + randMoveList2[rand()%randMoveNum];
	}

	virtual ~BaseMoveLibrary();
};

} /* namespace NSPpred */

#endif /* PRED_BASEMOVELIBRARY_H_ */
