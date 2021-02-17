/*
 * MoveMutator.h
 *
 *  Created on: Feb 7, 2019
 *      Author: s2982206
 */

#ifndef PRED_MOVEMUTATOR_H_
#define PRED_MOVEMUTATOR_H_
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "model/BaseDistanceMatrix.h"
#include "geometry/CsMove.h"
#include "geometry/localframe.h"
#include "tools/StringTool.h"
#include "dataio/datapaths.h"
#include "geometry/TransMatrix.h"

namespace NSPpred {

using namespace NSPgeometry;
using namespace NSPmodel;

class MoveMutator {
public:
	vector<CsMove*> moveList;
	vector<CsMove*> disturbanceList;
	double proportion;
	int moveNum;
	int dbNum;
	MoveMutator() {
		srand(time(NULL));
		this->moveNum = 0;
		this->dbNum = 0;
		this->proportion = 0.0;
	}

	MoveMutator(const string& type, int typeA, int typeB, double theta=2.0);

	MoveMutator(const string& type, CsMove* wildTypeMove, int typeA, int typeB, double cutoff);

	CsMove* randomMove();
	CsMove* randomDisturbance();

	virtual ~MoveMutator();
};


} /* namespace NSPpred */

#endif /* PRED_MOVEMUTATOR_H_ */
