/*
 * RNALoop.h
 *
 *  Created on: Feb 4, 2019
 *      Author: s2982206
 */

#ifndef PRED_RNALOOP_H_
#define PRED_RNALOOP_H_
#include <pred/BaseNode0.h>
#include "geometry/localframe.h"
#include "geometry/RMSD.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "model/RNABaseName.h"
#include "forcefield/NeighborConnectionGenerator.h"
#include "forcefield/EnergyTable.h"
#include <time.h>
#include <stdlib.h>

namespace NSPpred {
using namespace NSPforceField;
using namespace NSPmodel;
using namespace NSPgeometry;

class RNALoop {
public:
	int length;
	int* typeList;
	BaseNode** nodeList;
	NeighborMove** moveList;
	NeighborConnectionGenerator* ncg;
	EnergyTable* et;

	RNALoop(const string& pdbFile, EnergyTable* et);
	RNALoop(const string& seq, NeighborConnectionGenerator* ncg, EnergyTable* et);
	void sampleLoopWithTerminalConstrain(BaseDistanceMatrix& dm, RNALoop& ref);
	double getTotalEnergy();
	double rmsdTo(RNALoop& other);
	void printLoop(const string& file);
	virtual ~RNALoop();
};

} /* namespace NSPpred */

#endif /* PRED_RNALOOP_H_ */
