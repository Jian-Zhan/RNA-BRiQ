/*
 * NeighborConnectionGenerator.h
 *
 *  Created on: Feb 4, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_NEIGHBORCONNECTIONGENERATOR_H_
#define FORCEFIELD_NEIGHBORCONNECTIONGENERATOR_H_

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "geometry/CsMove.h"
#include "tools/StringTool.h"
#include "dataio/datapaths.h"

namespace NSPforcefield {
using namespace NSPgeometry;

class NeighborMove{
public:
	CsMove move;
	XYZ localP1;
	XYZ localP2;
	double energy;

	NeighborMove() {this->energy = 0.0;}
	NeighborMove(const string& line);
	NeighborMove(CsMove& move, XYZ& p1, XYZ& p2, double ene){
		this->move = move;
		this->localP1 = p1;
		this->localP2 = p2;
		this->energy = ene;
	}

	NeighborMove& operator=(NeighborMove& other) {
		this->move = other.move;
		this->localP1 = other.localP1;
		this->localP2 = other.localP2;
		this->energy = other.energy;
		return *this;
	}
	virtual ~NeighborMove();
};


class NeighborConnectionGenerator {
	vector<vector<NeighborMove*>> moveListList;
	NeighborMove* move0;
public:
	NeighborConnectionGenerator() {this->move0 = NULL;}
	NeighborConnectionGenerator(int type);
	NeighborMove* standardMove(){
		return move0;
	}
	NeighborMove* getRandomMove(int pairType);
	virtual ~NeighborConnectionGenerator();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_NEIGHBORCONNECTIONGENERATOR_H_ */
