/*
 * EnergyTable.cpp
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#include "forcefield/EnergyTable.h"

namespace NSPforcefield {

	EnergyTable::~EnergyTable(){
		delete this->atET;
		delete this->bpET;
		delete this->pb;
		delete this->roET;

	}

} /* namespace NSPtest */
