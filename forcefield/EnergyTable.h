/*
 * EnergyTable.h
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#ifndef FORCEFIELD_ENERGYTABLE_H_
#define FORCEFIELD_ENERGYTABLE_H_

#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"

#include "forcefield/BasePair6DEnergyTable.h"
#include "forcefield/AtomicEnergyTable.h"
#include "forcefield/RiboConnectHashMap.h"
#include "forcefield/RiboConnectToPO3.h"
#include "forcefield/PO3Builder.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "forcefield/SimpleRiboConnect.h"
#include "forcefield/BaseStackingEnergyTable.h"

#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;


class EnergyTable{
public:
	BasePair6DEnergyTable* bpET;
	AtomicEnergyTable* atET;
	PO3Builder* pb;
	RiboseOxygenEnergyTable* roET;
	NSPpara::Parameter para;

	EnergyTable(){

		bpET = new BasePair6DEnergyTable(&para);
		atET = new AtomicEnergyTable(&para);
		roET = new RiboseOxygenEnergyTable();
		pb = new PO3Builder(&para,atET,roET);
	}

	EnergyTable(const string& paraFile){
		para = NSPpara::Parameter(paraFile);
		bpET = new BasePair6DEnergyTable(&para);
		atET = new AtomicEnergyTable(&para);
		roET = new RiboseOxygenEnergyTable();
		pb = new PO3Builder(&para,atET,roET);
	}

	virtual ~EnergyTable();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_ENERGYTABLE_H_ */
