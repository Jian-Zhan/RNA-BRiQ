/*
 * EnergyTable2.h
 *
 *  Created on: Sep 30, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_ENERGYTABLE2_H_
#define FORCEFIELD_ENERGYTABLE2_H_

#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"

#include "forcefield/BasePair6DEnergyTable.h"
#include "forcefield/AtomicEnergyTable.h"
#include "forcefield/RiboConnectHashMap.h"
#include "forcefield/RiboConnectToPO3.h"
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


class EnergyTable2{
public:

	AtomicEnergyTable* atET;
	RiboseOxygenEnergyTable* roET;
	NSPpara::Parameter para;


	EnergyTable2(){
		atET = new AtomicEnergyTable(&para);
		roET = new RiboseOxygenEnergyTable();
	}

	EnergyTable2(const string& paraFile){
		para = NSPpara::Parameter(paraFile);
		atET = new AtomicEnergyTable(&para);
		roET = new RiboseOxygenEnergyTable();
	}


	virtual ~EnergyTable2();
};

} /* namespace NSPforceField */

#endif /* FORCEFIELD_ENERGYTABLE2_H_ */
