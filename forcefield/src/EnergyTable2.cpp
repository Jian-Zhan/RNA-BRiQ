/*
 * EnergyTable2.cpp
 *
 *  Created on: Sep 30, 2019
 *      Author: s2982206
 */

#include <forcefield/EnergyTable2.h>

namespace NSPforcefield {



EnergyTable2::~EnergyTable2() {
	// TODO Auto-generated destructor stub
	delete this->atET;
	delete this->roET;
}

} /* namespace NSPforcefield */
