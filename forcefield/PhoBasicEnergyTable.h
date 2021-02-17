/*
 * PhoBasicEnergyTable.h
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_PHOBASICENERGYTABLE_H_
#define FORCEFIELD_PHOBASICENERGYTABLE_H_
#include "geometry/xyz.h"
#include "dataio/datapaths.h"
#include <fstream>

namespace NSPforcefield {
using namespace NSPgeometry;
using namespace std;
using namespace NSPdataio;

class basePhoET0 {
	/*
	 * x from -3.8 to 2.0
	 * y from -6.0 to 5.4
	 * z from -4.8 to 5.0
	 */
public:
	double* bpET;

	basePhoET0();
	double getEnergy(const XYZ& t);
	virtual ~basePhoET0();
};

class basePhoET1 {
	/*
	 * x from -4.8 to 0.0
	 * y from -4.6 to 3.8
	 * z from -5.0 to 3.5
	 */
public:
	double* bpET;

	basePhoET1();
	double getEnergy(const XYZ& t);
	virtual ~basePhoET1();
};

class phoPhoET{

public:

	double* pp1;
	phoPhoET();
	double getPhoPhoEnergy(double d, int sep);
	virtual ~phoPhoET();
};


} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PHOBASICENERGYTABLE_H_ */
