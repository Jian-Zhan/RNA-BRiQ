/*
 * BasePairEnergyTable.h
 *
 *  Created on: Jul 12, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_BASEPAIRENERGYTABLE_H_
#define FORCEFIELD_BASEPAIRENERGYTABLE_H_

#include "model/BasePairComposition.h"
#include "model/BaseDistanceMatrix.h"
#include "model/PhoDistanceMatrix.h"
#include "math/KDTree.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"
#include "geometry/xyz.h"
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;
using namespace NSPmath;
using namespace NSPgeometry;


class BasePairEnergyTable {


public:
	map<string,double> nbMapList[16];
	map<string,XYZ> nbKeyToP[16];
	map<string,double> nnbMapList[16];
	map<string,double>::iterator it;
	map<string,XYZ>::iterator it2;

	BasePairEnergyTable(); //bw1.2

	void setWeight(double wtNb, double wtNnb);

	double getEnergyNb(const string& key, int typeA, int typeB) {

		it =  nbMapList[typeA*4+typeB].find(key);
		if(it == nbMapList[typeA*4+typeB].end())
			return 0;
		else
			return it->second;
	}

	double getEnergyNnb(const string& key, int typeA, int typeB) {
		it =  nnbMapList[typeA*4+typeB].find(key);
		if(it == nnbMapList[typeA*4+typeB].end())
			return 4.0;
		else
			return it->second;
	}

	XYZ getP(const string& key, int typeA, int typeB) {
		it2 =  nbKeyToP[typeA*4+typeB].find(key);
		if(it2 == nbKeyToP[typeA*4+typeB].end())
			return XYZ();
		else
			return it2->second;
	}

	virtual ~BasePairEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_BASEPAIRENERGYTABLE_H_ */
