/*
 * BaseStackingEnergyTable.h
 *
 *  Created on: Aug 16, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_BASESTACKINGENERGYTABLE_H_
#define FORCEFIELD_BASESTACKINGENERGYTABLE_H_
#include "geometry/xyz.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"
#include "tools/StringTool.h"
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>

namespace NSPforcefield {
using namespace std;
using namespace NSPgeometry;

class BaseStackingEnergyTable {
public:
	/*
	 * 0: A1
	 * 1: A2
	 * 2: U
	 * 3: G1
	 * 4: G2
	 * 5: C
	 */

	vector<double> etList[24];
	double wt;

	BaseStackingEnergyTable();

	void setWeight(double wt){
		this->wt = wt;
	}

	double getEnergy(int baseType, int centerType, const XYZ& localCoord){
		if(localCoord.x_ <= -2 || localCoord.x_ >= 8) return 0.0;
		if(localCoord.y_ <= -4 || localCoord.y_ >=6) return 0.0;
		if(localCoord.z_ <= -5 || localCoord.z_ >= 5) return 0.0;
		int idX = (int)((localCoord.x_+2.0)*5);
		int idY = (int)((localCoord.y_+4.0)*5);
		int idZ = (int)((localCoord.z_+5.0)*5);
		return etList[baseType*6+centerType][idX*2500+idY*50+idZ] * wt;
	}

	virtual ~BaseStackingEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_BASESTACKINGENERGYTABLE_H_ */
