/*
 * PO3Builder.h
 *
 *  Created on: Oct 1, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_PO3BUILDER_H_
#define FORCEFIELD_PO3BUILDER_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "forcefield/AtomicEnergyTable.h"
#include "forcefield/RiboseOxygenEnergyTable.h"
#include "model/PhophateGroup.h"
#include "para/Parameter.h"
#include "geometry/Angles.h"
#include <fstream>
#include <vector>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPpara;


class PO3Builder {
public:

	vector<XYZ> pList;
	vector<XYZ> op1List;
	vector<XYZ> op2List;
	vector<XYZ> o5List;
	vector<XYZ> op1SuppList;
	vector<XYZ> op2SuppList;
	vector<XYZ> o5SuppList;

	double len1 = 1.605;
	double len2 = 1.592;
	double len3 = 1.422;
	double ang1 = 120.1;
	double ang2 = 103.5;
	double ang3 = 120.7;
	double ang4 = 111.1;

	vector<double> eImpD1D2;
	vector<double> eImpD4D5;
	vector<double> eD2D4D3;

	vector<XYZ> d1d2Lib1A; //dihed1-dihed2 library level 1, impA < 0
	vector<vector<XYZ>> d1d2Lib2A; //dihed1-dihed2 library level 2 impA < 0

	vector<XYZ> d1d2Lib1B; //dihed1-dihed2 library level 1, impA > 0
	vector<vector<XYZ>> d1d2Lib2B; //dihed1-dihed2 library level 2, impA > 0

	vector<double> lib2ErrorA;
	vector<double> lib2ErrorB;

	vector<double> eDihed1;
	vector<double> eDihed2;
	vector<double> eDihed3;
	vector<double> eDihed4;
	vector<double> eDihed5;

	XYZ c2,c3,o3,p,op1,op2,o5,o5Supp, op1Supp, op2Supp, c5,c4,o4;

	Parameter* para;
	AtomicEnergyTable* aET;
	RiboseOxygenEnergyTable* rET;

	PO3Builder(Parameter* para, AtomicEnergyTable* aET, RiboseOxygenEnergyTable* rET);

	PhophateGroupLocal getPhoLocal(const LocalFrame& cs2A, XYZ* riboCoordsB, double impA, double impB);

	virtual ~PO3Builder();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_PO3BUILDER_H_ */
