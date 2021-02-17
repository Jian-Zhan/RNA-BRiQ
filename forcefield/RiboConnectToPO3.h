/*
 * RiboConnectToPO3.h
 *
 *  Created on: Aug 5, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_RIBOCONNECTTOPO3_H_
#define FORCEFIELD_RIBOCONNECTTOPO3_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "model/PhophateGroup.h"
#include "geometry/Angles.h"
#include <fstream>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include "para/Parameter.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;
using namespace NSPpara;

class RiboConnectToPO3 {

private:
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

	vector<XYZ> d1d2Lib1; //dihed1-dihed2 library level 1
	vector<vector<XYZ>> d1d2Lib2; //dihed1-dihed2 library level 2
	vector<double> lib2Error;

	vector<double> dihedList1;
	vector<double> dihedList2;


	XYZ c2,o2,c3,o3,p,op1,op2,o5,c5,c4,o4;

	Parameter* para;

public:

	RiboConnectToPO3();
	RiboConnectToPO3(Parameter* para);
	PhophateGroupLocal getPO3Local(XYZ* riboCoordsA, XYZ* riboCoordsB, double impA, double impB);
	PhophateGroupLocal getPO3Local(XYZ* riboCoordsA, XYZ* riboCoordsB, double impA, double impB, bool verbose);

	double printEnergyDetail(XYZ* riboCoordsA, XYZ* riboCoordsB, PhophateGroup& pho, double impA, double impB);

	virtual ~RiboConnectToPO3();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_RIBOCONNECTTOPO3_H_ */
