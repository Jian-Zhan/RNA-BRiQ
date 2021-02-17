/*
 * RiboseOxygenEnergyTable.h
 *
 *  Created on: Jul 26, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_
#define FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_

#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"
#include "tools/StringTool.h"
#include <time.h>
#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPpara;

class RiboseOxygenEnergyTable {
public:

	/*
	 * 0: O2'
	 * 1: O3'
	 * 2: O4'
	 * 3: O5'
	 * 4: OP
	 */

	vector<double> etList[20];
	vector<bool> directDep[20];
	vector<XYZ> hbAtom[20];
	vector<double> ang1List[20];
	vector<double> ang2List[20];
	vector<double> ang3List[20];

	/*
	 * sep=1
	 * base-OP
	 */
	vector<double> etListP1[4];
	vector<bool> directDepP1[4];
	vector<XYZ> hbAtomP1[4];
	vector<double> ang1ListP1[4];
	vector<double> ang2ListP1[4];
	vector<double> ang3ListP1[4];

	/*
	 * sep=-1
	 * base-O4
	 * base-OP
	 */
	vector<double> etListM1[8];
	vector<double> directDepM1[8];
	vector<XYZ> hbAtomM1[8];
	vector<double> ang1ListM1[8];
	vector<double> ang2ListM1[8];
	vector<double> ang3ListM1[8];

	RiboseOxygenEnergyTable();

	double getEnergy(int baseType, int oxygenType, const XYZ& localCoreCoord, const XYZ& localSupportCoord, int sep){
		if(sep == 0)
			return 0.0;
		else if(sep == 1){
			if(oxygenType == 4){
				if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
				if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
				if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
				int idX = (int)((localCoreCoord.x_+3.0)*5);
				int idY = (int)((localCoreCoord.y_+6.0)*5);
				int idZ = (int)((localCoreCoord.z_+5.0)*5);
				int type = baseType;
				int index = idX*4000+idY*50+idZ;

				double ene = etListP1[type][index];
				if(ene >= 0)
					return ene;
				if(!directDepP1[type][index])
					return ene;


				double ang = NSPgeometry::angleX(hbAtomP1[type][index], localCoreCoord, localSupportCoord);
				double ang1 = ang1ListP1[type][index];
				double ang2 = ang2ListP1[type][index];
				double ang3 = ang3ListP1[type][index];
				if(ang < ang1 || ang > ang3)
					return 0;

				if(ang < ang2){
					double u = (ang - ang2)/(ang1 - ang2);
					return ene * (1 - u*u);
				}
				else {
					double u = (ang - ang2)/(ang2 - ang3);
					return ene * (1 - u*u);
				}
			}
			else
				return 0.0;
		}
		else if(sep == -1){
			if(oxygenType == 2){
				if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
				if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
				if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
				int idX = (int)((localCoreCoord.x_+3.0)*5);
				int idY = (int)((localCoreCoord.y_+6.0)*5);
				int idZ = (int)((localCoreCoord.z_+5.0)*5);
				int type = baseType*2;
				int index = idX*4000+idY*50+idZ;

				double ene = etListM1[type][index];
				if(ene >= 0)
					return ene;
				if(!directDepM1[type][index])
					return ene;

				double ang = NSPgeometry::angleX(hbAtomM1[type][index], localCoreCoord, localSupportCoord);
				double ang1 = ang1ListM1[type][index];
				double ang2 = ang2ListM1[type][index];
				double ang3 = ang3ListM1[type][index];
				if(ang < ang1 || ang > ang3)
					return 0;

				if(ang < ang2){
					double u = (ang - ang2)/(ang1 - ang2);
					return ene * (1 - u*u);
				}
				else {
					double u = (ang - ang2)/(ang2 - ang3);
					return ene * (1 - u*u);
				}
			}
			else if(oxygenType == 4){
				if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
				if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
				if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
				int idX = (int)((localCoreCoord.x_+3.0)*5);
				int idY = (int)((localCoreCoord.y_+6.0)*5);
				int idZ = (int)((localCoreCoord.z_+5.0)*5);
				int type = baseType*2+1;
				int index = idX*4000+idY*50+idZ;

				double ene = etListM1[type][index];
				if(ene >= 0)
					return ene;
				if(!directDepM1[type][index])
					return ene;

				double ang = NSPgeometry::angleX(hbAtomM1[type][index], localCoreCoord, localSupportCoord);
				double ang1 = ang1ListM1[type][index];
				double ang2 = ang2ListM1[type][index];
				double ang3 = ang3ListM1[type][index];
				if(ang < ang1 || ang > ang3)
					return 0;

				if(ang < ang2){
					double u = (ang - ang2)/(ang1 - ang2);
					return ene * (1 - u*u);
				}
				else {
					double u = (ang - ang2)/(ang2 - ang3);
					return ene * (1 - u*u);
				}
			}
			else
				return 0.0;
		}
		else {
			if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
			if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
			if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
			int idX = (int)((localCoreCoord.x_+3.0)*5);
			int idY = (int)((localCoreCoord.y_+6.0)*5);
			int idZ = (int)((localCoreCoord.z_+5.0)*5);
			int type = baseType*5+oxygenType;
			int index = idX*4000+idY*50+idZ;

			double ene = etList[type][index];
			if(ene > 0)
				return ene;
			if(!directDep[type][index])
				return ene;

			double ang = NSPgeometry::angleX(hbAtom[type][index], localCoreCoord, localSupportCoord);
			double ang1 = ang1List[type][index];
			double ang2 = ang2List[type][index];
			double ang3 = ang3List[type][index];
			if(ang < ang1 || ang > ang3)
				return 0;

			if(ang < ang2){
				double u = (ang - ang2)/(ang1 - ang2);
				return ene * (1 - u*u);
			}
			else {
				double u = (ang - ang2)/(ang2 - ang3);
				return ene * (1 - u*u);
			}
		}
	}

	/*
	double getEnergy(int baseType, int oxygenType, const XYZ& localCoreCoord, const XYZ& localSupportCoord){
		if(localCoreCoord.x_ <= -3 || localCoreCoord.x_ >= 12) return 0.0;
		if(localCoreCoord.y_ <= -6 || localCoreCoord.y_ >= 10) return 0.0;
		if(localCoreCoord.z_ <= -5 || localCoreCoord.z_ >= 5) return 0.0;
		int idX = (int)((localCoreCoord.x_+3.0)*5);
		int idY = (int)((localCoreCoord.y_+6.0)*5);
		int idZ = (int)((localCoreCoord.z_+5.0)*5);
		int type = baseType*5+oxygenType;
		int index = idX*4000+idY*50+idZ;

		double ene = etList[type][index];
		if(ene > 0)
			return ene;
		if(!directDep[type][index])
			return ene;

		double ang = NSPgeometry::angleX(hbAtom[type][index], localCoreCoord, localSupportCoord);
		double ang1 = ang1List[type][index];
		double ang2 = ang2List[type][index];
		double ang3 = ang3List[type][index];
		if(ang < ang1 || ang > ang3)
			return 0;

		if(ang < ang2){
			double u = (ang - ang2)/(ang1 - ang2);
			return ene * (1 - u*u);
		}
		else {
			double u = (ang - ang2)/(ang2 - ang3);
			return ene * (1 - u*u);
		}
	}
	*/

	/*
	double getEnergy(int baseType, int oxygenType, const XYZ& localCoord){
		if(localCoord.x_ <= -3 || localCoord.x_ >= 12) return 0.0;
		if(localCoord.y_ <= -6 || localCoord.y_ >= 10) return 0.0;
		if(localCoord.z_ <= -5 || localCoord.z_ >= 5) return 0.0;
		int idX = (int)((localCoord.x_+3.0)*5);
		int idY = (int)((localCoord.y_+6.0)*5);
		int idZ = (int)((localCoord.z_+5.0)*5);
		return etList[baseType*5+oxygenType][idX*4000+idY*50+idZ];
	}
	*/
	virtual ~RiboseOxygenEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_RIBOSEOXYGENENERGYTABLE_H_ */
