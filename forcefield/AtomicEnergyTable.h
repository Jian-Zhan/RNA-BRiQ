/*
 * AtomicEnergyTable.h
 *
 *  Created on: Jul 27, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_ATOMICENERGYTABLE_H_
#define FORCEFIELD_ATOMICENERGYTABLE_H_

#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include "geometry/localframe.h"
#include "dataio/datapaths.h"
#include "para/Parameter.h"
#include "tools/StringTool.h"


namespace NSPforcefield {

using namespace std;
using namespace NSPpara;
using namespace NSPgeometry;
using namespace NSPtools;

class AtomicEnergyTable {
public:
	vector<vector<double>> et1;  //atomic energy sep 1
	vector<vector<double>> et2;  //atomic energy sep -1
	vector<vector<double>> et3;  //atomic energy sep 2

	vector<vector<double>> vdwEt;

	vector<double> O2O2HbondEnergy; //o2-o2 hbond energy
	vector<double> O2OPHbondEnergy; //o2-op hbond energy
	vector<double> OPOPHbondEnergy; //op-op hbond energy

	vector<double> baseAtomClashDistance[88];
	vector<double> baseAtomNearestAtomIndex[88];

	map<int, int> baseAtomIndexToClashAtomIndex;

	Parameter* para;

	AtomicEnergyTable(Parameter* para);

	double vdwEnergy(double d0, double d){
		double u, e, slop, e1;

		if(d < d0 - 0.4){
			u = -0.4*para->lamdaClash;
			slop = 4*u*u*u*para->lamdaClash;
			e1 = u*u*u*u;
			e = e1 + slop*(d-d0+0.4);
		}
		else if(d < d0) {
			u = (d-d0)*para->lamdaClash;
			e = u*u*u*u;
		}
		else
			e = 0;
		return e*para->wtClash;
	}

	double getBaseBaseClashEnergy(int baseTypeA, LocalFrame& csA, XYZ* coordsA, int baseTypeB, LocalFrame& csB, XYZ* coordsB, double shift){
		if(squareDistance(csA.origin_, csB.origin_) > 256.0) return 0.0;
		int i,indexA, indexB, idX, idY, idZ;
		double e1, e2, w1, w2;
		int nearestAtomIndex;
		double dd, clashDistance;
		XYZ t;
		double ene = 0.0;
		//baseA atom in csB
		for(i=1;i<10;i++){
			indexA = baseAtomIndexToClashAtomIndex[baseTypeA*100+i];
			if(indexA < 0) continue;

			t = global2local(csB, coordsA[i]);
			if(t.x_ <= -3 || t.x_ >= 12) continue;
			if(t.y_ <= -6 || t.y_ >= 10) continue;
			if(t.z_ <= -5 || t.z_ >= 5) continue;
			idX = (int)((t.x_+3.0)*5);
			idY = (int)((t.y_+6.0)*5);
			idZ = (int)((t.z_+5.0)*5);
			nearestAtomIndex = this->baseAtomNearestAtomIndex[baseTypeB*22 + indexA][idX*4000+idY*50+idZ];
			clashDistance = this->baseAtomClashDistance[baseTypeB*22 + indexA][idX*4000+idY*50+idZ] - shift;
			if(clashDistance < 0) continue;
			dd = squareDistance(coordsA[i], coordsB[nearestAtomIndex]);
			if(dd >= clashDistance*clashDistance) continue;
			ene += 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
		}

		//baseB atom in csA
		for(i=1;i<10;i++){
			indexB = baseAtomIndexToClashAtomIndex[baseTypeB*100+i];
			if(indexB < 0) continue;

			t = global2local(csA, coordsB[i]);
			if(t.x_ <= -3 || t.x_ >= 12) continue;
			if(t.y_ <= -6 || t.y_ >= 10) continue;
			if(t.z_ <= -5 || t.z_ >= 5) continue;
			idX = (int)((t.x_+3.0)*5);
			idY = (int)((t.y_+6.0)*5);
			idZ = (int)((t.z_+5.0)*5);
			nearestAtomIndex = this->baseAtomNearestAtomIndex[baseTypeA*22 + indexB][idX*4000+idY*50+idZ];
			clashDistance = this->baseAtomClashDistance[baseTypeA*22 + indexB][idX*4000+idY*50+idZ] - shift;
			if(clashDistance < 0) continue;
			dd = squareDistance(coordsB[i], coordsA[nearestAtomIndex]);
			if(dd >= clashDistance*clashDistance) continue;
			ene += 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
		}
		return ene;
	}


	double getBaseBaseClashEnergyVerbose(int baseTypeA, LocalFrame& csA, XYZ* coordsA, int baseTypeB, LocalFrame& csB, XYZ* coordsB, double shift){
		if(squareDistance(csA.origin_, csB.origin_) > 256.0) return 0.0;
		int i,indexA, indexB, idX, idY, idZ;
		double e1, e2, w1, w2;
		int nearestAtomIndex;
		double dd, clashDistance;
		XYZ t;
		double ene = 0.0;
		//baseA atom in csB
		for(i=1;i<10;i++){
			indexA = baseAtomIndexToClashAtomIndex[baseTypeA*100+i];
			if(indexA < 0) continue;

			t = global2local(csB, coordsA[i]);
			if(t.x_ <= -3 || t.x_ >= 12) continue;
			if(t.y_ <= -6 || t.y_ >= 10) continue;
			if(t.z_ <= -5 || t.z_ >= 5) continue;
			idX = (int)((t.x_+3.0)*5);
			idY = (int)((t.y_+6.0)*5);
			idZ = (int)((t.z_+5.0)*5);
			nearestAtomIndex = this->baseAtomNearestAtomIndex[baseTypeB*22 + indexA][idX*4000+idY*50+idZ];
			clashDistance = this->baseAtomClashDistance[baseTypeB*22 + indexA][idX*4000+idY*50+idZ] - shift;
			if(clashDistance < 0) continue;
			dd = squareDistance(coordsA[i], coordsB[nearestAtomIndex]);
			if(dd >= clashDistance*clashDistance) continue;
			ene += 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
			double e = 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
			if(e > 0.1){
				printf("%5.3f %d %5.3f %6.3f\n", clashDistance, nearestAtomIndex, sqrt(dd), e);
			}
		}


		//baseB atom in csA
		for(i=1;i<10;i++){
			indexB = baseAtomIndexToClashAtomIndex[baseTypeB*100+i];
			if(indexB < 0) continue;

			t = global2local(csA, coordsB[i]);
			if(t.x_ <= -3 || t.x_ >= 12) continue;
			if(t.y_ <= -6 || t.y_ >= 10) continue;
			if(t.z_ <= -5 || t.z_ >= 5) continue;
			idX = (int)((t.x_+3.0)*5);
			idY = (int)((t.y_+6.0)*5);
			idZ = (int)((t.z_+5.0)*5);
			nearestAtomIndex = this->baseAtomNearestAtomIndex[baseTypeA*22 + indexB][idX*4000+idY*50+idZ];
			clashDistance = this->baseAtomClashDistance[baseTypeA*22 + indexB][idX*4000+idY*50+idZ] - shift;
			if(clashDistance < 0) continue;
			dd = squareDistance(coordsB[i], coordsA[nearestAtomIndex]);
			if(dd >= clashDistance*clashDistance) continue;

			ene += 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
			double e = 0.5*vdwEt[(int)(clashDistance*100)][(int)(dd*100)];
			if(e > 0.1){
				printf("%5.3f %d %5.3f %6.3f\n", clashDistance, nearestAtomIndex, sqrt(dd), e);
			}
		}
		return ene;
	}

	double getBaseBaseEnergy(int baseTypeA, int atomTypeA, int baseTypeB, int atomTypeB, double dd, int sep){
		if(dd >= 16.0) return 0;
		if(sep == 0) return 0;

		int typeA, typeB;
		if(baseTypeA == 0)
			typeA = atomTypeA + 10;
		else if(baseTypeA == 1)
			typeA = atomTypeA + 20;
		else if(baseTypeA == 2)
			typeA = atomTypeA + 28;
		else
			typeA = atomTypeA + 39;

		if(baseTypeB == 0)
			typeB = atomTypeB + 10;
		else if(baseTypeB == 1)
			typeB = atomTypeB + 20;
		else if(baseTypeB == 2)
			typeB = atomTypeB + 28;
		else
			typeB = atomTypeB + 39;

		int index = typeA*47+typeB;
		if(index >= 2209){
			cout << baseTypeA << " " << atomTypeA << " " << baseTypeB << " " << atomTypeB << endl;
		}

		if(sep == 1)
			return et1[index][(int)(dd*100)];
		else if(sep == -1)
			return et2[index][(int)(dd*100)];
		else
			return et3[index][(int)(dd*100)];
	}

	double getBaseRiboseEnergy(int baseTypeA, int atomTypeA, int riboseAtomTypeB, double dd, int sep){

		if(dd >= 16.00) return 0;
		if(sep == 0) return 0;
		int typeA;
		if(baseTypeA == 0)
			typeA = atomTypeA + 10;
		else if(baseTypeA == 1)
			typeA = atomTypeA + 20;
		else if(baseTypeA == 2)
			typeA = atomTypeA + 28;
		else
			typeA = atomTypeA + 39;

		int index = typeA*47+riboseAtomTypeB;
		if(index >= 2209) {
			cout << "invalid type: " << baseTypeA << " " << atomTypeA << " " << riboseAtomTypeB << " " << endl;
			exit(1);
		}

		if(sep == 1)
			return et1[index][(int)(dd*100)];
		else if(sep == -1)
			return et2[index][(int)(dd*100)];
		else
			return et3[index][(int)(dd*100)];
	}

	double getRiboseRiboseEnergy(int riboseTypeA, int riboseTypeB, double dd, int sep){

		if(dd >= 16.00) return 0;
		if(sep == 0) return 0;
		int index = riboseTypeA*47+riboseTypeB;
		if(index >= 2209) {
			cout << "invalid type: " << riboseTypeA << " " << riboseTypeB << endl;
			exit(1);
		}
		if(sep == 1)
			return et1[index][(int)(dd*100)];
		else if(sep == -1)
			return et2[index][(int)(dd*100)];
		else
			return et3[index][(int)(dd*100)];
	}

	int getIndex(int baseTypeA, int atomTypeA, int phoAtomTypeB, double dd, int sep){

		int typeA;
		if(baseTypeA == 0)
			typeA = atomTypeA + 10;
		else if(baseTypeA == 1)
			typeA = atomTypeA + 20;
		else if(baseTypeA == 2)
			typeA = atomTypeA + 28;
		else
			typeA = atomTypeA + 39;

		int typeB = phoAtomTypeB + 7;
		if(typeB > 9) typeB = 9;

		int index = typeA*47+typeB;
		return index;
	}

	double getBasePhoEnergy(int baseTypeA, int atomTypeA, int phoAtomTypeB, double dd, int sep){

		if(dd >= 16.00) return 0;
		if(sep == 0) return 0;
		int typeA;
		if(baseTypeA == 0)
			typeA = atomTypeA + 10;
		else if(baseTypeA == 1)
			typeA = atomTypeA + 20;
		else if(baseTypeA == 2)
			typeA = atomTypeA + 28;
		else
			typeA = atomTypeA + 39;

		int typeB = phoAtomTypeB + 7;
		if(typeB > 9) typeB = 9;

		int index = typeA*47+typeB;
		if(index >= 2209) {
			cout << "invalid type: " << baseTypeA << " " << atomTypeA << " " << phoAtomTypeB << endl;
			exit(1);
		}
		if(sep == 1)
			return et1[index][(int)(dd*100)];
		else if(sep == -1)
			return et2[index][(int)(dd*100)];
		else
			return et3[index][(int)(dd*100)];
	}

	double getRibosePhoEnergy(int riboseTypeA, int phoTypeB, double dd, int sep){

		if(dd >= 16.00) return 0;
		if(sep == 0) return 0;
		int typeB = phoTypeB + 7;
		if(typeB > 9) typeB = 9;
		int index = riboseTypeA*47+typeB;
		if(index >= 2209) {
			cout << "invalid type: " << riboseTypeA << " " << phoTypeB << endl;
			exit(1);
		}
		if(sep == 1)
			return et1[index][(int)(dd*100)];
		else if(sep == -1)
			return et2[index][(int)(dd*100)];
		else
			return et3[index][(int)(dd*100)];
	}

	double getPhoPhoEnergy(int phoTypeA, int phoTypeB, double dd, int sep){

		double e1 = 0;
		if(dd < 16) {
			if(sep == 0) e1 =0;
			int typeA = phoTypeA + 7;
			if(typeA > 9) typeA = 9;
			int typeB = phoTypeB + 7;
			if(typeB > 9) typeB = 9;
			int index = typeA*36+typeB;
			if(index >= 2209){
				cout << "invalid type: " << phoTypeA << " " << phoTypeB << endl;
				exit(1);
			}
			if(sep == 1)
				e1 = et1[index][(int)(dd*100)];
			else if(sep == -1)
				e1 =  et2[index][(int)(dd*100)];
			else
				e1 = et3[index][(int)(dd*100)];
		}

		return e1;
	}

	double getRiboseRiboseHbondEnergy(XYZ* riboseA, XYZ* riboseB){

		XYZ a = riboseA[1];
		XYZ b = riboseA[5];
		XYZ c = riboseB[5];
		XYZ d = riboseB[1];
		double len = b.distance(c);
		if(len >= 4.0) return 0;
		double ang1 = angleX(a,b,c);
		double ang2 = angleX(b,c,d);
		double dihed = dihedral(a,b,c,d);
		if(dihed < 0) dihed += 360;

		int id1 = (int)(len*10);
		int id2 = (int)(ang1*0.33333333);
		int id3 = (int)(ang2*0.33333333);
		int id4 = (int)(dihed*0.1666666);
		return this->O2O2HbondEnergy[id1*216000+id2*3600+id3*60+id4];
	}

	double getRibosePhoHbondEnergy(XYZ* riboseA, XYZ* phoB){

		XYZ a = riboseA[1];
		XYZ b = riboseA[5];
		XYZ c;
		XYZ d = phoB[0];
		double len;
		double len1 = b.distance(phoB[2]);
		double len2 = b.distance(phoB[3]);
		if(len1 < len2 && len1 < 4.0) {
			c = phoB[2];
			len = len1;
		}
		else if(len1 >= len2 && len2 < 4.0) {
			len = len2;
			c = phoB[3];
		}
		else
			return 0;
		double ang1 = angleX(a,b,c);
		double ang2 = angleX(b,c,d);
		double dihed = dihedral(a,b,c,d);
		if(dihed < 0) dihed += 360;
		int id1 = (int)(len*10);
		int id2 = (int)(ang1*0.33333333);
		int id3 = (int)(ang2*0.33333333);
		int id4 = (int)(dihed*0.1666666);
		return this->O2OPHbondEnergy[id1*216000+id2*3600+id3*60+id4];
	}

	double getPhoPhoHbondEnergy(XYZ* phoA, XYZ* phoB){

		XYZ a = phoA[0];
		XYZ d = phoB[0];
		if(squareDistance(a,d) > 49) return 0;

		double minLen = 99.9;
		double len;
		int idA, idB;
		for(int i=2;i<4;i++){
			for(int j=2;j<4;j++){
				len = phoA[i].distance(phoB[j]);
				if(len < minLen){
					minLen = len;
					idA = i;
					idB = j;
				}
			}
		}

		double hbondEne = 0;
		double rep = 0;

		if(minLen < 4.0){
			XYZ b = phoA[idA];
			XYZ c = phoB[idB];
			double ang1 = angleX(a,b,c);
			double ang2 = angleX(b,c,d);
			double dihed = dihedral(a,b,c,d);
			if(dihed < 0) dihed += 360;
			int id1 = (int)(minLen*10);
			int id2 = (int)(ang1*0.33333333);
			int id3 = (int)(ang2*0.33333333);
			int id4 = (int)(dihed*0.1666666);
			hbondEne = this->OPOPHbondEnergy[id1*216000+id2*3600+id3*60+id4];
		}

		rep = 5.0/minLen;

		return hbondEne + rep;
	}

	virtual ~AtomicEnergyTable();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_ATOMICENERGYTABLE_H_ */
