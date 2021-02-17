/*
 * SimpleRiboConnect.h
 *
 *  Created on: Aug 7, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_SIMPLERIBOCONNECT_H_
#define FORCEFIELD_SIMPLERIBOCONNECT_H_

#include <vector>
#include "geometry/localframe.h"
#include "geometry/xyz.h"
#include "geometry/Angles.h"
#include <fstream>
#include <map>
#include "forcefield/PhoBasicEnergyTable.h"
#include "geometry/CsMove.h"
#include "pred/PhoBasic.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"

namespace NSPforcefield {
using namespace std;
using namespace NSPgeometry;
using namespace NSPpred;



class SimpleRiboConnect {
public:

	/*
	 * connection geometry energy table
	 * 40*36*36*36 table
	 * d= i1*0.3 C1'-C1'
	 * ang1 = i2*5.0 N-C1'-C1'
	 * ang2 = i3*5.0 C1'-C1'-N
	 * dihed = i4*10.0 N-C1'-C1'-N
	 */

	vector<double> eneList;

	basePhoET0 et0;
	basePhoET1 et1;

	vector<XYZ> phoLib1;
	vector<vector<XYZ>> phoLib2;

	vector<double> phoLib1Ene;
	vector<vector<double>> phoLib2Ene;
	double wt;

	map<int,PhoBasicLocal> nbKeyToP[16];
	map<int,PhoBasicLocal>::iterator it;

	SimpleRiboConnect();

	void setWeight(double wt){
		this->wt = wt;
	}

	double getConnectionEnergy(const LocalFrame& csA, const LocalFrame& csB){
		XYZ x(1.0, 0, 0);
		XYZ a = local2global(csA, x);
		XYZ b = csA.origin_;
		XYZ c = csB.origin_;
		XYZ d = local2global(csB, x);
		double ene;

		double dist = b.distance(c);
		double ang1 = angleX(a,b,c);
		double ang2 = angleX(b,c,d);
		double dihed=  dihedral(a,b,c,d);
		int i1 = (int)(dist*3.3333333);
		int i2 = (int)(ang1*0.2);
		int i3 = (int)(ang2*0.2);
		int i4 = (int)((dihed+180)*0.1);
		if(i2 > 35) i2 = 35;
		if(i3 > 35) i3 = 35;
		if(i4 > 35) i4 = 35;
		if(i1 < 40)
			ene = eneList[i1*46656+i2*1296+i3*36+i4];
		else {
			double e1 = eneList[39*46656+i2*1296+i3*36+i4];
			ene = e1 + (dist - 11.85)*(e1 - eneList[38*46656+i2*1296+i3*36+i4])*3.3333333;
		}
		return ene;
	}


	PhoBasicLocal getPhoBasicLocal(int bpIndex, int typeA, int typeB, const LocalFrame& csA, const LocalFrame& csB){
		it =  nbKeyToP[typeA*4+typeB].find(bpIndex);
		if(it != nbKeyToP[typeA*4+typeB].end()) {
			PhoBasicLocal pl = it->second;
			pl.e = pl.e*wt;
			//cout << "key: " << key << " " << pl.t.toString() << " " << pl.e << endl;

			return pl;
		}
		else {
			XYZ x(1.0, 0, 0);
			XYZ a = local2global(csA, x);
			XYZ b = csA.origin_;
			XYZ c = csB.origin_;
			XYZ d = local2global(csB, x);
			double ene;

			double dist = b.distance(c);
			double ang1 = angleX(a,b,c);
			double ang2 = angleX(b,c,d);
			double dihed=  dihedral(a,b,c,d);
			int i1 = (int)(dist*3.3333333);
			int i2 = (int)(ang1*0.2);
			int i3 = (int)(ang2*0.2);
			int i4 = (int)((dihed+180)*0.1);
			if(i2 > 35) i2 = 35;
			if(i3 > 35) i3 = 35;
			if(i4 > 35) i4 = 35;
			if(i1 < 40)
				ene = eneList[i1*46656+i2*1296+i3*36+i4];
			else {
				double e1 = eneList[39*46656+i2*1296+i3*36+i4];
				ene = e1 + (dist - 11.85)*(e1 - eneList[38*46656+i2*1296+i3*36+i4])*3.3333333;
			}


			double minPhoE = 99999.9;
			int minIndex = -1;
			XYZ ta, tb;
			double phoE;
			for(int i=0;i<10;i++){
				ta = local2global(csA, phoLib1[i]);
				phoE = phoLib1Ene[i];
				tb = global2local(csB, ta);
				phoE += et0.getEnergy(tb);
				if(phoE < minPhoE){
					minIndex = i;
					minPhoE = phoE;
				}
			}

			minPhoE = 99999.9;
			int minIndex2 = -1;
			for(int i=0;i<10;i++){
				ta = local2global(csA, phoLib2[minIndex][i]);
				phoE = phoLib2Ene[minIndex][i];
				tb = global2local(csB, ta);
				phoE += et0.getEnergy(tb);
				if(phoE < minPhoE){
					minIndex2 = i;
					minPhoE = phoE;
				}
			}

			//cout << "not in  map: " << minIndex << " " << minIndex2 << endl;

			return PhoBasicLocal(phoLib2[minIndex][minIndex2], (ene + minPhoE)*wt);
		}
	}


	PhoBasicLocal getPhoBasicLocal(const LocalFrame& csA, const LocalFrame& csB){

		XYZ x(1.0, 0, 0);
		XYZ a = local2global(csA, x);
		XYZ b = csA.origin_;
		XYZ c = csB.origin_;
		XYZ d = local2global(csB, x);
		double ene;

		double dist = b.distance(c);
		double ang1 = angleX(a,b,c);
		double ang2 = angleX(b,c,d);
		double dihed=  dihedral(a,b,c,d);
		int i1 = (int)(dist*3.3333333);
		int i2 = (int)(ang1*0.2);
		int i3 = (int)(ang2*0.2);
		int i4 = (int)((dihed+180)*0.1);
		if(i2 > 35) i2 = 35;
		if(i3 > 35) i3 = 35;
		if(i4 > 35) i4 = 35;
		if(i1 < 40)
			ene = eneList[i1*46656+i2*1296+i3*36+i4];
		else {
			double e1 = eneList[39*46656+i2*1296+i3*36+i4];
			ene = e1 + (dist - 11.85)*(e1 - eneList[38*46656+i2*1296+i3*36+i4])*3.3333333;
		}


		double minPhoE = 99999.9;
		int minIndex = -1;
		XYZ ta, tb;
		double phoE;
		for(int i=0;i<10;i++){
			ta = local2global(csA, phoLib1[i]);
			phoE = phoLib1Ene[i];
			tb = global2local(csB, ta);
			phoE += et0.getEnergy(tb);
			if(phoE < minPhoE){
				minIndex = i;
				minPhoE = phoE;
			}
		}

		int minIndex2 = -1;
		minPhoE = 99999.9;

		for(int i=0;i<10;i++){
			ta = local2global(csA, phoLib2[minIndex][i]);
			phoE = phoLib2Ene[minIndex][i];
			tb = global2local(csB, ta);
			phoE += et0.getEnergy(tb);
			if(phoE < minPhoE){
				minIndex2 = i;
				minPhoE = phoE;
			}
		}

		return PhoBasicLocal(phoLib2[minIndex][minIndex2], (ene + minPhoE)*wt);
	}


	virtual ~SimpleRiboConnect();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_SIMPLERIBOCONNECT_H_ */
