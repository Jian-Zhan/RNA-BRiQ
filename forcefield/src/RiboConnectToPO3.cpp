/*
 * RiboConnectToPO3.cpp
 *
 *  Created on: Aug 5, 2019
 *      Author: s2982206
 */

#include <forcefield/RiboConnectToPO3.h>

namespace NSPforcefield {

RiboConnectToPO3::RiboConnectToPO3() {
	// TODO Auto-generated constructor stub
	this->para = new Parameter();
	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	string fileName = path+"impD1D2.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	vector<string> spt;
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD1D2.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();

	fileName = path+"impD4D5.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD4D5.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();

	fileName = path+"D2D4D3.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eD2D4D3.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();


	c2 = XYZ(-2.008, -1.402, 0);
	c3 = XYZ(-1.418, 0, 0);
	o3 = XYZ(0, 0, 0);

	fileName = path+"d1d2Lib/rot-level1.txt";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	getline(file,s);
	double x,y;
	while(getline(file,s)){
		NSPtools::splitString(s," ",&spt);
		x = atof(spt[0].c_str());
		y = atof(spt[1].c_str());
		this->d1d2Lib1.push_back(XYZ(x,y,0));
	}
	file.close();

	for(int i=0;i<d1d2Lib1.size();i++){
		char ss[10];
		sprintf(ss, "%d", i);
		fileName = path + "d1d2Lib/rot-level2-"+string(ss)+".txt";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()){
			cout << "fail to open file: " << fileName << endl;
			exit(1);
		}
		getline(file,s);
		NSPtools::splitString(s," ",&spt);
		double error = atof(spt[2].c_str());
		this->lib2Error.push_back(error);
		double x,y;
		vector<XYZ> subList;
		while(getline(file,s)){
			NSPtools::splitString(s," ",&spt);
			x = atof(spt[0].c_str());
			y = atof(spt[1].c_str());
			subList.push_back(XYZ(x,y,0));
		}
		file.close();
		this->d1d2Lib2.push_back(subList);
	}

    dihedList1.push_back(200.0);
    dihedList1.push_back(230.0);
    dihedList1.push_back(260.0);
    dihedList1.push_back(290.0);
    dihedList1.push_back(320.0);
    dihedList1.push_back(40.0);
    dihedList1.push_back(60.0);
    dihedList1.push_back(70.0);
    dihedList1.push_back(80.0);
    dihedList1.push_back(90.0);
    dihedList1.push_back(100.0);
    dihedList1.push_back(110.0);
    dihedList1.push_back(120.0);
    dihedList1.push_back(130.0);
    dihedList1.push_back(140.0);
    dihedList1.push_back(160.0);


    dihedList2.push_back(160.0);
    dihedList2.push_back(130.0);
    dihedList2.push_back(100.0);
    dihedList2.push_back(70.0);
    dihedList2.push_back(40.0);
    dihedList2.push_back(320.0);
    dihedList2.push_back(300.0);
    dihedList2.push_back(290.0);
    dihedList2.push_back(280.0);
    dihedList2.push_back(270.0);
    dihedList2.push_back(260.0);
    dihedList2.push_back(250.0);
    dihedList2.push_back(240.0);
    dihedList2.push_back(230.0);
    dihedList2.push_back(220.0);
    dihedList2.push_back(200.0);


}

RiboConnectToPO3::RiboConnectToPO3(Parameter* para) {
	// TODO Auto-generated constructor stub
	this->para = para;
	string path = NSPdataio::datapath();
	ifstream file;
	string s;
	string fileName = path+"impD1D2.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	vector<string> spt;
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD1D2.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();

	fileName = path+"impD4D5.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eImpD4D5.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();

	fileName = path+"D2D4D3.ene";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	while(getline(file,s)){
		NSPtools::splitString(s, " ", &spt);
		for(int i=0;i<180;i++){
			this->eD2D4D3.push_back(atof(spt[i].c_str()));
		}
	}
	file.close();


	c2 = XYZ(-2.008, -1.402, 0);
	c3 = XYZ(-1.418, 0, 0);
	o3 = XYZ(0, 0, 0);

	fileName = path+"d1d2Lib/rot-level1.txt";
	file.open(fileName.c_str(), ios::in);
	if(! file.is_open()){
		cout << "fail to open file: " << fileName << endl;
		exit(1);
	}
	getline(file,s);
	double x,y;
	while(getline(file,s)){
		NSPtools::splitString(s," ",&spt);
		x = atof(spt[0].c_str());
		y = atof(spt[1].c_str());
		this->d1d2Lib1.push_back(XYZ(x,y,0));
	}
	file.close();

	for(int i=0;i<d1d2Lib1.size();i++){
		char ss[10];
		sprintf(ss, "%d", i);
		fileName = path + "d1d2Lib/rot-level2-"+string(ss)+".txt";
		file.open(fileName.c_str(), ios::in);
		if(! file.is_open()){
			cout << "fail to open file: " << fileName << endl;
			exit(1);
		}
		getline(file,s);
		NSPtools::splitString(s," ",&spt);
		double error = atof(spt[2].c_str());
		this->lib2Error.push_back(error);
		double x,y;
		vector<XYZ> subList;
		while(getline(file,s)){
			NSPtools::splitString(s," ",&spt);
			x = atof(spt[0].c_str());
			y = atof(spt[1].c_str());
			subList.push_back(XYZ(x,y,0));
		}
		file.close();
		this->d1d2Lib2.push_back(subList);
	}

    dihedList1.push_back(200.0);
    dihedList1.push_back(230.0);
    dihedList1.push_back(260.0);
    dihedList1.push_back(290.0);
    dihedList1.push_back(320.0);
    dihedList1.push_back(40.0);
    dihedList1.push_back(60.0);
    dihedList1.push_back(70.0);
    dihedList1.push_back(80.0);
    dihedList1.push_back(90.0);
    dihedList1.push_back(100.0);
    dihedList1.push_back(110.0);
    dihedList1.push_back(120.0);
    dihedList1.push_back(130.0);
    dihedList1.push_back(140.0);
    dihedList1.push_back(160.0);


    dihedList2.push_back(160.0);
    dihedList2.push_back(130.0);
    dihedList2.push_back(100.0);
    dihedList2.push_back(70.0);
    dihedList2.push_back(40.0);
    dihedList2.push_back(320.0);
    dihedList2.push_back(300.0);
    dihedList2.push_back(290.0);
    dihedList2.push_back(280.0);
    dihedList2.push_back(270.0);
    dihedList2.push_back(260.0);
    dihedList2.push_back(250.0);
    dihedList2.push_back(240.0);
    dihedList2.push_back(230.0);
    dihedList2.push_back(220.0);
    dihedList2.push_back(200.0);

}

double RiboConnectToPO3::printEnergyDetail(XYZ* riboCoordsA, XYZ* riboCoordsB, PhophateGroup& pho, double impA, double impB){

	LocalFrame csA(riboCoordsA[1], riboCoordsA[2], riboCoordsA[6]);
	LocalFrame cs0;

	c2 = riboCoordsA[1];
	c3 = riboCoordsA[2];
	o3 = riboCoordsA[6];
	o2 = riboCoordsA[5];
	c5 = riboCoordsB[7];
	c4 = riboCoordsB[3];
	o4 = riboCoordsB[4];
	LocalFrame cs1;
	LocalFrame cs2;
	double e, u;
	double eBond, eAng1, eAng2, eD12, eD45, eD3;

	double dihed1, dihed2;
	double o2o4, c2o4;
	o2o4 = o2.distance(o4);
	c2o4 = c2.distance(o4);
	int idA, idB;
	idA = (int)(o2o4*50);
	idB = (int)(c2o4*50);
	if(idA > 499) idA = 499;
	if(idB > 499) idB = 499;

	p = pho.tList[0];
	o5 = pho.tList[1];

	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;

	xd3 = o5.distance(c5);
	xang3 = angleX(p, o5, c5);
	xang4 = angleX(o5, c5, c4);
	dihed1 = dihedral(c2, c3, o3, p);
	dihed2 = dihedral(c3, o3, p, o5);
	xdihed3 = dihedral(o3, p, o5, c5);
	xdihed4 = dihedral(p, o5, c5, c4);
	xdihed5 = dihedral(o5, c5, c4, o4);
	if(dihed1 < 0) dihed1 += 360;
	if(dihed2 < 0) dihed2 += 360;
	if(xdihed3 < 0) xdihed3 += 360;
	if(xdihed4 < 0) xdihed4 += 360;
	if(xdihed5 < 0) xdihed5 += 360;
	e = 0;
	u = (xd3-len3)*para->kBond;
	if(u<1 && u>-1)
		eBond = u*u;
	else if(u > 1)
		eBond = (2*u-1);
	else
		eBond = (-2*u-1);

	u = (xang3-ang3)*para->kAng;
	if(u<1 && u>-1)
		eAng1 = u*u;
	else if(u > 1)
		eAng1 = (2*u-1);
	else
		eAng1 = (-2*u-1);
	u = (xang4-ang4)*para->kAng;
	if(u<1 && u>-1)
		eAng2 = u*u;
	else if(u > 1)
		eAng2 = (2*u-1);
	else
		eAng2 = (-2*u-1);


	int impIndexA = (int)((impA+60)*0.166666667);
	if(impIndexA < 0) impIndexA = 0;
	if(impIndexA > 19) impIndexA = 19;

	int impIndexB = (int)((impB+60)*0.166666667);
	if(impIndexB < 0) impIndexB = 0;
	if(impIndexB > 19) impIndexB = 19;

	eD12 = eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
	eD45 = eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
	eD3 = eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

	c2 = XYZ(-2.008, -1.402, 0);
	c3 = XYZ(-1.418, 0, 0);
	o3 = XYZ(0, 0, 0);

	//printf("bond: %7.3f ang1: %7.3f ang2: %7.3f D1: %6.3f D2: %6.3f D3: %6.3f D4: %6.3f D5: %6.3f ImpA: %6.3f ImpB: %6.3f\n", xd3, xang3, xang4, dihed1, dihed2, xdihed3, xdihed4, xdihed5, impA, impB);
	printf("eBond: %7.3f eAng1: %7.3f eAng2: %7.3f eD12: %7.3f eD45: %7.3f eD3: %7.3f\n", eBond, eAng1, eAng2, eD12, eD45, eD3);
	return  eBond + eAng1 + eAng2 + eD12 + eD45 + eD3;
}

PhophateGroupLocal RiboConnectToPO3::getPO3Local(XYZ* riboCoordsA, XYZ* riboCoordsB, double impA, double impB){

	double d0 = riboCoordsA[6].distance(riboCoordsB[7]);
	/*
	if(d0 > 4.5) {
		double u = (d0-3.9)*30;
		PhophateGroupLocal pho;
		pho.e = u*u;
		return pho;
	}
	*/

	LocalFrame csA(riboCoordsA[1], riboCoordsA[2], riboCoordsA[6]);
	LocalFrame cs0;

	o2 = global2local(csA, riboCoordsA[5]);
	c5 = global2local(csA, riboCoordsB[7]);
	c4 = global2local(csA, riboCoordsB[3]);
	o4 = global2local(csA, riboCoordsB[4]);
	LocalFrame cs1;
	LocalFrame cs2;

	LocalFrame bestCs1, bestCs2;
	double bestDihed1, bestDihed2;
	double localBestD1, localBestD2;


	double minE = 999999999.9;
	double localMinE = 9999.9;
	double e, u;
	double dihed1, dihed2;
	double o2o4, c2o4;
	int idA, idB;

	o2o4 = o2.distance(o4);
	c2o4 = c2.distance(o4);
	idA = (int)(o2o4*50);
	idB = (int)(c2o4*50);
	if(idA > 499) idA = 499;
	if(idB > 499) idB = 499;

	double xd3, xang3, xang4, xdihed3, xdihed4, xdihed5;

	int impIndexA = (int)((impA+60)*0.166666667);
	if(impIndexA < 0) impIndexA = 0;
	if(impIndexA > 19) impIndexA = 19;

	int impIndexB = (int)((impB+60)*0.166666667);
	if(impIndexB < 0) impIndexB = 0;
	if(impIndexB > 19) impIndexB = 19;
	int indexD1, indexD2, indexD3, indexD4, indexD5;

	int d1d2LibSize = this->d1d2Lib1.size();
	int bestIndex1=0, bestIndex2=0;

	for(int i=0;i<d1d2LibSize;i++){
		dihed1 = d1d2Lib1[i].x_;
		dihed2 = d1d2Lib1[i].y_;
		cs1 = cs0.csNext(len1, ang1, dihed1);
		cs2 = cs1.csNext(len2, ang2, dihed2);
		p = cs1.origin_;
		o5 = cs2.origin_;
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			bestIndex1 = i;
			minE = e;
		}
	}

	if(d0 > 4.0){
		cs1 = cs0.csNext(len1, ang1, bestDihed1);
		cs2 = cs1.csNext(len2, ang2, bestDihed2);
		return PhophateGroupLocal(cs1.origin_, cs2.origin_, local2global(cs2, XYZ(-2.062, 0.570, 1.278)), local2global(cs2, XYZ(-2.054, 0.589, -1.275)), minE);
	}

	for(int i=0;i<d1d2LibSize;i++){
		dihed1 = d1d2Lib2[bestIndex1][i].x_;
		dihed2 = d1d2Lib2[bestIndex1][i].y_;
		cs1 = cs0.csNext(len1, ang1, dihed1);
		cs2 = cs1.csNext(len2, ang2, dihed2);
		p = cs1.origin_;
		o5 = cs2.origin_;
		xd3 = o5.distance(c5);
		xang3 = angleX(p, o5, c5);
		xang4 = angleX(o5, c5, c4);
		xdihed3 = dihedral(o3, p, o5, c5);
		xdihed4 = dihedral(p, o5, c5, c4);
		xdihed5 = dihedral(o5, c5, c4, o4);
		if(xdihed3 < 0) xdihed3 += 360;
		if(xdihed4 < 0) xdihed4 += 360;
		if(xdihed5 < 0) xdihed5 += 360;
		e = 0;
		u = (xd3-len3)*para->kBond;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang3-ang3)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);
		u = (xang4-ang4)*para->kAng;
		if(u<1 && u>-1)
			e += u*u;
		else if(u > 1)
			e += (2*u-1);
		else
			e += (-2*u-1);

		e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
		e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
		e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

		if(e < minE){
			bestDihed1 = dihed1;
			bestDihed2 = dihed2;
			minE = e;
		}
	}

	if(bestDihed1 < 220) {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<9;i++){
			dihed1 = localBestD1 + (i-4);
			for(int j=0;j<9;j++){
				dihed2 = localBestD2 + (j-4);
				cs1 = cs0.csNext(len1, ang1, dihed1);
				cs2 = cs1.csNext(len2, ang2, dihed2);
				p = cs1.origin_;
				o5 = cs2.origin_;
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = 0;
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}

		cs1 = cs0.csNext(len1, ang1, bestDihed1);
		cs2 = cs1.csNext(len2, ang2, bestDihed2);
		return PhophateGroupLocal(cs1.origin_, cs2.origin_, local2global(cs2, XYZ(-2.062, 0.570, 1.278)), local2global(cs2, XYZ(-2.054, 0.589, -1.275)), minE);
	}
	else {
		localBestD1 = bestDihed1;
		localBestD2 = bestDihed2;
		for(int i=0;i<11;i++){
			dihed1 = localBestD1 + (i-5)*2.0;
			for(int j=0;j<11;j++){
				dihed2 = localBestD2 + (j-5)*2.0;
				cs1 = cs0.csNext(len1, ang1, dihed1);
				cs2 = cs1.csNext(len2, ang2, dihed2);
				p = cs1.origin_;
				o5 = cs2.origin_;
				xd3 = o5.distance(c5);
				xang3 = angleX(p, o5, c5);
				xang4 = angleX(o5, c5, c4);
				xdihed3 = dihedral(o3, p, o5, c5);
				xdihed4 = dihedral(p, o5, c5, c4);
				xdihed5 = dihedral(o5, c5, c4, o4);
				if(xdihed3 < 0) xdihed3 += 360;
				if(xdihed4 < 0) xdihed4 += 360;
				if(xdihed5 < 0) xdihed5 += 360;
				e = 0;
				u = (xd3-len3)*para->kBond;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang3-ang3)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);
				u = (xang4-ang4)*para->kAng;
				if(u<1 && u>-1)
					e += u*u;
				else if(u > 1)
					e += (2*u-1);
				else
					e += (-2*u-1);

				e += eImpD1D2[impIndexA*32400 + ((int)(dihed1*0.5))*180 + (int)(dihed2*0.5)] * para->wtDihed;
				e += eImpD4D5[impIndexB*32400 + ((int)(xdihed4*0.5))*180 + (int)(xdihed5*0.5)] * para->wtDihed;
				e += eD2D4D3[((int)(dihed2*0.166666666))*10800 + ((int)(xdihed4*0.166666666))*180 + (int)(xdihed3*0.5)] * para->wtDihed;

				if(e < minE){
					bestDihed1 = dihed1;
					bestDihed2 = dihed2;
					minE = e;
				}
			}
		}

		cs1 = cs0.csNext(len1, ang1, bestDihed1);
		cs2 = cs1.csNext(len2, ang2, bestDihed2);
		return PhophateGroupLocal(cs1.origin_, cs2.origin_, local2global(cs2, XYZ(-2.062, 0.570, 1.278)), local2global(cs2, XYZ(-2.054, 0.589, -1.275)), minE);
	}
}


RiboConnectToPO3::~RiboConnectToPO3() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
