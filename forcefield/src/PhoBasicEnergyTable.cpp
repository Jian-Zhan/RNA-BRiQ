/*
 * PhoBasicEnergyTable.cpp
 *
 *  Created on: Aug 8, 2019
 *      Author: s2982206
 */

#include <forcefield/PhoBasicEnergyTable.h>

namespace NSPforcefield {


basePhoET0::basePhoET0() {
	string path = NSPdataio::datapath();
	string densityFile = path + "/rnaDM/pho0.density";
	ifstream file;
	file.open(densityFile, ios::in);
	if(file.fail()){
		cout << "can't open file: " << densityFile << endl;
		exit(0);
	}

	this->bpET = new double[80*128*123];
	int index = 0;
	double ene;
	while( file >> ene) {
		this->bpET[index] = ene;
		index++;
	}
	file.close();

}

double basePhoET0::getEnergy(const XYZ& t){
	double x = t.x_;
	double y = t.y_;
	double z = t.z_;

	int ix = (x+5.0)/0.1;
	int iy = (y+6.5)/0.1;
	int iz = (z+6.0)/0.1;
	if(ix < 0 || ix >=80 || iy < 0 || iy >= 128 || iz < 0 || iz >= 123)
		return this->bpET[0];

	return this->bpET[ix*15744+iy*123+iz];
}

basePhoET0::~basePhoET0(){
	delete [] this->bpET;
}

basePhoET1::basePhoET1() {
	string path = NSPdataio::datapath();
	string densityFile = path + "/rnaDM/pho1.density";
	ifstream file;
	file.open(densityFile, ios::in);
	if(file.fail()){
		cout << "can't open file: " << densityFile << endl;
		exit(0);
	}

	this->bpET = new double[57*104*107];
	int index = 0;
	double ene;
	while( file >> ene) {
		this->bpET[index] = ene;
		index++;
	}
	file.close();
}

double basePhoET1::getEnergy(const XYZ& t){
	double x = t.x_;
	double y = t.y_;
	double z = t.z_;

	int ix = (x+5.0)/0.1;
	int iy = (y+5.3)/0.1;
	int iz = (z+5.4)/0.1;
	if(ix < 0 || ix >= 57 || iy < 0 || iy >= 104 || iz < 0 || iz >= 107)
		return this->bpET[0];

	return this->bpET[ix*11128+iy*107+iz];
}

basePhoET1::~basePhoET1(){
	delete [] this->bpET;
}

phoPhoET::phoPhoET(){
	string path = NSPdataio::datapath();
	string ppFile = path+"/rnaDM/pp1.ene";
	ifstream file;
	file.open(ppFile, ios::in);
	if(file.fail()){
		cout << "can't open file: " << ppFile << endl;
		exit(0);
	}
	this->pp1 = new double[800];
	int index=  0;
	double e;
	while(file >> e){
		this->pp1[index] = e;
		index++;
	}
	file.close();

}

double phoPhoET::getPhoPhoEnergy(double d, int sep){

	if(sep == 1) {
		if(d <= 2.0)
			return (424.583-140.519*d);
		else if(d > 10.0)
			return (189.5*d-1672.47);
		else {
			int id = (d-2.0)/0.01;
			return this->pp1[id];
		}
	}
	else {
		if(d > 5.0)
			return 0;
		double u = (d-5.0)*2;
		return u*u*u*u;
	}

}

phoPhoET::~phoPhoET(){
	delete [] pp1;
}

} /* namespace NSPforcefield */
