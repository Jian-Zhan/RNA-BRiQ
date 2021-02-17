/*
 * AtomicEnergyTable.cpp
 *
 *  Created on: Jul 27, 2019
 *      Author: s2982206
 */

#include <forcefield/AtomicEnergyTable.h>

namespace NSPforcefield {


AtomicEnergyTable::AtomicEnergyTable(Parameter* para) {
	// TODO Auto-generated constructor stub

	this->para = para;
	string path = NSPdataio::datapath();
	string fileName = path+"clashDistance.txt";


	vector<double> dList1;
	vector<double> dList2;
	vector<double> dList3;

	ifstream file;
	file.open(fileName, ios::in);
	if(file.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}

	string line;

	double e;
	vector<string> spt;
	while(getline(file, line)){
		splitString(line, " ", &spt);
		dList1.push_back(atof(spt[0].c_str()));
		dList2.push_back(atof(spt[1].c_str()));
		dList3.push_back(atof(spt[2].c_str()));
	}
	file.close();

	double d1, d2, d3;
	double d, e1, e2, e3;

	for(double d0=0.0;d0<=4.005;d0+=0.01){
		vector<double> eList;
		for(int j=0;j<1701;j++){
			d = sqrt(j*0.01);
			e1 = vdwEnergy(d0, d);
			eList.push_back(e1);
		}
		this->vdwEt.push_back(eList);
	}

	for(int i=0;i<2209;i++) {
		vector<double> eList1;
		vector<double> eList2;
		vector<double> eList3;
		d1 = dList1[i];
		d2 = dList2[i];
		d3 = dList3[i];
		//cout << d1 << " " << d2 << " " << d3 << endl;
		for(int j=0;j<3601;j++){
			d = sqrt(j*0.01);

			e1 = vdwEnergy(d1, d);
			e2 = vdwEnergy(d2, d);
			e3 = vdwEnergy(d3, d);

			eList1.push_back(e1);
			eList2.push_back(e2);
			eList3.push_back(e3);
		}
		//cout << eList1[1200] << " " << eList2[1200] << " " << eList3[1200] << endl;
		et1.push_back(eList1);
		et2.push_back(eList2);
		et3.push_back(eList3);
	}


	fileName = path + "backboneHbond/O2O2.ene";
	file.open(fileName, ios::in);
	if(file.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}
	while(file >> e){
		this->O2O2HbondEnergy.push_back(energyRescale(e + 2.5)); //minE: -3
	}
	file.close();

	fileName = path + "backboneHbond/O2OP.ene";
	file.open(fileName, ios::in);
	if(file.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}
	while(file >> e){
		this->O2OPHbondEnergy.push_back(energyRescale(e + 3.5)); //minE: -2
	}
	file.close();

	fileName = path + "backboneHbond/OPOP.ene";
	file.open(fileName, ios::in);
	if(file.fail()){
		cout << "can't open file: " << fileName << endl;
		exit(0);
	}
	while(file >> e){
		this->OPOPHbondEnergy.push_back(energyRescale(e + 2.0)); //minE: -1.5
	}
	file.close();

	vector<string> typeList;
	typeList.push_back("A_");
	typeList.push_back("U_");
	typeList.push_back("G_");
	typeList.push_back("C_");
	vector<string> baseAtoms;
	fileName = path + "baseAtomClash/baseAtoms";
	file.open(fileName.c_str(), ios::in);
	if(file.fail()) {
		cout << "can't open " + fileName << endl;
		exit(1);
	}
	while(getline(file, line)){
		baseAtoms.push_back(line);
	}
	file.close();

	int index;
	double clashDistance;
	for(int i=0;i<4;i++){
		for(int j=0;j<22;j++){
			string fn = path + "baseAtomClash/clash" + para->clashType + "/"+typeList[i]+baseAtoms[j]+".dist";
			file.open(fn.c_str(), ios::in);
			if(file.fail()){
				cout << "can't open " + fn << endl;
				exit(1);
			}
			while(file >> index >> clashDistance){

				this->baseAtomClashDistance[i*22+j].push_back(clashDistance);
				this->baseAtomNearestAtomIndex[i*22+j].push_back(index);
			}
			file.close();
		}
	}



	for(int i=0;i<11;i++)
		baseAtomIndexToClashAtomIndex[i] = -1;
	for(int i=0;i<11;i++)
		baseAtomIndexToClashAtomIndex[100+i] = -1;
	for(int i=0;i<11;i++)
		baseAtomIndexToClashAtomIndex[200+i] = -1;
	for(int i=0;i<11;i++)
		baseAtomIndexToClashAtomIndex[300+i] = -1;

	baseAtomIndexToClashAtomIndex[2] = 0;
	baseAtomIndexToClashAtomIndex[5] = 1;
	baseAtomIndexToClashAtomIndex[6] = 2;
	baseAtomIndexToClashAtomIndex[8] = 3;
	baseAtomIndexToClashAtomIndex[1] = 4;
	baseAtomIndexToClashAtomIndex[7] = 5;

	baseAtomIndexToClashAtomIndex[102] = 6;
	baseAtomIndexToClashAtomIndex[103] = 7;
	baseAtomIndexToClashAtomIndex[105] = 8;
	baseAtomIndexToClashAtomIndex[106] = 9;
	baseAtomIndexToClashAtomIndex[107] = 10;

	baseAtomIndexToClashAtomIndex[202] = 11;
	baseAtomIndexToClashAtomIndex[205] = 12;
	baseAtomIndexToClashAtomIndex[206] = 13;
	baseAtomIndexToClashAtomIndex[208] = 14;
	baseAtomIndexToClashAtomIndex[209] = 15;
	baseAtomIndexToClashAtomIndex[201] = 16;


	baseAtomIndexToClashAtomIndex[302] = 17;
	baseAtomIndexToClashAtomIndex[303] = 18;
	baseAtomIndexToClashAtomIndex[305] = 19;
	baseAtomIndexToClashAtomIndex[306] = 20;
	baseAtomIndexToClashAtomIndex[307] = 21;

}



AtomicEnergyTable::~AtomicEnergyTable() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
