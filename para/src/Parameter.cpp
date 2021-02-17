/*
 * Parameter.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#include "para/Parameter.h"

namespace NSPpara {

Parameter::Parameter(const string& paraFile) {
	ifstream input;
	input.open(paraFile.c_str(), ios::in);
	if(!input.is_open()) {
		cout << "fail to open file: " << paraFile << endl;
		exit(1);
	}
	string s;
	vector<string> spt;
	this->wtBp1 = 0.4;
	this->wtBp2 = 0.8;
	this->wtBp3 = 1.0;
	this->wtBaseOxygen = 1.2;
	this->wtPho = 0.7;
	this->wtRotamer = 0.7;
	this->kAng = 0.1;
	this->kBond = 5.0;
	this->bbClash = 0.4;
	this->wtDihed = 1.2;
	this->wtConnect = 1.0;

	this->dihedEneType = "wt10";
	this->rotEneType = "wt10";


	clashType = "1";
	atomicEneType = "1";
	this->wtClash = 1.0;

	this->T0 = 2.5;
	this->T1 = 0.1;
	this->anneal = 0.95;
	this->stepNum = 10;
	this->outFreq = 1000;

	this->lamdaClash = 2.8;

	this->initConnectWT = 0.1;
	this->initClashWT = 0.05;

	this->initShift = 0.6;
	this->dShift = 0.02;

	this->connectWTFactor = 1.05;
	this->clashWTFactor = 1.05;

	loopRiboConnectMove = true;
	ctRandMove = true;
	f3Move = true;
	singleBaseMove = true;
	reverseRotMove = true;

	this->spType = "sp2000";

	while(getline(input,s)) {
		NSPtools::splitString(s, " ", &spt);
		if(spt[0] == "bp1")
			this->wtBp1 = atof(spt[1].c_str());
		else if(spt[0] == "bp2")
			this->wtBp2 = atof(spt[1].c_str());
		else if(spt[0] == "bp3")
			this->wtBp3 = atof(spt[1].c_str());
		else if(spt[0] == "rot")
			this->wtRotamer = atof(spt[1].c_str());
		else if(spt[0] == "pho")
			this->wtPho = atof(spt[1].c_str());
		else if(spt[0] == "oxy")
			this->wtBaseOxygen = atof(spt[1].c_str());
		else if(spt[0] == "bond")
			this->kBond = atof(spt[1].c_str());
		else if(spt[0] == "ang")
			this->kAng = atof(spt[1].c_str());
		else if(spt[0] == "clash")
			this->wtClash = atof(spt[1].c_str());
		else if(spt[0] == "kclash")
			this->lamdaClash = atof(spt[1].c_str());
		else if(spt[0] == "dihed")
			this->wtDihed = atof(spt[1].c_str());
		else if(spt[0] == "dihedEneType")
			this->dihedEneType = spt[1];
		else if(spt[0] == "rotEneType")
			this->rotEneType = spt[1];
		else if(spt[0] == "con")
			this->wtConnect = atof(spt[1].c_str());
		else if(spt[0] == "T0")
			this->T0 = atof(spt[1].c_str());
		else if(spt[0] == "T1")
			this->T1 = atof(spt[1].c_str());
		else if(spt[0] == "step")
			this->stepNum = atoi(spt[1].c_str());
		else if(spt[0] == "atomic")
			this->atomicEneType = spt[1].c_str();
		else if(spt[0] == "initConnect")
			this->initConnectWT = atof(spt[1].c_str());
		else if(spt[0] == "initClash")
			this->initClashWT = atof(spt[1].c_str());
		else if(spt[0] == "connectFactor")
			this->connectWTFactor = atof(spt[1].c_str());
		else if(spt[0] == "clashFactor")
			this->clashWTFactor = atof(spt[1].c_str());
		else if(spt[0] == "bbClash")
			this->bbClash = atof(spt[1].c_str());
		else if(spt[0] == "clashType")
			this->clashType = spt[1];
		else if(spt[0] == "anneal")
			this->anneal = atof(spt[1].c_str());
		else if(spt[0] == "loopRiboConnectMove"){
			if(spt[1] == "false")
				this->loopRiboConnectMove = false;
			else
				this->loopRiboConnectMove = true;
		}
		else if(spt[0] == "ctRandMove"){
			if(spt[1] == "false")
				this->ctRandMove = false;
			else
				this->ctRandMove = true;
		}
		else if(spt[0] == "f3Move"){
			if(spt[1] == "false")
				this->f3Move = false;
			else
				this->f3Move = true;
		}
		else if(spt[0] == "singleBaseMove"){
			if(spt[1] == "false")
				this->singleBaseMove = false;
			else
				this->singleBaseMove = true;
		}
		else if(spt[0] == "reverseRotMove"){
			if(spt[1] == "false")
				this->reverseRotMove = false;
			else
				this->reverseRotMove = true;
		}
		else if(spt[0] == "shift")
			this->initShift = atof(spt[1].c_str());
		else if(spt[0] == "dshift")
			this->dShift = atof(spt[1].c_str());
		else if(spt[0] == "sp")
			this->spType = spt[1];
		else{
			cout << "unknow tag " << spt[0] << endl;
		}
	}
}

Parameter::~Parameter() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPtest */
