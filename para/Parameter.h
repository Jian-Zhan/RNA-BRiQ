/*
 * Parameter.h
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#ifndef PARA_PARAMETER_H_
#define PARA_PARAMETER_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "tools/StringTool.h"

namespace NSPpara {

using namespace std;

class Parameter {

public:
	double wtBp1;
	double wtBp2;
	double wtBp3;

	double wtRotamer;
	double wtBaseOxygen;
	double kBond;
	double kAng;
	double wtDihed;
	double wtPho;
	double wtConnect;

	string dihedEneType;
	string rotEneType;
	string atomicEneType;
	string clashType;

	double wtClash;
	double lamdaClash;
	double bbClash;

	double T0;
	double T1;
	double anneal;

	double initShift;
	double dShift;

	int stepNum;

	int outFreq;

	double initConnectWT;
	double initClashWT;
	double connectWTFactor;
	double clashWTFactor;

	bool loopRiboConnectMove;
	bool ctRandMove;
	bool f3Move;
	bool singleBaseMove;
	bool reverseRotMove;

	string spType;

	/*
	int smcModelNum;
	double smcT0;
	double smcT1;
	int smcStepNum;
	int smcClusterNum;

	double omcT0;
	double omcT1;
	int omcStepNum;
	*/

	Parameter(){
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
	}


	Parameter(const string& paraFile);

	virtual ~Parameter();
};

inline double energyRescale(double e){
	if(e < 0.81)
		return e;
	else
		return 1.8*sqrt(e) - 0.81;
}

inline double energyRescale8(double e){
	//return 0;
	if(e < 8.0)
		return e;
	else {
		return 8.0*log(e) - 8.635532333;
	}
}

} /* namespace NSPtest */

#endif /* PARA_PARAMETER_H_ */
