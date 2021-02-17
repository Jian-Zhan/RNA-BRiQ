/*
 * BuildBackbone2.cpp
 *
 *  Created on: Oct 9, 2019
 *      Author: s2982206
 */


#include "pred/BackboneModelingTemplate.h"
#include "pred/MCRun.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "forcefield/NeighborConnectionGenerator.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;

int main(int argc, char** argv){

	string pdbFile = string(argv[1]);
	string paraFile = string(argv[2]);
	string output = string(argv[3]);

	cout << "start: " << endl;



	Parameter* para = new Parameter(paraFile);

	cout << "init bm" << endl;
	BackboneModelingTemplate bm(pdbFile, paraFile);

	cout << "runMC" << endl;
	//bm.printEnergyDetail();
	//bm.debug();
	bm.runMC();

	//bm.printEnergyDetail();

	//bm.printDiheds(output);

	BRTreeInfo* info = bm.toTreeInfo();
	info->printPDB(output);



}


