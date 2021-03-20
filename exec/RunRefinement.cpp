/*
 * RunRefinement.cpp
 *
 *  Created on: 2020Äê12ÔÂ23ÈÕ
 *      Author: pengx
 */

#include "pred/BRFoldingTree.h"
#include "pred/BRFoldingTreeBasic.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "forcefield/NeighborConnectionGenerator.h"
#include "pred/MCRun.h"
#include <iostream>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;

void printHelp(){
	cout << "Usage: BRiQ_Refinement $INPUT_FILE $OUTPUT_PDB $RANDOM_SEED" << endl;
}

int main(int argc, char** argv){

	if(argc != 4 || argv[1] == "-h")
	{
		printHelp();
		exit(0);
	}

	string inputFile = string(argv[1]);
	string outPDB = string(argv[2]);
	srand(atoi(argv[3]));
	EnergyTable* et = new EnergyTable();
	double t0 = 0.5;
	double kStep = 1.0;

	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	MCRun mc(ft);
	mc.optimize(t0, kStep);

	BRTreeInfo* info = ft->getTreeInfo();
	info->printPDB(outPDB);

}


