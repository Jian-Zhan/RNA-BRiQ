/*
 * RunPredict.cpp
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
#include <time.h>
#include <stdlib.h>
#include <iostream>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;


void printHelp(){
	cout << "Usage: BRiQ_Predict $INPUT_FILE $OUTPUT_PDB $RANDOM_SEED" << endl;
}

int main(int argc, char** argv){

	//Usage: BRiQ_Predict $INPUT_FILE $OUTPUT_PDB $RANDOM_SEED
	if(argc != 4 || argv[1] == "-h")
	{
		printHelp();
		exit(0);
	}

	string inputFile = string(argv[1]);
	string outputPDB = string(argv[2]);
	int randSeed = atoi(argv[3]);
	srand(randSeed);

	cout << "init energy table:" << endl;
	EnergyTable* et = new EnergyTable();

	cout << "init folding tree" << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);

	BRTreeInfo* info = ft->getTreeInfo();

	cout << "run mc: " << endl;

	MCRun mc(ft);
	mc.simpleMC(outputPDB, false);

	delete et;
	delete ft;

}


