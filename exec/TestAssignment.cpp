/*
 * TestAssignment.cpp
 *
 *  Created on: 2020Äê11ÔÂ18ÈÕ
 *      Author: pengx
 */

#include "pred/BRFoldingTree.h"
#include "pred/BRFoldingTreeBasic.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "model/AssignRNASS.h"
#include "forcefield/NeighborConnectionGenerator.h"
#include "pred/MCRun.h"
#include <time.h>
#include <stdlib.h>
#include <iostream>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;

int main(int argc, char** argv){
	/*
	 * input: $PDBFILE
	 */

	if(argc != 3 || argv[1] == "-h")
	{
		cout << "Usage: BRiQ_assignSS $PDBFILE $OUTFILE" << endl;
		exit(0);
	}

	//cout << "start" << endl;
	string pdbFile = string(argv[1]);
	string outFile = string(argv[2]);

	RNAPDB* pdb = new RNAPDB(pdbFile, "xxxx");
	RnaAtomLib* atLib = new RnaAtomLib();

	//cout << "init ar" << endl;
	AssignRNASS* ar = new AssignRNASS(pdb, atLib);
	//cout << "assign " << endl;
	ar->printInfo(outFile);

	delete pdb;
	delete atLib;
	delete ar;
}


