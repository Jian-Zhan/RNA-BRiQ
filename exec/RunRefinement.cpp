/*
 * RunRefinement.cpp
 *
 *  Created on: 2020��12��23��
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

void checkInputFile(const string& inputFile){
	ifstream file;
	cout << "check input file: " << endl;
	file.open(inputFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open inputFile: " << inputFile << endl;
		exit(0);
	}
	file.close();

	NSPtools::InputParser input(inputFile);
	string pdbFile = input.getValue("pdb");
	string baseSeq = input.getValue("seq");
	string baseSec = input.getValue("sec");
	string nwcSec = input.getValue("nwc");

	if(pdbFile == ""){
		cout << "invalid inputFile: initial pdb required" << endl;
		exit(0);
	}

	if(baseSeq == ""){
		cout << "invalid inputFile: base sequence required" << endl;
		exit(0);
	}

	if(baseSec == ""){
		cout << "invalid inputFile: Watson-Crick pairing sequence required" << endl;
		exit(0);
	}

	if(nwcSec == ""){
		cout << "invalid inputFile: non-Watson-Crick pairing sequence required" << endl;
		exit(0);
	}

	if(baseSec.length() != baseSeq.length()){
		cout << "the length of Watson-Crick pairing sequence must equal to the length of base sequence" << endl;
		exit(0);
	}

	if(nwcSec.length() != baseSeq.length()){
		cout << "the length of non-Watson-Crick pairing sequence must equal to the length of base sequence" << endl;
		exit(0);
	}


	file.open(pdbFile.c_str(), ios::in);
	if(!file.is_open()){
		cout << "fail to open pdb file: " << pdbFile << endl;
		exit(0);
	}
}

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

	checkInputFile(inputFile);

	cout << "init energy table:" << endl;
	EnergyTable* et = new EnergyTable();
	double t0 = 0.5;
	double kStep = 1.0;

	cout << "run refinement: " << endl;
	BRFoldingTree* ft = new BRFoldingTree(inputFile, et);
	MCRun mc(ft);
	mc.optimize(t0, kStep);

	BRTreeInfo* info = ft->getTreeInfo();
	info->printPDB(outPDB);

}


