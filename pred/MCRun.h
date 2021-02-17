/*
 * MCRun.h
 *
 *  Created on: Jul 30, 2019
 *      Author: s2982206
 */

#ifndef PRED_MCRUN_H_
#define PRED_MCRUN_H_

#include "pred/BRFoldingTree.h"

namespace NSPpred {

using namespace NSPmodel;
using namespace std;

class MCRun {

	//BRFoldingTreeBasic* simpFt;


	//BRTreeInfoBasic* initBasic;

public:
	BRFoldingTree* ft;
	BRTreeInfo* init;
	MCRun(BRFoldingTree* ft){
		this->ft = ft;
		init = ft->getTreeInfo();
	}

	void generateDecoysRandInit(const string& output);

	void optimizeFromInit(const string& keyFile, const string& outFilePrefix, const int startID);

	void simpleMC(const string& outpdb, bool traj);

	void optimize(double t0, double kStep);

	void optimizeBackbone(const string& output);


	void debug();

	void test();

	virtual ~MCRun();
};

} /* namespace NSPpred */

#endif /* PRED_MCRUN_H_ */
