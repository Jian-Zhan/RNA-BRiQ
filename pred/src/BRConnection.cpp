/*
 * BRConnection.cpp
 *
 *  Created on: Apr 4, 2019
 *      Author: s2982206
 */

#include <pred/BRConnection.h>

namespace NSPpred {

BRConnection::BRConnection() {
	// TODO Auto-generated constructor stub
	this->fatherNode = NULL;
	this->childNode = NULL;
	this->treeSize = 0;
	this->childOrNotChild = new bool[1];
	this->riboseInfo = new int[1];
	this->phoInfo = new int[1];
	this->ctPhoInfo = new int[1];
	this->hasThreeBaseFragment = false;
	this->f3BaseInfo = new int[1];
	this->f3RiboseInfo = new int[1];
	this->f3PhoInfo = new int[1];
	this->fixed = false;
	this->fragLib = NULL;
	this->f2Lib = NULL;
	this->f2Frag = NULL;
	this->f2FragTmp = NULL;
	this->f3Lib = NULL;
	this->childRot = NULL;
	this->childRotTmp = NULL;
	this->mvRotType = 0;
	this->tbFragmentIndex = 0;
	this->ctType = "null";

}

BRConnection::BRConnection(FragmentLibrary* fragLib, const string& connectType, BRNode* fatherNode, BRNode* childNode, int seqLen){

	this->fatherNode = fatherNode;
	this->childNode = childNode;
	this->treeSize = seqLen;
	this->childOrNotChild = new bool[treeSize];
	this->riboseInfo = new int[treeSize];
	this->phoInfo = new int[treeSize];
	this->ctPhoInfo = new int[treeSize];
	this->f3BaseInfo = new int[treeSize];
	this->f3RiboseInfo = new int[treeSize];
	this->f3PhoInfo = new int[treeSize];
	this->hasThreeBaseFragment = false;
	this->fixed = false;
	this->mvRotType = 0;
	this->cm = childNode->cs1 - fatherNode->cs1;
	this->cmTmp = this->cm;
	this->ctType = connectType;

	this->f2Lib = NULL;
	this->f3Lib = NULL;
	this->f2Frag = NULL;
	this->f2FragTmp = NULL;
	this->fragLib = fragLib;
	this->childRot = childNode->rot;
	this->childRotTmp = childNode->rot;

	if(connectType == "wc")
		this->f2Lib = fragLib->wcPair[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "nwc")
		this->f2Lib = fragLib->nwcPair[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "wcNb")
		this->f2Lib = fragLib->wcNb[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "nwcNb")
		this->f2Lib = fragLib->nwcNb[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "loopNb")
		this->f2Lib = fragLib->loopNb[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "revNb")
		this->f2Lib = fragLib->revNb[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "bulge13")
		this->f2Lib = fragLib->bulge13[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "AG")
		this->f2Lib = fragLib->agLib;
	else if(connectType == "GA")
		this->f2Lib = fragLib->gaLib;
	/*
	else if(connectType == "bulge14")
		this->f2Lib = fragLib->bulge14[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "revBulge13")
		this->f2Lib = fragLib->revBulge13[fatherNode->baseType*4 + childNode->baseType];
	else if(connectType == "revBulge14")
		this->f2Lib = fragLib->revBulge14[fatherNode->baseType*4 + childNode->baseType];
	*/
	else if(connectType == "jump")
		this->f2Lib = fragLib->jumpLib[fatherNode->baseType*4 + childNode->baseType];


	this->f2Frag = this->f2Lib->fragListLevel2[0];
	this->f2FragTmp = this->f2Frag;

	this->tbFragmentIndex = 0;
	//this->mutator = new MoveMutator(mutatorType, fatherNode->baseType, childNode->baseType);
}

BRConnection::~BRConnection() {
	this->childOrNotChild = new bool[treeSize];
	this->riboseInfo = new int[treeSize];
	this->phoInfo = new int[treeSize];
	this->ctPhoInfo = new int[treeSize];
	this->f3BaseInfo = new int[treeSize];
	this->f3RiboseInfo = new int[treeSize];
	this->f3PhoInfo = new int[treeSize];

	delete [] riboseInfo;
	delete [] phoInfo;
	delete [] ctPhoInfo;
	delete [] f3BaseInfo;
	delete [] f3RiboseInfo;
	delete [] f3PhoInfo;
	delete [] childOrNotChild;
}

} /* namespace NSPpred */
