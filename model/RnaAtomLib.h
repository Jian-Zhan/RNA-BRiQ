/*
 * RnaAtomLib.h
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#ifndef MODEL_RNAATOMLIB_H_
#define MODEL_RNAATOMLIB_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "model/AtomProperty.h"
#include "model/RNABaseName.h"
#include "dataio/datapaths.h"

namespace NSPmodel {


class RnaAtomLib {
private:
	vector<string> uniqueAtomNames;
	map<string,int> uniqueNamesToUniqueID;
	vector<AtomProperty*> properties;
	vector<vector<int>*> baseUniqueIDs;
	vector<vector<string>*> baseAtomNames;
	vector<vector<string>*> baseScAtomNames;
	RNABaseName rn;

public:
	RnaAtomLib();
	AtomProperty* getAtomProperty(string& uniqueName) const;
	vector<string>* getAtomNames(int type) const;
	vector<string>* getSidechainAtoms(int type) const;
	string getUniqueName(int id);
	int getUniqueID(string& uniqueName);

	virtual ~RnaAtomLib();
};

} /* namespace NSPmodel */

#endif /* MODEL_RNAATOMLIB_H_ */
