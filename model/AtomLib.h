/*
 * AtomLib.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ATOMLIB_H_
#define DESIGNSEQ_ATOMLIB_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "model/AtomProperty.h"
#include "model/ResName.h"
#include "dataio/datapaths.h"

namespace NSPmodel {

using namespace std;

class AtomLib {
private:
	vector<string> uniqueAtomNames;
	vector<AtomProperty*> properties;
	map<string,int> uniqueNameToUniqueID;
	vector<vector<int>*> aaUniqueIDs;
	vector<vector<string>*> aaAtomNames;
	vector<vector<string>*> aaScAtomNames;
	ResName rn;
public:
	AtomLib();
	bool atomDefined(const string& uniqueName) const;
	int uniqueNameToID(const string& uniqueName) const;
	string uniqueIDToName(int uniqueID) const;
	vector<int>* getAminoAcidAtomIDs(const string& triName);
	vector<string>* getAminoAcidAtomNames(const string& triName);
	vector<string>* getAminoAcidSidechainAtomNames(const string& triName);
	AtomProperty* getAtomProperty(const string& uniqueName);
	virtual ~AtomLib();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ATOMLIB_H_ */
