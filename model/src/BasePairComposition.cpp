/*
 * BasePairComposition.cpp
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#include <model/BasePairComposition.h>

namespace NSPmodel {

BasePairComposition::BasePairComposition() {
	for(int i=0;i<16;i++){
		counts[i] = 0.0;
		pa[i] = 0.0;
		totalCount = 0.0;
	}
	bpMap["AA"] = 0;
	bpMap["AU"] = 1;
	bpMap["AG"] = 2;
	bpMap["AC"] = 3;
	bpMap["UA"] = 4;
	bpMap["UU"] = 5;
	bpMap["UG"] = 6;
	bpMap["UC"] = 7;
	bpMap["GA"] = 8;
	bpMap["GU"] = 9;
	bpMap["GG"] = 10;
	bpMap["GC"] = 11;
	bpMap["CA"] = 12;
	bpMap["CU"] = 13;
	bpMap["CG"] = 14;
	bpMap["CC"] = 15;
}

BasePairComposition::BasePairComposition(const string& line){
	vector<string> list;
	NSPtools::splitString(line, " ", &list);
	for(int i=0;i<16;i++){
		counts[i] = 0.0;
		pa[i] = atof(list[i].c_str());
		totalCount = 0.0;
	}
	bpMap["AA"] = 0;
	bpMap["AU"] = 1;
	bpMap["AG"] = 2;
	bpMap["AC"] = 3;
	bpMap["UA"] = 4;
	bpMap["UU"] = 5;
	bpMap["UG"] = 6;
	bpMap["UC"] = 7;
	bpMap["GA"] = 8;
	bpMap["GU"] = 9;
	bpMap["GG"] = 10;
	bpMap["GC"] = 11;
	bpMap["CA"] = 12;
	bpMap["CU"] = 13;
	bpMap["CG"] = 14;
	bpMap["CC"] = 15;
}

BasePairComposition::~BasePairComposition() {
	// TODO Auto-generated destructor stub
}

SingleBaseComposition::SingleBaseComposition() {
	for(int i=0;i<4;i++){
		counts[i] = 0.0;
		pa[i] = 0.0;
		totalCount = 0.0;
	}
	bpMap['A'] = 0;
	bpMap['U'] = 1;
	bpMap['G'] = 2;
	bpMap['C'] = 3;
}

SingleBaseComposition::SingleBaseComposition(const string& line){
	vector<string> list;
	NSPtools::splitString(line, " ", &list);
	for(int i=0;i<4;i++){
		counts[i] = 0.0;
		pa[i] = atof(list[i].c_str());
		totalCount = 0.0;
	}
	bpMap['A'] = 0;
	bpMap['U'] = 1;
	bpMap['G'] = 2;
	bpMap['C'] = 3;
}

SingleBaseComposition::~SingleBaseComposition() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPtest */
