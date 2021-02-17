/*
 * BasePairComposition.h
 *
 *  Created on: Nov 15, 2018
 *      Author: s2982206
 */

#ifndef MODEL_BASEPAIRCOMPOSITION_H_
#define MODEL_BASEPAIRCOMPOSITION_H_

#include <map>
#include "tools/StringTool.h"
namespace NSPmodel {

using namespace std;

class BasePairComposition {
public:
	double counts[16];
	double pa[16];
	double totalCount;
	map<string, int> bpMap;
	BasePairComposition();
	BasePairComposition(const string& line);
	void addCount(const string& key, double wt){
		counts[bpMap[key]] += wt;
		totalCount += wt;
	}
	void updateDistribution() {
		if(totalCount == 0) return;
		double bg[] = {0.0584, 0.0581, 0.0515, 0.0369, 0.0577, 0.0259, 0.0544, 0.0334, 0.0444, 0.0554, 0.0923, 0.1514, 0.0342, 0.0361, 0.1513, 0.0586};
		for(int i=0;i<16;i++) {
			pa[i] = (counts[i]+bg[i]*1)/(totalCount+1);
		}
	}
	double getP(const string& key){
		return pa[bpMap[key]];
	}

	virtual ~BasePairComposition();
};

class SingleBaseComposition {
public:
	double counts[4];
	double pa[4];
	double totalCount;
	map<char, int> bpMap;
	SingleBaseComposition();
	SingleBaseComposition(const string& line);
	void addCount(char key, double wt){
		counts[bpMap[key]] += wt;
		totalCount += wt;
	}
	void updateDistribution() {
		if(totalCount == 0) return;
		double bg[] = {0.25, 0.25, 0.25, 0.25};
		for(int i=0;i<4;i++) {
			pa[i] = (counts[i]+bg[i]*1)/(totalCount+1);
		}
	}
	double getP(char key){
		return pa[bpMap[key]];
	}

	virtual ~SingleBaseComposition();
};

} /* namespace NSPtest */

#endif /* MODEL_BASEPAIRCOMPOSITION_H_ */
