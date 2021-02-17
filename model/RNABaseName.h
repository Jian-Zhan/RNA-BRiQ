/*
 * RNABaseName.h
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#ifndef MODEL_RNABASENAME_H_
#define MODEL_RNABASENAME_H_

#include <string>
#include <vector>
#include <map>

namespace NSPmodel {

using namespace std;
class RNABaseName {
private:
	string rnaSeq;
	map<char,int> sinToIntMap;
public:
	RNABaseName();
	char intToSin(int i) const;
	int sinToInt(char sin) const;
	virtual ~RNABaseName();
};

} /* namespace NSPmodel */

#endif /* MODEL_RNABASENAME_H_ */
