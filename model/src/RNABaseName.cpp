/*
 * RNABaseName.cpp
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#include "model/RNABaseName.h"

namespace NSPmodel {

RNABaseName::RNABaseName() {
	// TODO Auto-generated constructor stub
	this->rnaSeq = "AUGCN";
	this->sinToIntMap['A'] = 0;
	this->sinToIntMap['U'] = 1;
	this->sinToIntMap['G'] = 2;
	this->sinToIntMap['C'] = 3;
	this->sinToIntMap['T'] = 1;
	this->sinToIntMap['a'] = 0;
	this->sinToIntMap['u'] = 1;
	this->sinToIntMap['g'] = 2;
	this->sinToIntMap['t'] = 1;
	this->sinToIntMap['c'] = 3;
	this->sinToIntMap['N'] = 4;
	this->sinToIntMap['X'] = 4;
}

char RNABaseName::intToSin(int i) const{
	if(i>=0 && i <=4)
		return rnaSeq[i];
	return 'N';
}

int RNABaseName::sinToInt(char sin) const{
	map<char,int>::const_iterator it = sinToIntMap.find(sin);
	if(it != sinToIntMap.end())
		return it->second;
	return 4;
}

RNABaseName::~RNABaseName() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
