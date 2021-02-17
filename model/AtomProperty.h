/*
 * AtomProperty.h
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ATOMPROPERTY_H_
#define DESIGNSEQ_ATOMPROPERTY_H_

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

namespace NSPmodel {
using namespace std;

class AtomProperty {
public:
	string atomUniqueName;
	string atomName;
	int connectNum;
	float vdwRadius;
	float sp2Radius[91];
	float cosSquareToRadius[101];
	float polarRadius;
	bool isHDonor;
	bool isHAcceptor;
	bool isPolar;
	bool isAromatic;

	AtomProperty(const string& line);
	float getSp2Radius(float angleToNorm);
	float getSp2Radius2(float cosSquare)
	{
		int i = (int)(cosSquare*100);
		if(i<0 || i > 100)
		{
			cerr << "invalid cosSquare: " << cosSquare << endl;
		}
		return cosSquareToRadius[i];
	}
	virtual ~AtomProperty();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ATOMPROPERTY_H_ */
