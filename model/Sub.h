/*
 * Sub.h
 *
 *  Created on: Nov 14, 2018
 *      Author: s2982206
 */

#ifndef MODEL_SUB_H_
#define MODEL_SUB_H_
#include "geometry/localframe.h"
#include "geometry/CsMove.h"
#include "geometry/xyz.h"

namespace NSPmodel {
using namespace NSPgeometry;

class Sub {
public:
	int seqID;
	LocalFrame cs;
	XYZ tPho;
	char type;
	int fatherIndex;
	bool isRoot;
	int childIndex[3];

	Sub();
	Sub(LocalFrame& cs, char type, int seqID);
	Sub& operator=(const Sub& other){
		this->seqID = other.seqID;
		this->cs = other.cs;
		this->tPho = other.tPho;
		this->type = other.type;
		this->fatherIndex = other.fatherIndex;
		this->isRoot = other.isRoot;
		this->childIndex[0] = other.childIndex[0];
		this->childIndex[1] = other.childIndex[1];
		this->childIndex[2] = other.childIndex[2];
		return *this;
	}
	Sub applyMove(CsMove& move);
	Sub copy();

	virtual ~Sub();
};

} /* namespace NSPtest */

#endif /* MODEL_SUB_H_ */
