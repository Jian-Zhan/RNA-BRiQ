/*
 * BaseRotamerLib.h
 *
 *  Created on: Jul 22, 2019
 *      Author: s2982206
 */

#ifndef MODEL_BASEROTAMERLIB_H_
#define MODEL_BASEROTAMERLIB_H_

#include <time.h>
#include "model/BaseRotamer.h"
#include "para/Parameter.h"

namespace NSPmodel {

using namespace NSPpara;

class BaseRotamerLib {

public:
	BaseRotamer* rotLib[4][1500];

	BaseRotamerLib();


	BaseRotamer* getLowestEnergyRotamer(int type){
		if(type == 0)
			return rotLib[0][568];
		else if(type == 1)
			return rotLib[1][90];
		else if(type == 2)
			return rotLib[2][37];
		else
			return rotLib[3][41];
	}

	BaseRotamer* getRandomRotamerLv1(int baseType){
		return rotLib[baseType][rand()%1500];
	}

	BaseRotamer* getFlipRotamer(int baseType){
		/*
		 * type0: imp < 0 && chi > -55 && chi < 160    600
		 * type1: imp < 0 && (chi < -55 || chi > 160)  300
		 * type2: imp > 0 && chi > -55 && chi < 160    300
		 * type3: imp > 0 && (chi < -55 || chi > 160)  300
		 */

		int r1 = rand()%10;
		if(r1 < 2)
			return rotLib[baseType][rand()%600];
		else if(r1 < 6)
			return rotLib[baseType][rand()%300 + 600];
		else if(r1 < 7)
			return rotLib[baseType][rand()%300 + 900];
		else
			return rotLib[baseType][rand()%300 + 1200];

	}


	BaseRotamer* getNearestRotamer(RNABase* base){
		int type = base->baseTypeInt;
		double minDist = 9999.9;

		if(!base->backboneComplete())
			return getLowestEnergyRotamer(base->baseTypeInt);

		BaseRotamer* rot = new BaseRotamer(base);
		BaseRotamer* nearest;
		int id = 0;
		for(int i=0;i<1500;i++){
			double d = rotLib[type][i]->distanceTo(rot);
			if(d < minDist){
				minDist = d;
				id = i;
				nearest = rotLib[type][i];
			}
		}

		delete rot;
		return nearest;
	}

	virtual ~BaseRotamerLib();
};

} /* namespace NSPforcefield */

#endif /* MODEL_BASEROTAMERLIB_H_ */
