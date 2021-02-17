/*
 * BasePair.cpp
 *
 *  Created on: 2020Äê11ÔÂ17ÈÕ
 *      Author: pengx
 */

#include "model/BasePair.h"

namespace NSPmodel {

BasePair::BasePair(RNABase* baseA, RNABase* baseB, RnaAtomLib* atLib) {
	// TODO Auto-generated constructor stub

	if(!baseA->sidechainComplete(atLib) || !baseB->sidechainComplete(atLib))
	{
		cout << "sidechain not complete" << endl;
		return;
	}

	this->baseA = baseA;
	this->baseB = baseB;
	this->dm = BaseDistanceMatrix(*baseA, *baseB);
	this->hbNum = 0;
	char xx[20];
	sprintf(xx, "%c%c", baseA->baseType, baseB->baseType);
	this->type = string(xx);


	vector<string>* baseAtomListA = atLib->getAtomNames(baseA->baseTypeInt);
	vector<string>* baseAtomListB = atLib->getAtomNames(baseB->baseTypeInt);

	for(int i=0;i<baseAtomListA->size();i++){
		Atom* a = baseA->getAtom(baseAtomListA->at(i));
		if(a == NULL){
			continue;
		}
		if(a->name[0] == 'C') continue;
		PolarAtom pa(baseA, a->name);

		if(pa.isEmpty()) continue;

		for(int j=0;j<baseAtomListB->size();j++){

			Atom* b = baseB->getAtom(baseAtomListB->at(j));
			if(b == NULL){
				continue;
			}
			if(b->name[0] == 'C') continue;

			PolarAtom pb(baseB, b->name);
			if(pa.isEmpty()) continue;
			if(pa.hbondedTo(&pb))
				this->hbNum++;
		}
	}

	if(hbNum>0) this->isHbondPair = true;
	else this->isHbondPair = false;

}

bool BasePair::isWCPair(){
	char typeA = baseA->baseType;
	char typeB = baseB->baseType;
	if((typeA == 'A' && typeB == 'U') || (typeA == 'U' && typeB == 'A') || (typeA == 'G' && typeB == 'C') || (typeA == 'C' && typeB == 'G')) {
		BaseDistanceMatrix dm0("G C  3.7980  6.9500  9.2390  4.6079  7.1241 12.2823 13.0212  8.3017  8.7787 12.5239 10.7388  5.6006  3.8695  7.9905  6.0148  3.5253");
		if(dm.distanceTo(dm0) < 1.5)
			return true;
		else
			return false;
	}
	else if(typeA == 'G' && typeB == 'U') {
		BaseDistanceMatrix dm0("G U  4.5095  6.9273 10.4250  6.2696  7.1630 12.0360 13.7964  9.2787  7.4334 11.5636 10.6688  5.2313  1.8353  6.8238  6.6813  4.1656");
		if(dm.distanceTo(dm0) < 1.5)
			return true;
		else
			return false;
	}
	else if(typeA == 'U' && typeB == 'G') {
		BaseDistanceMatrix dm0("U G  4.5095  7.1630  7.4334  1.8353  6.9273 12.0360 11.5636  6.8238 10.4250 13.7964 10.6688  6.6813  6.2696  9.2787  5.2313  4.1656");
		if(dm.distanceTo(dm0) < 1.5)
			return true;
		else
			return false;
	}
	return false;
}

bool BasePair::isHbondedPair(){
	if(this->hbNum > 0) return true;
	else return false;
}

double BasePair::distanceToWCPair(){
	char typeA = baseA->baseType;
	char typeB = baseB->baseType;

	if(typeA == 'G' && typeB == 'U') {
		BaseDistanceMatrix dm0("G U  4.5095  6.9273 10.4250  6.2696  7.1630 12.0360 13.7964  9.2787  7.4334 11.5636 10.6688  5.2313  1.8353  6.8238  6.6813  4.1656");
		return dm.distanceTo(dm0);
	}
	else if(typeA == 'U' && typeB == 'G') {
		BaseDistanceMatrix dm0("U G  4.5095  7.1630  7.4334  1.8353  6.9273 12.0360 11.5636  6.8238 10.4250 13.7964 10.6688  6.6813  6.2696  9.2787  5.2313  4.1656");
		return dm.distanceTo(dm0);
	}
	else {
		BaseDistanceMatrix dm0("G C  3.7980  6.9500  9.2390  4.6079  7.1241 12.2823 13.0212  8.3017  8.7787 12.5239 10.7388  5.6006  3.8695  7.9905  6.0148  3.5253");
		return dm.distanceTo(dm0);
	}
}

BasePair::~BasePair() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
