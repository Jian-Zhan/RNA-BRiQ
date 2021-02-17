/*
 * BaseConnection.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: s2982206
 */

#include "pred/BaseConnection.h"

namespace NSPpred {



BaseConnection::~BaseConnection() {
	// TODO Auto-generated destructor stub
	if(this->mutator != NULL)
		delete mutator;
	if(this->childOrNotChild != NULL)
		delete [] childOrNotChild;

}

} /* namespace NSPpred */
