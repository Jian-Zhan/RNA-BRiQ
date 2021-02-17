/*
 * CsMoveTo6DKey.cpp
 *
 *  Created on: Sep 6, 2019
 *      Author: s2982206
 */

#include <pred/CsMoveTo6DKey.h>

namespace NSPpred {

CsMoveTo6DKey::CsMoveTo6DKey() {
	// TODO Auto-generated constructor stub

	ifstream file;
	string fileName = NSPdataio::datapath() + "sphere500.index";
	file.open(fileName, ios::in);
	if(!file.is_open()) {
		cout << "can't open file: " << fileName << endl;
	}
	int keyIndex;
	int spIndex;
	while(file >> keyIndex >> spIndex){
		this->sphereKeyMap[keyIndex] = spIndex;
	}
	file.close();

}

string CsMoveTo6DKey::toKey(const LocalFrame& csA, const LocalFrame& csB){
	double len = csA.origin_.distance(csB.origin_);
	double xLen = 1.0/len;

	XYZ tBInA = global2local(csA, csB.origin_)*xLen;
	XYZ tAInB = global2local(csB, csA.origin_)*(-xLen);

	char ss[7];
	ss[6] = '\0';

	int idX1 = (int)((tBInA.x_ + 1.0)*100);
	if(idX1 == 200) idX1 = 199;
	int idY1 = (int)((tBInA.y_ + 1.0)*100);
	if(idY1 == 200) idY1 = 199;
	int idZ1 = (int)((tBInA.z_ + 1.0)*100);
	if(idZ1 == 200) idZ1 = 199;

	int idX2 = (int)((tAInB.x_ + 1.0)*100);
	if(idX2 == 200) idX2 = 199;
	int idY2 = (int)((tAInB.y_ + 1.0)*100);
	if(idY2 == 200) idY2 = 199;
	int idZ2 = (int)((tAInB.z_ + 1.0)*100);
	if(idZ2 == 200) idZ2 = 199;

	int spIndexA = sphereKeyMap[idX1*40000+idY1*200+idZ1];
	int spIndexB = sphereKeyMap[idX2*40000+idY2*200+idZ2];

	XYZ a = csA.x * (-1.0);
	XYZ b = (csB.origin_ - csA.origin_)*xLen;
	XYZ c = csB.x;
	XYZ n1 = a^b;
	XYZ n2 = b^c;
	XYZ m = n1^b;

	double ang = atan2(m*n2, n1*n2);
	//double ang = 0;
	if(ang < 0)
		ang = -ang * 57.29577952;
	else
		ang = 360 - ang * 57.29577952;


	int indexAng = (int)(ang*0.125 + 0.5);
	if(indexAng == 45) indexAng = 0;

	ss[0] = spIndexA/30 + '!';
	ss[1] = spIndexA%30 + '!';
	ss[2] = spIndexB/30 + '!';
	ss[3] = spIndexB%30 + '!';
	ss[4] = (int)(len*3.33333333333+0.5) + '!';
	ss[5] = indexAng + '!';
	return string(ss);
}

CsMoveTo6DKey::~CsMoveTo6DKey() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPforcefield */
