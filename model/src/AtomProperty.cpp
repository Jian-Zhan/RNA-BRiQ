/*
 * AtomProperty.cpp
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#include "model/AtomProperty.h"
#include "tools/StringTool.h"

namespace NSPmodel {

using namespace NSPtools;

AtomProperty::AtomProperty(const string& line) {
	vector<string> spt;
	string s;
	splitString(line, " ", &spt);
	this->atomUniqueName = spt.at(0);
	this->atomName = spt.at(1);
	this->vdwRadius = atof(spt.at(2).c_str());
	this->connectNum = atoi(spt.at(3).c_str());
	s = spt.at(4);
	if(s == "d")
	{
		this->isHDonor = true;
		this->isHAcceptor = false;
		this->isPolar = true;
	}
	else if(s == "a")
	{
		this->isHDonor = false;
		this->isHAcceptor = true;
		this->isPolar = true;
	}
	else if(s == "b")
	{
		this->isHDonor = true;
		this->isHAcceptor = true;
		this->isPolar = true;
	}
	else
	{
		this->isHDonor = false;
		this->isHAcceptor = false;
		this->isPolar = false;
	}

	s = spt.at(5);
	if(s == "a")
	{
		this->isAromatic = true;

		float r0 = atof(spt.at(6).c_str());
		float r1 = atof(spt.at(7).c_str());

        for(int i=0;i<91;i++)
        {
                float u = (90-i)/52.0;
                sp2Radius[i] = r0+(r1-r0)*exp(-u*u);
        }


		for(int i=0;i<101;i++)
		{
			float cosAng = sqrt(i*0.01);
			float ang = acos(cosAng)*180/3.14159265;
			float u = (90-ang)/52.0;
			cosSquareToRadius[i] = r0+(r1-r0)*exp(-u*u);
		}


	}
	else
	{
		this->isAromatic = false;
		for(int i=0;i<101;i++)
		{
			cosSquareToRadius[i]  = this->vdwRadius;
		}
        for(int i=0;i<91;i++)
        {
                sp2Radius[i] = this->vdwRadius;
        }

		//this->vdwRadiusM = vdwRadius;
		//this->vdwRadiusL = vdwRadius;
	}

	if(isHDonor || isHAcceptor)
	{
		/*
		 * to be modified
		 */
		char type = atomName.at(0);
		if(type == 'O')
			this->polarRadius = 1.35;
		else
			this->polarRadius = 1.45;
	}
	else
		this->polarRadius = vdwRadius;



	// TODO Auto-generated constructor stub

}




float AtomProperty::getSp2Radius(float angleToNorm)
{
	int a = (int)angleToNorm;
	if(a > 90)
		a = 180 - a;
	if(a < 0 || a >90)
	{
		cerr << "invalid angle: " << angleToNorm << endl;
		exit(1);
	}
	return sp2Radius[a];
}

AtomProperty::~AtomProperty() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
