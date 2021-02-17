/*
 * RiboConnectHashMap.h
 *
 *  Created on: Jul 23, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_RIBOCONNECTHASHMAP_H_
#define FORCEFIELD_RIBOCONNECTHASHMAP_H_

#include <vector>
#include <map>
#include <fstream>
#include "geometry/CsMove.h"
#include "dataio/datapaths.h"
#include "tools/StringTool.h"
#include "model/PhophateGroup.h"

namespace NSPforcefield {

using namespace std;
using namespace NSPgeometry;
using namespace NSPmodel;


class RiboConnectHashMap {
public:
	map<string,int> keyToIndex;
	vector<XYZ> PList;
	vector<XYZ> O5List;
	vector<XYZ> OP1List;
	vector<XYZ> OP2List;
	vector<double> eneList;

	int repNum;
	vector<string> repKeys;
	vector<double> repEnes;
	vector<XYZ> repPList;
	vector<XYZ> repO5List;
	vector<XYZ> repOP1List;
	vector<XYZ> repOP2List;

	double connectionEnergyFactor = 1.0;

	RiboConnectHashMap();
	void setWeight(double wt);
	void setConnectionEnergyFactor(double f){
		this->connectionEnergyFactor = f;
	}

	PhophateGroupLocal findPhophate(const string& key, double d0){
		map<string,int>::iterator it;
		it = keyToIndex.find(key);
		if(it != keyToIndex.end()) {
			int index = it->second;
			return PhophateGroupLocal(PList[index], O5List[index], OP1List[index], OP2List[index], eneList[index]);
		}

		int minD = 999.9;
		int minIndex = 0;
		int d, dd, i, j;

		for(i=0;i<repNum;i++){
			dd = 0;
			for(j=0;j<9;j++){
				d = key[j] - repKeys[i][j];
				dd += d*d;
			}
			if(dd < minD){
				minD = dd;
				minIndex = i;
			}
		}

		double e = repEnes[minIndex] + dd*0.1*connectionEnergyFactor;
		if(d0 < 2.6)
			e = e + (2.6-d0)*10*connectionEnergyFactor;
		else if(d0 > 4.0)
			e = e + (d0-4.0)*10*connectionEnergyFactor;


		return PhophateGroupLocal(repPList[minIndex], repO5List[minIndex], repOP1List[minIndex], repOP2List[minIndex], e);
	}


	virtual ~RiboConnectHashMap();
};

} /* namespace NSPforcefield */

#endif /* FORCEFIELD_RIBOCONNECTHASHMAP_H_ */
