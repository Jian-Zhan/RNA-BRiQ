/*
 * DMSearcher.h
 *
 *  Created on: Apr 13, 2019
 *      Author: s2982206
 */

#ifndef FORCEFIELD_DMSEARCHER_H_
#define FORCEFIELD_DMSEARCHER_H_

#include "model/BaseDistanceMatrix.h"
#include "dataio/datapaths.h"
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>

namespace NSPforcefield {

using namespace std;
using namespace NSPmodel;

class DMSearcher {
public:
	vector<BaseDistanceMatrix*> dmList;
	int x;
	int y;
	int z;
	int indexTable[2744][1000];

	DMSearcher(const vector<BaseDistanceMatrix*>& list) {
		int n = list.size();
		double cutoff = 1.0;
		x = 0;
		y = 1;
		z = 2;
		int maxElementSize = 1000;


		int a,b,c,i,j,k,l;
		int ix,iy,iz;
		double d, dx, dy, dz, dxMin, dxMax, dyMin, dyMax, dzMin, dzMax;

		for(i=0;i<14;i++){
			for(j=0;j<maxElementSize;j++){
				this->indexTable[i][j] = 0;
			}
		}

		for(i=0;i<n;i++) {
			this->dmList.push_back(list[i]);
		}

		int minLocal = n;
		int max;
		for(a=0;a<9;a++) {
			for(b=a+1;b<9;b++) {
				for(c=b+1;c<9;c++){
					int max = 0;
					int count[14][14][14];
					for(i=0;i<14;i++){
						for(j=0;j<14;j++){
							for(k=0;k<14;k++) {
								count[i][j][k] = 0;
							}
						}
					}
					for(i=0;i<n;i++){
						d = dmList[i]->dm[a];
						ix = (int)d - 3;
						if(ix < 0) ix = 0;
						if(ix > 13) ix = 13;

						d = dmList[i]->dm[b];
						iy = (int)d - 3;
						if(iy < 0) iy = 0;
						if(iy > 13) iy = 13;

						d = dmList[i]->dm[c];
						iz = (int)d - 3;
						if(iz < 0) iz = 0;
						if(iz > 13) iz = 13;

						count[ix][iy][iz] ++;
					}
					for(i=0;i<14;i++) {
						for(int j=0;j<14;j++){
							for(int k=0;k<14;k++){
								if(count[i][j][k] > max) max = count[i][j][k];
							}
						}
					}
					if(max < minLocal) {
						minLocal = max;
						x = a;
						y = b;
						z = c;
					}
				}
			}
		}

		cout << x << " " << y << " " << z << endl;
		cout << "minLocal " << minLocal << endl;

		double minX,minY,minZ,maxX,maxY,maxZ;
		int currentIndex;

		for(l=0;l<n;l++) {
			for(i=0;i<14;i++){
				minX = i+2.0;
				dxMin = dmList[l]->dm[x] - minX;
				if(dmList[l]->dm[x] < 3.0 && i==0) dxMin = 0;
				if(minX < dmList[l]->dm[x] - 5.0) continue;
				maxX = i+3.0;
				dxMax = dmList[l]->dm[x] - maxX;
				if(dmList[l]->dm[x] > 15.0 && i==13) dxMax = 0;
				if(maxX > dmList[l]->dm[x] + 5.0 && i>0) break;
				for(j=0;j<14;j++){
					minY = j+2.0;
					dyMin = dmList[l]->dm[y] - minY;
					if(dmList[l]->dm[y] < 3.0 && j==0) dyMin = 0;
					if(minY < dmList[l]->dm[y]-5.0) continue;
					maxY = j+3.0;
					dyMax = dmList[l]->dm[y] - maxY;
					if(dmList[l]->dm[y] > 15.0 && j==13) dyMax = 0;
					if(maxY > dmList[l]->dm[y] + 5.0 && j>0) break;

					for(k=0;k<14;k++){
						minZ = k+2.0;
						dzMin = dmList[l]->dm[z] - minZ;
						if(dmList[l]->dm[z] < 3.0 && k==0) dzMin = 0;
						if(minZ < dmList[l]->dm[z]-5.0) continue;
						maxZ = k+3.0;
						dzMax = dmList[l]->dm[z] - maxZ;
						if(dmList[l]->dm[z] > 15.0 && k==13) dzMax = 0;
						if(maxZ > dmList[l]->dm[z] + 5.0 && k>0) break;

						if(dxMin*dxMin + dyMin*dyMin + dzMin*dzMin < cutoff ||
						   dxMin*dxMin + dyMin*dyMin + dzMax*dzMax < cutoff ||
						   dxMin*dxMin + dyMax*dyMax + dzMin*dzMin < cutoff ||
						   dxMin*dxMin + dyMax*dyMax + dzMax*dzMax < cutoff ||
					       dxMax*dxMax + dyMin*dyMin + dzMin*dzMin < cutoff ||
						   dxMax*dxMax + dyMin*dyMin + dzMax*dzMax < cutoff ||
						   dxMax*dxMax + dyMax*dyMax + dzMin*dzMin < cutoff ||
						   dxMax*dxMax + dyMax*dyMax + dzMax*dzMax < cutoff) {
							currentIndex = this->indexTable[i*196+j*14+k][0];
							if(currentIndex+ 2 > maxElementSize) {
								continue;
								cerr << "max element size not enough: " << maxElementSize << endl;
							}
							this->indexTable[i*196+j*14+k][0] ++;
							this->indexTable[i*196+j*14+k][currentIndex+1] = l;
						}
					}
				}
			}
		}

		int maxN = 0;
		int maxNIndex = 0;
		for(i=0;i<2744;i++){
			if(this->indexTable[i][0] > 0)
			{
				ix = i/196;
				iy = (i-ix*14)/14;
				iz = i%14;
				//cout << ix << " " << iy << " " << iz << " " << this->indexTable[i][0] << endl;
				if(this->indexTable[i][0] > maxN) {
					maxN = this->indexTable[i][0];
					maxNIndex = i;
				}
			}
		}
		cout << "maxN: " << maxN << endl;
		ix = maxNIndex/196;
		iy = (maxNIndex-ix*14)/14;
		iz = maxNIndex%14;
		cout << ix << " " << iy << " " << iz << " " << this->indexTable[maxNIndex][0] << endl;
	}



	int getNearestIndex(BaseDistanceMatrix* dm){
		return 0;
	}

	virtual ~DMSearcher();
};



} /* namespace NSPforcefield */

#endif /* FORCEFIELD_DMSEARCHER_H_ */
