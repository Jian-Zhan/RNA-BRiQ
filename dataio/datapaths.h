/*
 * datapaths.h
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */

#ifndef DATAIO_DATAPATHS_H_
#define DATAIO_DATAPATHS_H_
#include <string>
namespace NSPdataio{
/*!get path string from a system ENVIRONMENTAL variable
 * returns an empty string if envvar is not a defined environmental variable name.
 */
std::string getenvpath(const std::string & envvar);

std::string datapath();

}

#endif /* DATAIO_DATAPATHS_H_ */
