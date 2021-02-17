/*
 * datapaths.cpp
 *
 *      Author: hyliu
 */

#include "dataio/datapaths.h"
#include <cstdlib>
std::string NSPdataio::getenvpath(const std::string & envvar) {
	char *envpath=std::getenv(envvar.c_str());
	std::string res="";
	if(envpath) res=std::string(envpath);
//	char sep= boost::filesystem::path::preferred_separator;
        char sep='/';
	if(res !="" && res.back() != sep) res +=sep;
	return res;
}

std::string NSPdataio::datapath(){
	return getenvpath("BRiQ_DATAPATH");

}

