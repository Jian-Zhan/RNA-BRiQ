/*
 * ResName.cpp
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#include "model/ResName.h"

namespace NSPmodel {

ResName::ResName() {
	// TODO Auto-generated constructor stub
	aaSeq = "ACDEFGHIKLMNPQRSTVWYBJX";

	this->triNameList.reserve(23);
	this->triNameList.push_back("ALA");
	this->triNameList.push_back("CYS");
	this->triNameList.push_back("ASP");
	this->triNameList.push_back("GLU");
	this->triNameList.push_back("PHE");
	this->triNameList.push_back("GLY");
	this->triNameList.push_back("HIS");
	this->triNameList.push_back("ILE");
	this->triNameList.push_back("LYS");
	this->triNameList.push_back("LEU");
	this->triNameList.push_back("MET");
	this->triNameList.push_back("ASN");
	this->triNameList.push_back("PRO");
	this->triNameList.push_back("GLN");
	this->triNameList.push_back("ARG");
	this->triNameList.push_back("SER");
	this->triNameList.push_back("THR");
	this->triNameList.push_back("VAL");
	this->triNameList.push_back("TRP");
	this->triNameList.push_back("TYR");
	this->triNameList.push_back("MSE");
	this->triNameList.push_back("PSD");
	this->triNameList.push_back("UNK");

	this->triToIntMap["ALA"] = 0;
	this->triToIntMap["CYS"] = 1;
	this->triToIntMap["ASP"] = 2;
	this->triToIntMap["GLU"] = 3;
	this->triToIntMap["PHE"] = 4;
	this->triToIntMap["GLY"] = 5;
	this->triToIntMap["HIS"] = 6;
	this->triToIntMap["ILE"] = 7;
	this->triToIntMap["LYS"] = 8;
	this->triToIntMap["LEU"] = 9;
	this->triToIntMap["MET"] = 10;
	this->triToIntMap["ASN"] = 11;
	this->triToIntMap["PRO"] = 12;
	this->triToIntMap["GLN"] = 13;
	this->triToIntMap["ARG"] = 14;
	this->triToIntMap["SER"] = 15;
	this->triToIntMap["THR"] = 16;
	this->triToIntMap["VAL"] = 17;
	this->triToIntMap["TRP"] = 18;
	this->triToIntMap["TYR"] = 19;
	this->triToIntMap["MSE"] = 20;
	this->triToIntMap["PSD"] = 21;
	this->triToIntMap["UNK"] = 22;

	this->sinToIntMap['A'] = 0;
	this->sinToIntMap['C'] = 1;
	this->sinToIntMap['D'] = 2;
	this->sinToIntMap['E'] = 3;
	this->sinToIntMap['F'] = 4;
	this->sinToIntMap['G'] = 5;
	this->sinToIntMap['H'] = 6;
	this->sinToIntMap['I'] = 7;
	this->sinToIntMap['K'] = 8;
	this->sinToIntMap['L'] = 9;
	this->sinToIntMap['M'] = 10;
	this->sinToIntMap['N'] = 11;
	this->sinToIntMap['P'] = 12;
	this->sinToIntMap['Q'] = 13;
	this->sinToIntMap['R'] = 14;
	this->sinToIntMap['S'] = 15;
	this->sinToIntMap['T'] = 16;
	this->sinToIntMap['V'] = 17;
	this->sinToIntMap['W'] = 18;
	this->sinToIntMap['Y'] = 19;
	this->sinToIntMap['B'] = 20;
	this->sinToIntMap['J'] = 21;
	this->sinToIntMap['X'] = 22;

	this->rnaNameToInt["0A"] = 0;
	this->rnaNameToInt["1MA"] = 0;
	this->rnaNameToInt["2AD"] = 0;
	this->rnaNameToInt["2BA"] = 0;
	this->rnaNameToInt["2MA"] = 0;
	this->rnaNameToInt["31H"] = 0;
	this->rnaNameToInt["31M"] = 0;
	this->rnaNameToInt["365"] = 0;
	this->rnaNameToInt["3AD"] = 0;
	this->rnaNameToInt["3AT"] = 0;
	this->rnaNameToInt["3DA"] = 0;
	this->rnaNameToInt["45A"] = 0;
	this->rnaNameToInt["4BW"] = 0;
	this->rnaNameToInt["574"] = 0;
	this->rnaNameToInt["5AA"] = 0;
	this->rnaNameToInt["5AD"] = 0;
	this->rnaNameToInt["6HA"] = 0;
	this->rnaNameToInt["6IA"] = 0;
	this->rnaNameToInt["6MZ"] = 0;
	this->rnaNameToInt["6NW"] = 0;
	this->rnaNameToInt["84T"] = 0;
	this->rnaNameToInt["8AN"] = 0;
	this->rnaNameToInt["a"] = 0;
	this->rnaNameToInt["A"] = 0;
	this->rnaNameToInt["  A"] = 0;
	this->rnaNameToInt["  U"] = 1;
	this->rnaNameToInt["  G"] = 2;
	this->rnaNameToInt["  C"] = 3;
	this->rnaNameToInt["A23"] = 0;
	this->rnaNameToInt["A2M"] = 0;
	this->rnaNameToInt["A2P"] = 0;
	this->rnaNameToInt["A44"] = 0;
	this->rnaNameToInt["A6A"] = 0;
	this->rnaNameToInt["A9Z"] = 0;
	this->rnaNameToInt["ACP"] = 0;
	this->rnaNameToInt["ADN"] = 0;
	this->rnaNameToInt["ADP"] = 0;
	this->rnaNameToInt["ADS"] = 0;
	this->rnaNameToInt["AET"] = 0;
	this->rnaNameToInt["AF2"] = 0;
	this->rnaNameToInt["AMO"] = 0;
	this->rnaNameToInt["AMP"] = 0;
	this->rnaNameToInt["ANP"] = 0;
	this->rnaNameToInt["APC"] = 0;
	this->rnaNameToInt["AT7"] = 0;
	this->rnaNameToInt["ATP"] = 0;
	this->rnaNameToInt["AVC"] = 0;
	this->rnaNameToInt["DA"] = 0;
	this->rnaNameToInt["DJF"] = 0;
	this->rnaNameToInt["DTP"] = 0;
	this->rnaNameToInt["EEM"] = 0;
	this->rnaNameToInt["F3N"] = 0;
	this->rnaNameToInt["F3O"] = 0;
	this->rnaNameToInt["GAP"] = 0;
	this->rnaNameToInt["GOM"] = 0;
	this->rnaNameToInt["GSU"] = 0;
	this->rnaNameToInt["ILA"] = 0;
	this->rnaNameToInt["LCA"] = 0;
	this->rnaNameToInt["LMS"] = 0;
	this->rnaNameToInt["M3O"] = 0;
	this->rnaNameToInt["MA6"] = 0;
	this->rnaNameToInt["MIA"] = 0;
	this->rnaNameToInt["MSP"] = 0;
	this->rnaNameToInt["N6G"] = 0;
	this->rnaNameToInt["N79"] = 0;
	this->rnaNameToInt["PPU"] = 0;
	this->rnaNameToInt["QSI"] = 0;
	this->rnaNameToInt["SAH"] = 0;
	this->rnaNameToInt["SAM"] = 0;
	this->rnaNameToInt["SFG"] = 0;
	this->rnaNameToInt["T6A"] = 0;
	this->rnaNameToInt["TSB"] = 0;
	this->rnaNameToInt["VAA"] = 0;
	this->rnaNameToInt["XAR"] = 0;
	this->rnaNameToInt["YMP"] = 0;
	this->rnaNameToInt["0U"] = 1;
	this->rnaNameToInt["0U1"] = 1;
	this->rnaNameToInt["2AU"] = 1;
	this->rnaNameToInt["2KH"] = 1;
	this->rnaNameToInt["2MU"] = 1;
	this->rnaNameToInt["5BU"] = 1;
	this->rnaNameToInt["5FU"] = 1;
	this->rnaNameToInt["5GS"] = 1;
	this->rnaNameToInt["5IU"] = 1;
	this->rnaNameToInt["5MU"] = 1;
	this->rnaNameToInt["PSU"] = 1;
	this->rnaNameToInt["6FU"] = 1;
	this->rnaNameToInt["6GS"] = 1;
	this->rnaNameToInt["6HT"] = 1;
	this->rnaNameToInt["6OP"] = 1;
	this->rnaNameToInt["75B"] = 1;
	this->rnaNameToInt["A6U"] = 1;
	this->rnaNameToInt["BRU"] = 1;
	this->rnaNameToInt["CM0"] = 1;
	this->rnaNameToInt["DT"] = 1;
	this->rnaNameToInt["DU"] = 1;
	this->rnaNameToInt["DUT"] = 1;
	this->rnaNameToInt["F2T"] = 1;
	this->rnaNameToInt["H2U"] = 1;
	this->rnaNameToInt["IU"] = 1;
	this->rnaNameToInt["4SU"] = 1;
	this->rnaNameToInt["NTT"] = 1;
	this->rnaNameToInt["OMU"] = 1;
	this->rnaNameToInt["SSU"] = 1;
	this->rnaNameToInt["T2T"] = 1;
	this->rnaNameToInt["T5S"] = 1;
	this->rnaNameToInt["TLN"] = 1;
	this->rnaNameToInt["TM2"] = 1;
	this->rnaNameToInt["u"] = 1;
	this->rnaNameToInt["U"] = 1;
	this->rnaNameToInt["U33"] = 1;
	this->rnaNameToInt["U34"] = 1;
	this->rnaNameToInt["U36"] = 1;
	this->rnaNameToInt["U37"] = 1;
	this->rnaNameToInt["U5M"] = 1;
	this->rnaNameToInt["U5R"] = 1;
	this->rnaNameToInt["UBD"] = 1;
	this->rnaNameToInt["UD5"] = 1;
	this->rnaNameToInt["UDP"] = 1;
	this->rnaNameToInt["UFT"] = 1;
	this->rnaNameToInt["UMO"] = 1;
	this->rnaNameToInt["UMS"] = 1;
	this->rnaNameToInt["UR3"] = 1;
	this->rnaNameToInt["URI"] = 1;
	this->rnaNameToInt["URU"] = 1;
	this->rnaNameToInt["UTP"] = 1;
	this->rnaNameToInt["UVP"] = 1;
	this->rnaNameToInt["UZR"] = 1;
	this->rnaNameToInt["XTR"] = 1;
	this->rnaNameToInt["0G"] = 2;
	this->rnaNameToInt["1MG"] = 2;
	this->rnaNameToInt["23G"] = 2;
	this->rnaNameToInt["2MG"] = 2;
	this->rnaNameToInt["2SG"] = 2;
	this->rnaNameToInt["5CG"] = 2;
	this->rnaNameToInt["5GP"] = 2;
	this->rnaNameToInt["6HG"] = 2;
	this->rnaNameToInt["7MG"] = 2;
	this->rnaNameToInt["A6G"] = 2;
	this->rnaNameToInt["BGM"] = 2;
	this->rnaNameToInt["C2E"] = 2;
	this->rnaNameToInt["DG"] = 2;
	this->rnaNameToInt["DGP"] = 2;
	this->rnaNameToInt["g"] = 2;
	this->rnaNameToInt["G"] = 2;
	this->rnaNameToInt["G46"] = 2;
	this->rnaNameToInt["G48"] = 2;
	this->rnaNameToInt["G7M"] = 2;
	this->rnaNameToInt["GCP"] = 2;
	this->rnaNameToInt["GDO"] = 2;
	this->rnaNameToInt["GDP"] = 2;
	this->rnaNameToInt["GF2"] = 2;
	this->rnaNameToInt["GH3"] = 2;
	this->rnaNameToInt["GMP"] = 2;
	this->rnaNameToInt["GNG"] = 2;
	this->rnaNameToInt["GNP"] = 2;
	this->rnaNameToInt["GRB"] = 2;
	this->rnaNameToInt["GTP"] = 2;
	this->rnaNameToInt["QUO"] = 2;
	this->rnaNameToInt["LCG"] = 2;
	this->rnaNameToInt["M2G"] = 2;
	this->rnaNameToInt["M7G"] = 2;
	this->rnaNameToInt["MGT"] = 2;
	this->rnaNameToInt["OMG"] = 2;
	this->rnaNameToInt["PGN"] = 2;
	this->rnaNameToInt["PGP"] = 2;
	this->rnaNameToInt["XGR"] = 2;
	this->rnaNameToInt["XUG"] = 2;
	this->rnaNameToInt["YG"] = 2;
	this->rnaNameToInt["YYG"] = 2;
	this->rnaNameToInt["0C"] = 3;
	this->rnaNameToInt["1SC"] = 3;
	this->rnaNameToInt["4OC"] = 3;
	this->rnaNameToInt["5CF"] = 3;
	this->rnaNameToInt["5CM"] = 3;
	this->rnaNameToInt["5HC"] = 3;
	this->rnaNameToInt["5HM"] = 3;
	this->rnaNameToInt["5IC"] = 3;
	this->rnaNameToInt["5MC"] = 3;
	this->rnaNameToInt["6FC"] = 3;
	this->rnaNameToInt["6HC"] = 3;
	this->rnaNameToInt["6OO"] = 3;
	this->rnaNameToInt["73W"] = 3;
	this->rnaNameToInt["A5M"] = 3;
	this->rnaNameToInt["A6C"] = 3;
	this->rnaNameToInt["c"] = 3;
	this->rnaNameToInt["C"] = 3;
	this->rnaNameToInt["C43"] = 3;
	this->rnaNameToInt["C5P"] = 3;
	this->rnaNameToInt["CBR"] = 3;
	this->rnaNameToInt["CBV"] = 3;
	this->rnaNameToInt["CCC"] = 3;
	this->rnaNameToInt["CDP"] = 3;
	this->rnaNameToInt["CFL"] = 3;
	this->rnaNameToInt["CFZ"] = 3;
	this->rnaNameToInt["CH1"] = 3;
	this->rnaNameToInt["CSG"] = 3;
	this->rnaNameToInt["CSL"] = 3;
	this->rnaNameToInt["CTP"] = 3;
	this->rnaNameToInt["DC"] = 3;
	this->rnaNameToInt["DCP"] = 3;
	this->rnaNameToInt["DCT"] = 3;
	this->rnaNameToInt["DCZ"] = 3;
	this->rnaNameToInt["DOC"] = 3;
	this->rnaNameToInt["LCC"] = 3;
	this->rnaNameToInt["N5C"] = 3;
	this->rnaNameToInt["N5M"] = 3;
	this->rnaNameToInt["NCU"] = 3;
	this->rnaNameToInt["O2C"] = 3;
	this->rnaNameToInt["OMC"] = 3;
	this->rnaNameToInt["RPC"] = 3;
	this->rnaNameToInt["RSQ"] = 3;
	this->rnaNameToInt["S4C"] = 3;
	this->rnaNameToInt["XCR"] = 3;

	augc.push_back("A");
	augc.push_back("U");
	augc.push_back("G");
	augc.push_back("C");
}

char ResName::intToSin(int i) const{
	if(i>=0 && i<=22)
			return aaSeq.at(i);
		return 'X';
}

int ResName::sinToInt(char sin) const{
	map<char,int>::const_iterator it = sinToIntMap.find(sin);
	if(it != sinToIntMap.end())
		return it->second;
	else
		return 22;
}

int ResName::triToInt(const string& tri) const{
	map<string,int>::const_iterator it = triToIntMap.find(tri);
	if(it != triToIntMap.end())
		return it->second;

	if(tri == "HID" || tri == "HIE" || tri == "HIP")
		return 6;
	return 22;
}

string ResName::intToTri(int i) const{
	if(i>=0 && i<=22)
		return triNameList.at(i);
	else
		return "UNK";
}

char ResName::triToSin(const string& tri) const{
	return intToSin(triToInt(tri));
}

string ResName::sinToTri(char sin) const{
	return intToTri(sinToInt(sin));
}

bool ResName::isStandardAminoAcid(const string& tri){
	return triToInt(tri) < 20;
}

bool ResName::isAminoAcid(const string& tri){
    int aaInt = triToInt(tri);
    if(aaInt < 22) return true;
    else if(tri == "GLX" || tri == "LEF" || tri == "PCA" || tri == "OCS") return true;
    else return false;

}

int ResName::chiNum(const string& tri) const
{
	int chiNum[] = {0, 1, 2, 3, 2, 0, 2, 2, 4, 2, 3, 2, 1, 3, 4, 1, 1, 1, 2, 2};
	int type = triToInt(tri);
	if(type > 19)
		return 0;
	else
		return chiNum[type];

}

bool ResName::isRNABase(const string& baseName) const{
	map<string,int>::const_iterator it = this->rnaNameToInt.find(baseName);
	if(it != rnaNameToInt.end())
		return true;
	else
		return false;
}

string ResName::toStandardBase(const string& baseName) const{
	map<string,int>::const_iterator it = this->rnaNameToInt.find(baseName);

	if(it != rnaNameToInt.end())
		return this->augc[it->second];
	else
		return "UNK";
}


ResName::~ResName() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
