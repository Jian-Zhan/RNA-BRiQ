/*
 * RNABaseLib.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: s2982206
 */

#include <model/RNABaseLib.h>

namespace NSPmodel {

RNABaseLib::RNABaseLib() {
	// TODO Auto-generated constructor stub

}

RNABase* RNABaseLib::getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs){
	string augc = "AUGC";
	RNABase* base = new RNABase(baseID, chainID, augc[baseType]);
    vector<Atom*> list;
    vector<string>* names = atLib.getSidechainAtoms(baseType);
    vector<XYZ> tList;
    if(baseType == 0){
            XYZ a = XYZ(1.468,   0.000,   0.000); //A-N9
            XYZ b = XYZ(2.306,  -1.084,   0.000); //A-C8
            XYZ c = XYZ(3.577,  -0.770,   0.000); //A-N7
            XYZ d = XYZ(3.576,   0.615,   0.000); //A-C5
            XYZ e = XYZ(4.614,   1.556,   0.000); //A-C6
            XYZ f = XYZ(5.904,   1.230,   0.000); //A-N6
            XYZ g = XYZ(4.276,   2.861,   0.000); //A-N1
            XYZ h = XYZ(2.980,   3.184,   0.000); //A-C2
            XYZ i = XYZ(1.914,   2.391,   0.000); //A-N3
            XYZ j = XYZ(2.285,   1.103,   0.000); //A-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
    }
    else if(baseType == 1){
            XYZ a = XYZ(1.478,   0.000,   0.000); //U-N1
            XYZ b = XYZ(2.122,   1.221,   0.000); //U-C2
            XYZ c = XYZ(1.528,   2.282,   0.000); //U-O2
            XYZ d = XYZ(3.491,   1.159,   0.000); //U-N3
            XYZ e = XYZ(4.265,   0.020,   0.000); //U-C4
            XYZ f = XYZ(5.490,   0.123,   0.000); //U-O4
            XYZ g = XYZ(3.526,  -1.204,   0.000); //U-C5
            XYZ h = XYZ(2.191,  -1.173,   0.000); //U-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }
    else if(baseType == 2) {
            XYZ a = XYZ(1.468,   0.000,   0.000); //G-N9
            XYZ b = XYZ(2.295,  -1.094,   0.000); //G-C8
            XYZ c = XYZ(3.560,  -0.779,   0.000); //G-N7
            XYZ d = XYZ(3.570,   0.606,   0.000); //G-C5
            XYZ e = XYZ(4.655,   1.514,   0.000); //G-C6
            XYZ f = XYZ(5.864,   1.265,   0.000); //G-O6
            XYZ g = XYZ(4.221,   2.832,   0.000); //G-N1
            XYZ h = XYZ(2.909,   3.225,   0.000); //G-C2
            XYZ i = XYZ(2.690,   4.543,   0.000); //G-N2
            XYZ j = XYZ(1.886,   2.389,   0.000); //G-N3
            XYZ k = XYZ(2.287,   1.103,   0.000); //G-C4
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
            tList.push_back(i);
            tList.push_back(j);
            tList.push_back(k);
    }
    else if(baseType == 3) {
            XYZ a = XYZ(1.478,   0.000,   0.000); //C-N1
            XYZ b = XYZ(2.151,   1.224,   0.000); //C-C2
            XYZ c = XYZ(1.490,   2.271,   0.000); //C-O2
            XYZ d = XYZ(3.503,   1.239,   0.000); //C-N3
            XYZ e = XYZ(4.178,   0.091,   0.000); //C-C4
            XYZ f = XYZ(5.508,   0.150,   0.000); //C-N4
            XYZ g = XYZ(3.519,  -1.170,   0.000); //C-C5
            XYZ h = XYZ(2.181,  -1.170,   0.000); //C-C6
            tList.push_back(a);
            tList.push_back(b);
            tList.push_back(c);
            tList.push_back(d);
            tList.push_back(e);
            tList.push_back(f);
            tList.push_back(g);
            tList.push_back(h);
    }

    for(int i=0;i<tList.size();i++){
    	base->addAtom(new Atom(names->at(i), local2global(cs, tList[i])));
    }
    return base;
}

RNABase* RNABaseLib::getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, BaseRotamer* rot){
	RNABase* base = getBase(baseID, chainID, baseType, cs);
    vector<string> names;
	names.push_back("C1'");
    names.push_back("C2'");
    names.push_back("C3'");
    names.push_back("C4'");
    names.push_back("O4'");
    names.push_back("O2'");
    names.push_back("O3'");
    names.push_back("C5'");

    for(int i=0;i<8;i++){
    	base->addAtom(new Atom(names[i], local2global(cs, rot->tList1[i])));
    }

	return base;
}

RNABase* RNABaseLib::getBase(const string& baseID, const string& chainID, int baseType, LocalFrame& cs, BaseRotamer* rot, PhophateGroupLocal* pl){
	RNABase* base = getBase(baseID, chainID, baseType, cs);
    vector<string> names;
	names.push_back("C1'");
    names.push_back("C2'");
    names.push_back("C3'");
    names.push_back("C4'");
    names.push_back("O4'");
    names.push_back("O2'");
    names.push_back("O3'");
    names.push_back("C5'");

    for(int i=0;i<8;i++){
    	base->addAtom(new Atom(names[i], local2global(cs, rot->tList1[i])));
    }
    LocalFrame cs2 = cs + rot->mv12;
    PhophateGroup p = PhophateGroup(*pl, cs2);
	base->addAtom(new Atom("P", p.tList[0]));
	base->addAtom(new Atom("O5'", p.tList[1]));
	base->addAtom(new Atom("OP1", p.tList[2]));
	base->addAtom(new Atom("OP2", p.tList[3]));
	return base;
}

RNABase* RNABaseLib::toStandardBase(RNABase* base){
	vector<string> names;
	vector<XYZ> tListStd;
	vector<XYZ> tListOld;
	names.push_back("C1'");
	if(base->baseTypeInt == 0){
		XYZ o = XYZ(0.0, 0.0, 0.0); //C1'
        XYZ a = XYZ(1.468,   0.000,   0.000); //A-N9
        XYZ b = XYZ(2.306,  -1.084,   0.000); //A-C8
        XYZ c = XYZ(3.577,  -0.770,   0.000); //A-N7
        XYZ d = XYZ(3.576,   0.615,   0.000); //A-C5
        XYZ e = XYZ(4.614,   1.556,   0.000); //A-C6
        XYZ f = XYZ(5.904,   1.230,   0.000); //A-N6
        XYZ g = XYZ(4.276,   2.861,   0.000); //A-N1
        XYZ h = XYZ(2.980,   3.184,   0.000); //A-C2
        XYZ i = XYZ(1.914,   2.391,   0.000); //A-N3
        XYZ j = XYZ(2.285,   1.103,   0.000); //A-C4
        tListStd.push_back(o);
        tListStd.push_back(a);
        tListStd.push_back(b);
        tListStd.push_back(c);
        tListStd.push_back(d);
        tListStd.push_back(e);
        tListStd.push_back(f);
        tListStd.push_back(g);
        tListStd.push_back(h);
        tListStd.push_back(i);
        tListStd.push_back(j);

        names.push_back("N9");
        names.push_back("C8");
        names.push_back("N7");
        names.push_back("C5");
        names.push_back("C6");
        names.push_back("N6");
        names.push_back("N1");
        names.push_back("C2");
        names.push_back("N3");
        names.push_back("C4");
	}
	else if(base->baseTypeInt == 1){
		XYZ o = XYZ(0.0, 0.0, 0.0); //C1'
        XYZ a = XYZ(1.478,   0.000,   0.000); //U-N1
        XYZ b = XYZ(2.122,   1.221,   0.000); //U-C2
        XYZ c = XYZ(1.528,   2.282,   0.000); //U-O2
        XYZ d = XYZ(3.491,   1.159,   0.000); //U-N3
        XYZ e = XYZ(4.265,   0.020,   0.000); //U-C4
        XYZ f = XYZ(5.490,   0.123,   0.000); //U-O4
        XYZ g = XYZ(3.526,  -1.204,   0.000); //U-C5
        XYZ h = XYZ(2.191,  -1.173,   0.000); //U-C6
        tListStd.push_back(o);
        tListStd.push_back(a);
        tListStd.push_back(b);
        tListStd.push_back(c);
        tListStd.push_back(d);
        tListStd.push_back(e);
        tListStd.push_back(f);
        tListStd.push_back(g);
        tListStd.push_back(h);

        names.push_back("N1");
        names.push_back("C2");
        names.push_back("O2");
        names.push_back("N3");
        names.push_back("C4");
        names.push_back("O4");
        names.push_back("C5");
        names.push_back("C6");
	}
	else if(base->baseTypeInt == 2){
		XYZ o = XYZ(0.0, 0.0, 0.0); //C1'
        XYZ a = XYZ(1.468,   0.000,   0.000); //G-N9
        XYZ b = XYZ(2.295,  -1.094,   0.000); //G-C8
        XYZ c = XYZ(3.560,  -0.779,   0.000); //G-N7
        XYZ d = XYZ(3.570,   0.606,   0.000); //G-C5
        XYZ e = XYZ(4.655,   1.514,   0.000); //G-C6
        XYZ f = XYZ(5.864,   1.265,   0.000); //G-O6
        XYZ g = XYZ(4.221,   2.832,   0.000); //G-N1
        XYZ h = XYZ(2.909,   3.225,   0.000); //G-C2
        XYZ i = XYZ(2.690,   4.543,   0.000); //G-N2
        XYZ j = XYZ(1.886,   2.389,   0.000); //G-N3
        XYZ k = XYZ(2.287,   1.103,   0.000); //G-C4
        tListStd.push_back(o);
        tListStd.push_back(a);
        tListStd.push_back(b);
        tListStd.push_back(c);
        tListStd.push_back(d);
        tListStd.push_back(e);
        tListStd.push_back(f);
        tListStd.push_back(g);
        tListStd.push_back(h);
        tListStd.push_back(i);
        tListStd.push_back(j);
        tListStd.push_back(k);
        names.push_back("N9");
        names.push_back("C8");
        names.push_back("N7");
        names.push_back("C5");
        names.push_back("C6");
        names.push_back("O6");
        names.push_back("N1");
        names.push_back("C2");
        names.push_back("N2");
        names.push_back("N3");
        names.push_back("C4");
	}
	else if(base->baseTypeInt == 3){
		XYZ o = XYZ(0.0, 0.0, 0.0); //C1'
        XYZ a = XYZ(1.478,   0.000,   0.000); //C-N1
        XYZ b = XYZ(2.151,   1.224,   0.000); //C-C2
        XYZ c = XYZ(1.490,   2.271,   0.000); //C-O2
        XYZ d = XYZ(3.503,   1.239,   0.000); //C-N3
        XYZ e = XYZ(4.178,   0.091,   0.000); //C-C4
        XYZ f = XYZ(5.508,   0.150,   0.000); //C-N4
        XYZ g = XYZ(3.519,  -1.170,   0.000); //C-C5
        XYZ h = XYZ(2.181,  -1.170,   0.000); //C-C6
        tListStd.push_back(o);
        tListStd.push_back(a);
        tListStd.push_back(b);
        tListStd.push_back(c);
        tListStd.push_back(d);
        tListStd.push_back(e);
        tListStd.push_back(f);
        tListStd.push_back(g);
        tListStd.push_back(h);
        names.push_back("N1");
        names.push_back("C2");
        names.push_back("O2");
        names.push_back("N3");
        names.push_back("C4");
        names.push_back("N4");
        names.push_back("C5");
        names.push_back("C6");
	}

	for(int i=0;i<names.size();i++){
		Atom* at = base->getAtom(names[i]);
		if(at == NULL){
			RNABase* baseB = new RNABase(base->baseID, base->chainID, base->baseType);
			vector<Atom*>* atomList = base->getAtomList();
			for(int j=0;j<atomList->size();j++){
				Atom* a = atomList->at(j);
				baseB->addAtom(new Atom(a->name, a->coord));
			}
			baseB->baseID = base->baseID;
			baseB->baseSeqID = base->baseSeqID;
			baseB->chainID = base->chainID;
			baseB->baseTypeInt = base->baseTypeInt;
			baseB->baseType = base->baseType;
			baseB->hasLocalFrame = base->hasLocalFrame;
			baseB->coordSys = base->coordSys;
			return baseB;
		}
		tListOld.push_back(at->coord);
	}

	XYZ Acog = NSPgeometry::getCOG(tListOld);
	XYZ Bcog = NSPgeometry::getCOG(tListStd);

	int len = tListStd.size();
	if(tListStd.size() != tListOld.size()){
		cout << "atom number not equal: " << tListStd.size() << " " << tListOld.size() << " " << base->baseType << endl;
		exit(1);
	}

	vector<XYZ> listA;
	vector<XYZ> listB;
	for(int i=0;i<len;i++){
		XYZ a = tListOld[i] - Acog;
		XYZ b = tListStd[i] - Bcog;
		listA.push_back(a);
		listB.push_back(b);
	}


	TransForm tf =  NSPgeometry::buildRotation(listA, listB);

	vector<XYZ> newListB;
	for(int i=0;i<len;i++) {
		XYZ c = tf.transform(listB[i]);
		newListB.push_back(c);
	}

	vector<XYZ> newListB2;
	for(int i=0;i<len;i++){
		XYZ c = newListB[i] + Acog;
		newListB2.push_back(c);
	}

	RNABase* baseB = new RNABase(base->baseID, base->chainID, base->baseType);
	vector<Atom*>* atomList = base->getBackboneAtoms();
	for(int j=0;j<atomList->size();j++){
		Atom* a = atomList->at(j);
		if(a->name == "C1'") continue;
		baseB->addAtom(new Atom(a->name, a->coord));
	}

	for(int i=0;i<names.size();i++){
		baseB->addAtom(new Atom(names[i], newListB2[i]));
	}
	baseB->baseID = base->baseID;
	baseB->baseSeqID = base->baseSeqID;
	baseB->chainID = base->chainID;
	baseB->baseTypeInt = base->baseTypeInt;
	baseB->baseType = base->baseType;
	baseB->hasLocalFrame = base->hasLocalFrame;
	baseB->coordSys = base->coordSys;
	return baseB;
}

RNABaseLib::~RNABaseLib() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPmodel */
