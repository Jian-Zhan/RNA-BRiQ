/*
 * FragmentLibrary.cpp
 *
 *  Created on: Nov 15, 2019
 *      Author: s2982206
 */

#include <pred/FragmentLibrary.h>

namespace NSPpred {

RiboConnectLib::RiboConnectLib(){
	string path = NSPdataio::datapath();
	string file = path + "fragLib/riboMove.txt";
	ifstream f;
	f.open(file.c_str(), ios::in);
	string line;
	while(getline(f, line)){
		CsMove cm(line);
		this->cmList.push_back(cm);
	}
	this->num = cmList.size();
}

CsMove RiboConnectLib::getRandomMove(){
	return cmList[rand()%this->num];
}

RiboConnectLib::~RiboConnectLib(){

}


F2Fragment::F2Fragment(int typeA, int typeB){
	this->typeA = typeA;
	this->typeB = typeB;
	this->isReverse = false;
	this->hasRibose = false;
	this->baseFlip = false;
	this->level1ID = 0;
	this->level2ID = 0;
	this->mvType = "";
	cm.oriMove = XYZ(4.0, 4.0, 4.0);
	this->cm.oriMove = XYZ(4.0, 4.0, 4.0);
}

F2Fragment::F2Fragment(const string& line,int typeA, int typeB, BaseRotamerLib* rotLib){

	this->typeA = typeA;
	this->typeB = typeB;
	this->isReverse = false;
	this->hasRibose = false;
	this->baseFlip = false;
	this->level1ID = 0;
	this->level2ID = 0;
	this->mvType = "";
	this->cm = CsMove(line.substr(4, 131));
}

double F2Fragment::distanceTo(F2Fragment* other){
	BaseDistanceMatrix dm1(this->cm);
	BaseDistanceMatrix dm2(other->cm);
	return dm1.distanceTo(dm2);
}

void F2Fragment::printPDB(const string& outfile){
	RNABaseLib bLib;
	RNAChain* rc = new RNAChain("A");
	LocalFrame csA;
	LocalFrame csB = csA + cm;
	rc->addBase(bLib.getBase("1", "A", typeA, csA));
	rc->addBase(bLib.getBase("2", "A", typeB, csB));


	ofstream out;
	out.open(outfile.c_str(), ios::out);
	rc->printPDBFormat(out, 1);
	out.close();
	delete rc;
}

F2Fragment::~F2Fragment(){

}


F3Fragment::F3Fragment(const string& line,int typeA, int typeB, int typeC, BaseRotamerLib* rotLib){
	this->cmAB = CsMove(line.substr(0, 131));
	this->cmBC = CsMove(line.substr(132, 131));
	this->typeA = typeA;
	this->typeB = typeB;
	this->typeC = typeC;
	vector<string> spt;
	splitString(line, " ", &spt);
	int a = atoi(spt[24].c_str());
	int b = atoi(spt[25].c_str());
	int c = atoi(spt[26].c_str());

	this->rotA = rotLib->rotLib[typeA][a];
	this->rotB = rotLib->rotLib[typeB][b];
	this->rotC = rotLib->rotLib[typeC][c];
	this->plA = PhophateGroupLocal(line.substr(277, 101));
	this->plB = PhophateGroupLocal(line.substr(381, 101));
}

void F3Fragment::printPDB(const string& outfile){
	RNABaseLib bLib;
	RNAChain* rc = new RNAChain("A");
	LocalFrame csA;
	LocalFrame csB = csA + cmAB;
	LocalFrame csC = csB + cmBC;
	rc->addBase(bLib.getBase("1", "A", typeA, csA, rotA, &plA));
	rc->addBase(bLib.getBase("2", "A", typeB, csB, rotB, &plB));
	rc->addBase(bLib.getBase("3", "A", typeC, csC, rotC));
	ofstream out;
	out.open(outfile.c_str(), ios::out);
	rc->printPDBFormat(out, 1);
	out.close();
	delete rc;

}

F3Fragment::~F3Fragment(){

}

F2FragmentLib::F2FragmentLib(const string& tag, BaseRotamerLib* rotLib){
	string path = NSPdataio::datapath();
	this->hasRibose = false;
	this->isReverse = false;

	if(tag == "AG"){
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "GA") {
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = false;
	}

	string libPath = path + "fragLib/"+tag +"/";
	ifstream file;
	string line;
	int level1ID = 0;


	file.open(libPath + "level1/" +tag+".rot", ios::in);
	if(!file.is_open()){
		cout << "fail to open: " << libPath + "level1/" +tag+".rot" << endl;
		exit(1);
	}

	int typeA, typeB;
	if(tag == "AG"){
		typeA = 0;
		typeB = 2;
	}
	else{
		typeA = 2;
		typeB = 0;
	}

	level1ID = 0;
	getline(file,line);
	while(getline(file, line)){
		F2Fragment* frag = new F2Fragment(line, typeA, typeB, rotLib);
		double minD = 999.9;

		frag->mvType = tag;

		frag->level1ID = level1ID;
		frag->level2ID = 0;
		level1ID ++;
		frag->hasRibose = this->hasRibose;
		frag->isReverse = this->isReverse;
		this->fragListLevel1.push_back(frag);
	}
	file.close();

	char yy[20];
	for(int i=0;i<rotNum;i++){
		sprintf(yy, "%d", i);
		file.open(libPath + "level2/" +tag+"-" + string(yy)+".rot", ios::in);
		if(!file.is_open()){
			cout << "fail to open: " << libPath + "level2/" +tag+"-" + string(yy)+".rot" << endl;
			exit(1);
		}
		int level2ID = 0;
		getline(file, line);
		while(getline(file, line)){
			F2Fragment* frag = new F2Fragment(line, typeA, typeB, rotLib);

			double minD = 999.9;

			frag->mvType = tag;

			frag->level1ID = i;
			frag->level2ID = level2ID;
			level2ID++;
			frag->hasRibose = this->hasRibose;
			frag->isReverse = this->isReverse;
			this->fragListLevel2.push_back(frag);
		}
		file.close();
	}
}

F2FragmentLib::F2FragmentLib(const string& tag, int typeA, int typeB, BaseRotamerLib* rotLib){
	string path = NSPdataio::datapath();
	this->hasRibose = false;
	this->isReverse = false;

	if(tag == "wc") {
		this->rotNum = 60;
		this->hasRibose = false;
		this->isReverse = false;
	}
	else if(tag == "nwc"){
		this->rotNum = 60;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "wcNb"){
		this->rotNum = 60;
		this->hasRibose = false;
		this->isReverse = false;
	}
	else if(tag == "nwcNb"){
		this->rotNum = 60;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "revNb"){
		this->rotNum = 60;
		this->hasRibose = true;
		this->isReverse = true;
	}
	else if(tag == "loopNb"){
		this->rotNum = 60;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "bulge13"){
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "bulge14"){
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = false;
	}
	else if(tag == "revBulge13"){
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = true;
	}
	else if(tag == "revBulge14"){
		this->rotNum = 20;
		this->hasRibose = true;
		this->isReverse = true;
	}
	else if(tag == "jump"){
		this->rotNum = 1;
		this->hasRibose = false;
		this->isReverse = false;
		F2Fragment* frag0 = new F2Fragment(typeA, typeB);
		frag0->mvType = "A";
		F2Fragment* frag1 = new F2Fragment(typeA, typeB);
		frag1->mvType = "A";
		F2Fragment* frag2 = new F2Fragment(typeA, typeB);
		frag2->mvType = "A";
		this->fragListLevel0.push_back(frag0);
		this->fragListLevel1.push_back(frag1);
		this->fragListLevel2.push_back(frag2);
		return;
	}

	string libPath = path + "fragLib/"+tag +"/";
	ifstream file;
	char xx[20];
	char yy[20];
	string line;
	int level1ID = 0;

	sprintf(xx, "%d", typeA*4+typeB);


	file.open(libPath + "level0/" + tag + ".rot-"+string(xx), ios::in);
	if(!file.is_open()){
		cout << "fail to open: " << libPath + "level0/" +tag+".rot-"+string(xx) << endl;
		exit(1);
	}

	getline(file, line);
	while(getline(file, line)){
		F2Fragment* frag = new F2Fragment(line, typeA, typeB, rotLib);
		frag->level1ID = 0;
		frag->level2ID = 0;
		frag->hasRibose = this->hasRibose;
		frag->isReverse = this->isReverse;
		this->fragListLevel0.push_back(frag);
	}

	int level0RotNum = this->fragListLevel0.size();
	file.close();


	file.open(libPath + "level1/" +tag+".rot-"+string(xx), ios::in);
	if(!file.is_open()){
		cout << "fail to open: " << libPath + "level1/" +tag+".rot-"+string(xx) << endl;
		exit(1);
	}

	level1ID = 0;
	getline(file,line);
	while(getline(file, line)){
		F2Fragment* frag = new F2Fragment(line, typeA, typeB, rotLib);
		double minD = 999.9;
		int lv0ID = 0;
		for(int i=0;i<level0RotNum;i++){
			double d = frag->distanceTo(this->fragListLevel0[i]);
			if(d < minD){
				minD = d;
				lv0ID = i;
			}
		}

		char mt[2];
		mt[0] = 'A'+lv0ID;
		mt[1] = '\0';
		frag->mvType = string(mt);

		frag->level1ID = level1ID;
		frag->level2ID = 0;
		level1ID ++;
		frag->hasRibose = this->hasRibose;
		frag->isReverse = this->isReverse;
		this->fragListLevel1.push_back(frag);
	}
	file.close();

	for(int i=0;i<rotNum;i++){
		sprintf(xx, "%d", typeA*4+typeB);
		sprintf(yy, "%d", i);
		file.open(libPath + "level2/" +tag+".rot-"+string(xx) + "-" + string(yy), ios::in);
		if(!file.is_open()){
			cout << "fail to open: " << libPath + "level2/" +tag+".rot-"+string(xx) + "-" + string(yy) << endl;
			exit(1);
		}
		int level2ID = 0;
		getline(file, line);
		while(getline(file, line)){
			F2Fragment* frag = new F2Fragment(line, typeA, typeB, rotLib);
			if(i > 43 && (tag == "loopNb" || tag == "revNb"))
				frag->baseFlip = true;
			double minD = 999.9;
			int lv0ID = 0;
			for(int i=0;i<level0RotNum;i++){
				double d = frag->distanceTo(this->fragListLevel0[i]);
				if(d < minD){
					minD = d;
					lv0ID = i;
				}
			}
			char mt[2];
			mt[0] = 'A'+lv0ID;
			mt[1] = '\0';
			frag->mvType = string(mt);

			frag->level1ID = i;
			frag->level2ID = level2ID;
			level2ID++;
			frag->hasRibose = this->hasRibose;
			frag->isReverse = this->isReverse;
			this->fragListLevel2.push_back(frag);
		}
		file.close();
	}
}

F2FragmentLib::~F2FragmentLib(){
	for(int i=0;i<fragListLevel0.size();i++){
		delete fragListLevel0[i];
	}

	for(int i=0;i<fragListLevel1.size();i++){
		delete fragListLevel1[i];
	}

	for(int i=0;i<fragListLevel2.size();i++){
		delete fragListLevel2[i];
	}
}

F3FragmentLib::F3FragmentLib(const string& tag, int typeA, int typeB, int typeC, BaseRotamerLib* rotLib){
	string path = NSPdataio::datapath();
	string libPath = path + "fragLib/"+tag +"/";
	if(tag == "acF3")
		rotNum = 50;
	else if(tag == "allF3")
		rotNum = 100;
	ifstream file;
	string fileName;
	char xx[20];
	sprintf(xx, "%d", typeA*16+typeB*4+typeC);
	if(tag == "acF3")
		fileName = libPath+"f3-ac-"+string(xx)+".frag";
	else if(tag == "allF3")
		fileName = libPath+"f3-all-"+string(xx)+".frag";
	file.open(fileName, ios::in);
	if(!file.is_open()){
		cout << "fail to open: " << fileName << endl;
		exit(1);
	}
	string line;
	while(getline(file, line)){
		F3Fragment* frag = new F3Fragment(line, typeA, typeB, typeC, rotLib);
		this->f3List.push_back(frag);
	}
	file.close();
}

F3FragmentLib::~F3FragmentLib(){
	for(int i=0;i<f3List.size();i++){
		delete f3List[i];
	}
}

FragmentLibrary::FragmentLibrary(BaseRotamerLib* rotLib) {
	// TODO Auto-generated constructor stub
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			wcNb.push_back(new F2FragmentLib("wcNb", i, j, rotLib));
			nwcNb.push_back(new F2FragmentLib("nwcNb", i, j, rotLib));
			loopNb.push_back(new F2FragmentLib("loopNb", i, j, rotLib));

			revNb.push_back(new F2FragmentLib("revNb", i, j, rotLib));
			wcPair.push_back(new F2FragmentLib("wc", i, j, rotLib));
			nwcPair.push_back(new F2FragmentLib("nwc", i, j, rotLib));

			bulge13.push_back(new F2FragmentLib("bulge13", i, j, rotLib));

			jumpLib.push_back(new F2FragmentLib("jump", i, j, rotLib));
			for(int k=0;k<4;k++){
				acF3Lib.push_back(new F3FragmentLib("acF3", i, j, k, rotLib));
				allF3Lib.push_back(new F3FragmentLib("allF3", i, j, k, rotLib));
			}
		}
	}

	this->agLib = new F2FragmentLib("AG", rotLib);
	this->gaLib = new F2FragmentLib("GA", rotLib);

	this->ribLib = new RiboConnectLib();
}

FragmentLibrary::~FragmentLibrary() {
	// TODO Auto-generated destructor stub

	for(int i=0;i<16;i++){
		delete wcNb[i];
		delete nwcNb[i];
		delete loopNb[i];
		delete revNb[i];
		delete wcPair[i];
		delete nwcPair[i];
		delete bulge13[i];
		delete jumpLib[i];
	}

	for(int i=0;i<64;i++){
		delete acF3Lib[i];
		delete allF3Lib[i];
	}

	delete ribLib;
}

} /* namespace NSPforcefield */
