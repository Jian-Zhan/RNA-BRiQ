/*
 * InputParser.h
 *
 *  Created on: Feb 22, 2019
 *      Author: s2982206
 */

#ifndef TOOLS_INPUTPARSER_H_
#define TOOLS_INPUTPARSER_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "tools/StringTool.h"

namespace NSPtools {

using namespace std;

class InputOption{
public:
	string optionTag;
	string optionValue;

	InputOption(const string& tag, const string& value){
		this->optionTag = tag;
		this->optionValue = value;
	}

	virtual ~InputOption();
};

class InputParser {
private:
	vector<InputOption> options;
public:
	InputParser();
	InputParser(const string file){
		ifstream f;
		f.open(file.c_str(), ios::in);
		if(!f.is_open()){
			cout << "fail to open file: " << file << endl;
			exit(1);
		}
		string s;
		while(getline(f,s)){
			vector<string> spt;
			splitString(s, " ", &spt);
			if(spt.size() < 2) continue;
			string name = spt[0];
			string value = s.substr(name.length(), s.length()-name.length());
			value = trimString(value);
			options.push_back(InputOption(name,value));
		}
	}

	bool specifiedOption(const string& tag){
		for(int i=0;i<options.size();i++){
			if(options.at(i).optionTag.compare(tag) == 0)
				return true;
		}
		return false;
	}

	string getValue(const string& tag){
		for(int i=0;i<options.size();i++){
			//cout << options[i].optionTag << " " << options[i].optionValue << " " << tag << endl;
			if(options.at(i).optionTag.compare(tag) == 0) {
				return options.at(i).optionValue;
			}
		}
		return "";
	}

	void printOptions() {
		for(int i=0;i<options.size();i++) {
			cout << options[i].optionTag << " " << options[i].optionValue << endl;
		}
	}

	virtual ~InputParser();
};

} /* namespace NSPtools */

#endif /* TOOLS_INPUTPARSER_H_ */
