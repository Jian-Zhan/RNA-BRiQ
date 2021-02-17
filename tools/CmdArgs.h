/*
 * CmdArgs.h
 *
 *  Created on: Feb 12, 2019
 *      Author: s2982206
 */

#ifndef TOOLS_CMDARGS_H_
#define TOOLS_CMDARGS_H_

#include <string>
#include <vector>
#include <iostream>

namespace NSPtools {

using namespace std;

class CmdOption{
public:
	string optionTag;
	string optionValue;

	CmdOption(string& tag, string& value){
		this->optionTag = tag;
		this->optionValue = value;
	}

	virtual ~CmdOption();
};

class CmdArgs {

private:
	vector<CmdOption> options;
public:
	CmdArgs(int argc, char** args){
		for(int i=1;i<argc;i++){
			if(args[i][0] == '-'){
				string tag = string(args[i]);
				string value = "";
				if(i+1<argc)
					value = string(args[i+1]);
				options.push_back(CmdOption(tag, value));
			}
		}
	}

	bool specifiedOption(string tag){
		for(int i=0;i<options.size();i++){
			if(options.at(i).optionTag.compare(tag) == 0)
				return true;
		}
		return false;
	}

	string getValue(string tag){
		for(int i=0;i<options.size();i++){
			if(options.at(i).optionTag.compare(tag) == 0)
				return options.at(i).optionValue;
		}
		cerr << "Can't find option: " << tag << endl;
		exit(1);
	}

	virtual ~CmdArgs();
};

} /* namespace NSPtools */

#endif /* TOOLS_CMDARGS_H_ */
