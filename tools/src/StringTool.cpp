/*
 * StringTool.cpp
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#include "tools/StringTool.h"

namespace NSPtools {

string trimString(const string& s)
{
     if (s.empty()) {
             return s;
     }
     int begin=0, end=s.length()-1;
     for (int i = 0; i < s.length(); i++)
     {
             if (s.at(i) == ' ')
                     begin++;
             else
                     break;
     }
     for (int i = s.length() - 1; i >= 0; i--)
     {
             if (s.at(i) == ' ')
                     end--;
             else
                     break;
     }
     return s.substr(begin, end - begin + 1);
 }

void splitString(const string& s, const string& delim, vector<string >* ret)
{
	size_t last = 0;
	size_t len = s.length();
	size_t index = s.find_first_of(delim, last);
	ret->clear();
	while (index != std::string::npos)
	{
	     if(index > last)
	    	 ret->push_back(s.substr(last, index - last));
	     last = index + 1;
	     index = s.find_first_of(delim, last);
	}
	if (last < len)
	{
	       ret->push_back(s.substr(last, len - last));
	}


}


} /* namespace NSPdesignseq */
