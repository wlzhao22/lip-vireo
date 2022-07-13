/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**/

#if !defined SCRIPTPARSER_H
#define SCRIPTPARSER_H

#include <string>
#include <vector>
#include <list>
#include <map>


using namespace std;

typedef map<const char *,const char *>::value_type paraType;

class ScriptParser
{
public:
	ScriptParser()
	{
	}
	static void getDtTask(map<string,const char *> &paras,const char *fn);
	static void getKeyFrames(vector<const char*> &kflst,const char fname[]);
	virtual ~ScriptParser()
	{

	}
	static void test();

};

#endif
