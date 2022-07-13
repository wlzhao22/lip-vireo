#include "scriptparser.h"
#include <iostream>
#include <cstring>
#include <cstdio>

#include "vstring.h"

using namespace std;

void ScriptParser::getDtTask(map<string,const char *> &paras,const char *fn)
{
    FILE *fp=fopen(fn,"r");
    if(!fp)
    {
		printf("Script %s can not be open.\n",fn);
		return ;
	}

    int nargs = 0;
	char result[100] = "";
    string fst_str = "", sec_str = "";


	while(!feof(fp))
	{
		nargs = fscanf(fp, "%s", result);

		if(nargs != 1)
		continue;

        if(!strcmp(result, "") || strlen(result) < 3)
           continue;

		VString::parsePair(result,fst_str,sec_str,'=');
		char *pairf = new char[fst_str.length()+1];
		strcpy(pairf, fst_str.c_str());
		char *pairs = new char[sec_str.length()+1];
		strcpy(pairs, sec_str.c_str());
		paras.insert(paraType(pairf, pairs));
		fst_str.erase(fst_str.begin(), fst_str.end());
		sec_str.erase(sec_str.begin(), sec_str.end());
                strcpy(result,"");
	}
	fclose(fp);
	return ;
}

void ScriptParser::getKeyFrames(vector<const char*> &kflst,const char fname[])
{
	printf("Read KeyFrames set...\n");

	FILE *fp=fopen(fname,"r");
	if(!fp)
	{
		printf("%s File can not be open.\n",fname);
		return ;
	}

	char result[60];
	char *crntKeyFrm;
	int i = 0, nargs = 0;
	while(!feof(fp))
	{
		nargs = fscanf(fp, "%s", result);

		if(nargs != 1)
		continue;

		crntKeyFrm = new char[60];
		strcpy(crntKeyFrm, result);
		//out<<result<<endl;
		kflst.push_back(crntKeyFrm);
		i++;
	}
	//cout<<"Count: "<<i<<endl;
	fclose(fp);
	return ;
}

void ScriptParser::test()
{

}
