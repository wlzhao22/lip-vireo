#include "iotool.h"
#include "vstring.h"
#include "config.h"

#include <sys/stat.h>
#include <dirent.h>
#include <dirent.h>
#include <cstring>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <cmath>

#ifndef LINUX
#define LINUX
#endif

using namespace std;

bool IOTool::IS_Dir(const char *dir)
{
    DIR *_mydir = opendir(dir);
    if(_mydir != NULL)
    {
        closedir(_mydir);
        return true;
    }else
    {
        return false;
    }
}
bool IOTool::IS_FILE(const char *fn)
{
    ifstream in;
    in.open(fn);
    if(in.fail())
    {
        return false;
    }else
    {
        in.close();
    }
    return true;
}


bool IOTool::creatDIR(const char *dirname)
{
    DIR *dirdp;
    if(dirname != NULL)
    {
        if((dirdp  = opendir(dirname)) == NULL)
        {
#ifdef LINUX
            if(mkdir(dirname, ACCESS_MODE) != 0)
            {
                cout<<"Creating directory '"<<dirname<<"' failed!\n";
                exit(1);
            }
#else
            if(mkdir(dirname) != 0)
            {
                cout<<"Creating directory '"<<dirname<<"' failed!\n";
                exit(1);
            }
#endif
        }
        else
        {
            closedir(dirdp);
        }
    }
	return true;
}

bool IOTool::getKeypoint(const char *keyfn, vector<KeyPoint *> &kplst)
 {
	int col, num;
	KeyPoint *crnt_pt;
	ifstream *keyfile=new ifstream();
	keyfile->open(keyfn,ios::in);

	if(keyfile->fail())
	{
		cout<<"Fail to read "<<keyfn<<endl;
		return false;
	}

	(*keyfile)>>num;
	(*keyfile)>>col;
    unsigned int counter = 0;
	unsigned int pt_num = num;

	assert(col >= 6);
	float x = 0, y = 0;

	while(!keyfile->eof())
	{
		if(counter == pt_num)
			break;

		crnt_pt = new KeyPoint();
		crnt_pt->KP = true;
		if(col == 6)
		{
		    (*keyfile)>>x;
		    crnt_pt->x = (int)round(x);
		    (*keyfile)>>y;
		    crnt_pt->y = (int)round(y);
		    (*keyfile)>>crnt_pt->a;
		    (*keyfile)>>crnt_pt->b;
		    (*keyfile)>>crnt_pt->c;
		    (*keyfile)>>crnt_pt->iscale;
		    crnt_pt->ori = 0;
		    crnt_pt->funcVal = 1;
		}else if(col == 7)
		{
            (*keyfile)>>x;
		    crnt_pt->x = (int)round(x);
		    (*keyfile)>>y;
		    crnt_pt->y = (int)round(y);
		    (*keyfile)>>crnt_pt->a;
		    (*keyfile)>>crnt_pt->b;
		    (*keyfile)>>crnt_pt->c;
		    (*keyfile)>>crnt_pt->iscale;
            (*keyfile)>>crnt_pt->ori;
            crnt_pt->funcVal = 1;
		}else if(col == 8){
            (*keyfile)>>x;
		    crnt_pt->x = (int)round(x);
		    (*keyfile)>>y;
		    crnt_pt->y = (int)round(y);
		    (*keyfile)>>crnt_pt->a;
		    (*keyfile)>>crnt_pt->b;
		    (*keyfile)>>crnt_pt->c;
		    (*keyfile)>>crnt_pt->iscale;
            (*keyfile)>>crnt_pt->ori;
            (*keyfile)>>crnt_pt->funcVal;
		}
        kplst.push_back(crnt_pt);
		counter++;
	}//end while
	keyfile->close();
	delete keyfile;
	return true;
}


float *IOTool::load_pcaMat(const char *pcaMatFn, vector<float> &means,
                           vector<float> &variances, unsigned int &row, unsigned int &col)
{
    float *matrix = NULL;
    ifstream *inStrm = new ifstream(pcaMatFn);

    if(!inStrm->is_open())
    {
        row = 0;
        col = 0;
        cout<<"File "<<pcaMatFn<<" not found!\n";
        return NULL;
    }

    (*inStrm)>>row;
    (*inStrm)>>col;

    if(row > 0)
    {
        matrix = new float[row*col];
    }
    else
    {
        inStrm->close();
        return NULL;
    }

    unsigned int counter = 0, i = 0, j = 0;
    float val;
    for(i = 0; !inStrm->eof() && i < col; i++)
    {
        (*inStrm)>>val;
        means.push_back(val);
    }

    for(i = 0; !inStrm->eof() && i < row; i++)
    {
        (*inStrm)>>val;
        val = sqrt(val);
        variances.push_back(val);
    }
    i = 0;
    while(!inStrm->eof() && i < row)
    {
        for(j = 0; j < col; j++)
        {
            (*inStrm)>>matrix[counter];
            counter++;
        }
        i++;
    }

    inStrm->close();
    return matrix;
}

float *IOTool::loadMat(const char *matFn, unsigned int &row, unsigned int &col)
{
    float *matrix = NULL;
    ifstream *inStrm = new ifstream(matFn);

    if(!inStrm->is_open())
    {
        row = 0;
        col = 0;
        cout<<"File "<<matFn<<" not found!\n";
        return NULL;
    }

    (*inStrm)>>row;
    (*inStrm)>>col;

    if(row > 0)
    {
        matrix = new float[row*col];
    }
    else
    {
        inStrm->close();
        return NULL;
    }

    unsigned int idx = 0, i = 0, j = 0;
    while(!inStrm->eof() && i < row)
    {
        for(j = 0; j < col; j++)
        {
            (*inStrm)>>matrix[idx];
            idx++;
        }
        i++;
    }

    inStrm->close();
    return matrix;
}
