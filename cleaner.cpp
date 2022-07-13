#include "cleaner.h"
#include <iostream>

using namespace std;


bool Cleaner::releaseScales(vector<Image*> &scales)
{
    vector<Image *>::iterator it;
    Image *tempim = NULL;
    for(it = scales.begin(); it != scales.end(); it++)
    {
        tempim = *it;
        delete tempim;
    }
    scales.clear();
    return true;
}

bool Cleaner::releaseBoards(vector<Board*> &scales)
{
    vector<Board *>::iterator it;
    Board *tempim = NULL;
    for(it = scales.begin(); it != scales.end(); it++)
    {
        tempim = *it;
        delete tempim;
    }
    scales.clear();
    return true;
}

bool Cleaner::releaseOctaves(vector<vector<Image*> > &octaves)
{
    vector<Image *>::iterator it;
    Image *tempim;
    for (unsigned i = 0; i < octaves.size(); i++)
    {
        for(it=octaves[i].begin(); it!=octaves[i].end(); it++)
        {
            tempim = *it;
            delete tempim;
        }
        octaves[i].clear();
    }
    octaves.clear();
    return true;
}

bool Cleaner::releaseOctaves(vector<vector<ImageSymmMat*> > &octaves)
{
    vector<ImageSymmMat *>::iterator it;
	ImageSymmMat *tempim;
    for (unsigned i = 0; i < octaves.size(); i++)
	{
	    for(it = octaves[i].begin(); it!=octaves[i].end(); it++)
		{
			tempim = *it;
			delete tempim;
		}
		octaves[i].clear();
	}
	octaves.clear();
	return true;
}

bool Cleaner::clearParas(map<string,const char*> &paras)
{
    map<string, const char*>::iterator mit;
    string crnt_str;
    for(mit =  paras.begin(); mit != paras.end(); mit++)
    {
        crnt_str = mit->first;
        crnt_str.clear();
        const char *val = mit->second;
        delete [] val;
    }
    paras.clear();
    return true;
}

void Cleaner::clearCurves(vector<HCPoint*> &curves)
{
    vector<HCPoint*>::iterator it;
    CPoint  *crnt_pt = NULL, *nxt_pt;
    HCPoint *crnt_head = NULL;

    /**1. calc curvature**/
    for(it = curves.begin(); it != curves.end(); it++)
    {
        crnt_head = *it;
        crnt_pt   = crnt_head->next;
        while(crnt_pt != NULL)
        {
            nxt_pt  = crnt_pt->next;
            delete crnt_pt;
            crnt_pt = nxt_pt;
        }
        delete crnt_head;
    }

    curves.clear();
    return ;
}

void Cleaner::clearSegment(vector<LPoint*> *segment)
{
    vector<LPoint*>::iterator vit;
    LPoint* crnt_pt = NULL;
    for(vit = segment->begin(); vit != segment->end(); vit++)
    {
        crnt_pt = *vit;
        delete crnt_pt;
    }
    segment->clear();
    return ;
}

void Cleaner::clearAdjTab(map<unsigned int, list<unsigned int>*> &AdjTable)
{
    map<unsigned int, list<unsigned int>*>::iterator mit;
    list<unsigned int>* crnt_lst;
    list<unsigned int>::iterator lit;

    for(mit = AdjTable.begin(); mit != AdjTable.end(); mit++)
    {
        crnt_lst = mit->second;
        crnt_lst->clear();
        delete crnt_lst;
    }
    AdjTable.clear();
    return ;
}

void Cleaner::clearSegTab(map<unsigned int, map<unsigned int, vector<LPoint*>* >* >   &SegTable)
{
    LPoint* crnt_pt;
    map<unsigned int, vector<LPoint*>* >* crnt_map;
    vector<LPoint*>* crnt_seg;
    vector<LPoint*>::iterator it;
    map<unsigned int, map<unsigned int, vector<LPoint*>* >* >::iterator mit;
    map<unsigned int, vector<LPoint*>* >::iterator cmit;
    for(mit = SegTable.begin(); mit != SegTable.end(); mit++)
    {
        crnt_map = mit->second;
        for(cmit = crnt_map->begin(); cmit != crnt_map->end(); cmit++)
        {
            crnt_seg = cmit->second;
            for(it = crnt_seg->begin(); it != crnt_seg->end(); it++)
            {
                crnt_pt = *it;
                delete crnt_pt;
            }
            crnt_seg->clear();
        }
        crnt_map->clear();
    }
    SegTable.clear();
}

void Cleaner::clearDistMap(map<unsigned int, map<unsigned int, float >* > &distMap)
{
    map<unsigned int, map<unsigned int, float >* >::iterator mit;
    map<unsigned int, float > *DMap;
    for(mit = distMap.begin(); mit != distMap.end(); mit++)
    {
        DMap = mit->second;
        DMap->clear();
        delete DMap;
    }
    distMap.clear();
}

void Cleaner::clearKeyPoint(vector<KeyPoint*> &mykps)
{
    vector<KeyPoint*>::iterator itkp;
    KeyPoint* crntkp;

    for(itkp = mykps.begin(); itkp != mykps.end(); itkp++)
    {
        crntkp = *itkp;
        delete crntkp;
    }
    mykps.clear();
    return ;
}

void Cleaner::clearRangeMap(map<int, vector<int> *> &rangeTab)
{
    map<int, vector<int> *>::iterator mit;
    vector<int> *crnt_pts;

    for(mit = rangeTab.begin(); mit != rangeTab.end(); mit++)
    {
        crnt_pts = mit->second;
        crnt_pts->clear();
        delete crnt_pts;
    }
    return ;
}

void Cleaner::clearXPoint(vector<XPoint*> &segment)
{
    vector<XPoint*>::iterator vit;
    XPoint *crnt_pt;
    for(vit = segment.begin(); vit != segment.end(); vit++)
    {
        crnt_pt = *vit;
        delete crnt_pt;
    }
    segment.clear();
}

void Cleaner::clear2DArray(vector<vector<float> > &array2D)
{
    vector<vector<float> >::iterator vit;
    for(vit = array2D.begin(); vit != array2D.end(); vit++)
    {
        vector<float> array = *vit;
        array.clear();
    }
    array2D.clear();
    return ;
}

void Cleaner::clear1DArray(vector<float> &array1D)
{
    array1D.clear();
    return ;
}
