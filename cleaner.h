#ifndef CLEANER_H
#define CLEANER_H

#include "imagesymmmat.h"
#include "keypoint.h"
#include "board.h"
#include "image.h"

#include <string>
#include <vector>
#include <list>
#include <map>

using namespace std;

class Cleaner {

public:
    static bool releaseBoards(vector<Board*> &scales);
    static bool releaseScales(vector<Image*> &scales);
    static bool releaseOctaves(vector<vector<Image*> > &octaves);
    static bool releaseOctaves(vector<vector<ImageSymmMat*> > &octaves);
    static bool clearParas(map<string, const char*> &paras);
    static int  cleanDupKP(vector<KeyPoint *> &kps, const int width, const int height);
    static void clearCurves(vector<HCPoint*> &curves);
    static void clearSegment(vector<LPoint*> *segment);
    static void clearAdjTab(map<unsigned int, list<unsigned int>*> &AdjTable);
    static void clearSegTab(map<unsigned int, map<unsigned int, vector<LPoint*>* >* >   &SegTable);
    static void clearDistMap(map<unsigned int, map<unsigned int, float >* > &distMap);
    static void clearKeyPoint(vector<KeyPoint*> &kps);
    static void clearXPoint(vector<XPoint*> &segment);
    static void clearRangeMap(map<int, vector<int> *> &rangeTab);
    static void clear2DArray(vector<vector<float> > &array2D);
    static void clear1DArray(vector<float> &array1D);
};

#endif
