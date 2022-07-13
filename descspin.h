#ifndef DESCSPIN_H
#define DESCSPIN_H

#include "abstractdescriptor.h"
#include "keypoint.h"
#include "filter.h"
#include "image.h"

#include <vector>

using namespace std;

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

class DescSPIN :public AbstractDescriptor
{

private:
    //for descriptor
    static const int Ints;
    static const int Dist;
    static const float alpha;
    static const float belta;
    static const int PatchMag;
    static const int PatchSize;
    static const int CLR_DEPTH;

public:
    DescSPIN(DESC desc);
    int getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size);
    int getOctavedDescPatch(KeyPoint * key,const Image * blur,float *descWin);

    int buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
    int buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate);

    static int getWinSize(KeyPoint *key,float &sizeratio);

    int getSPINDescriptor(const float *myWin);
    int getSPINDescriptor(Image *win,const float sizeratio);

    ~DescSPIN();

};
#endif
