#ifndef DESCRIFT_H
#define DESCRIFT_H

#include "abstractdescriptor.h"
#include "keypoint.h"
#include "filter.h"
#include "nnitem.h"
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

class DescRIFT :public AbstractDescriptor
{

private:
    static const int Ori;
    static const int BLCK;

    static const int Dist;
    static const int DistF;
    static const int PatchMag;
    static const int PatchSize;
    static const float Theta_per_bin;
    static const float Blck_Size;
    static const float Dist_per_bin;
    static const int factor;
    static const float alpha;
    //float *Sigmas2;
    //float *Sigmas3;

public:
    DescRIFT(DESC desc);
    int getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size);
    int (DescRIFT::*computeDesc)(const float *win);
    int buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate);
    int buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate);

    static int getWinSize(KeyPoint *key,float &sizeratio);
    int getERIFTDescriptor(const float *myWin);
    int getRIFTDescriptor(const float *myWin);

    ~DescRIFT();

};

#endif
