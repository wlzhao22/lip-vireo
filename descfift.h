#ifndef DESCFIFT_H
#define DESCFIFT_H

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

class DescFIFT :public AbstractDescriptor{

    private:
        ///for descriptor
        static const int PatchSize;
        static const int NumOrient;
        static const int GRID;
        static const int DSize;
        static const int PSIFTLen;
        static const int PatchMag;
        vector<vector<float> > gMat;
        vector<vector<float> > cgmat;

    public:
        DescFIFT(DESC desc);
        int        buildDescriptor(const int kpnum,vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
        int        buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate);
        static int getWinSize(KeyPoint *key,float &sizeratio);
        int        calcDescPatch(KeyPoint *keyp,float *myWin,const int Size);
        int        calcPatchCurl(KeyPoint *keyp,float *myWin,const int Size);
        int        getFIFTDescriptor(const float *myWin, const float sc);
        int        getFIFTDescriptor(Image *win, const float sizeratio);
        static vector<vector<float> > GaussianWeight2D(const int winSize);
        static vector<vector<float> > GaussianWeight2D(const float sigma, const int radius);

        ~DescFIFT();

};
#endif
