#ifndef DESCSIFT_H
#define DESCSIFT_H

#include <vector>

#include "abstractdescriptor.h"
#include "keypoint.h"
#include "image.h"

using namespace std;

/**********************************************************************
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao and the patent holder
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
************************************************************************/

class DescSIFT :public AbstractDescriptor
{

private:
    vector<vector<float> > gMat;

    ///for descriptor
    static const int NumOrient;
    static const int GRID;
    static const int DSize;
    static const int PatchSize;
    static const int PSIFTLen;

    static const int PatchMag;

public:
    DescSIFT(DESC desc);

    int buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
    int buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);

    int getSIFTDescriptor(const float *myWin, const float sigma);
    int getSIFTDescriptor(Image *win,const float sizeratio);
    int getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size);

    static void test();
    static vector<vector<float> > GaussianWeight2D(const float sigma, const int radius);
    static vector<vector<float> > GaussianWeight2D(const int winSize);
    static void getPatch(const char*srcimg, const char *srcfn, const char *dstpatch);
    static bool getPixGradient(const int x,const int y,Image *win,const float sizeratio,float &mag,float &theta);
    static bool getPixGradient(const int x,const int y,const float *pixels,const int size,float &mag,float &theta);

    ~DescSIFT();

};
#endif
