#ifndef DESCFIND_H
#define DESCFIND_H

#include "abstractdescriptor.h"

class DescFIND: public AbstractDescriptor
{

private:
    vector<vector<float> > gmat;

    ///for descriptor
    static const int ORIENTATION;
    static const int GRID;
    static const int DSize;
    static const int DEGREE;
    static const int PatchSize;

    static const int PSIFTLen;
    static const int SIFT_MAG;
    static const int PatchMag;

    static const float inv_sqrt2;


public:
    DescFIND(DESC desc);

    int getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size);
    int buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate);
    int buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate);
    int getFINDDescriptor(const float *myWin, const int flip);
    static void test();
    static vector<vector<float> > GaussianWeight2D(const int winSize);

    ~DescFIND();

};
#endif
