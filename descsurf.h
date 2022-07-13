#ifndef DESCSURF_H
#define DESCSURF_H

#include "abstractdescriptor.h"
#include "intimage.h"

/**
AoD:  Aggregation on Derivatives
SURF: Speed-Up Robust Features

Please always normalize the surf descriptor, it looks significantly better

@author: Wan-lei Zhao
@date:   Jan-15-2012

I found it performs less than getDx and getDy when I use getHaarX and getHaarY,
which is a little bit surperising! The latter is supposed to be the precise
way to calculate Haar wavelet

***/
/**
Codes from OpenSURF doesn't work well (I tried a lot of times)
http://code.google.com/p/opensurf1/
**/

class DescSURF: public AbstractDescriptor
{
private:
    static const int PatchSize;
    vector<vector<float> > gmat;

    ///for descriptor
    static const int GRID;
    static const int CHANNEL;
    static const int iScale0;
    static const int Radius;
    static const float dScale0;

    int rChnl;
    IntImage *intImg;
public:
    DescSURF(DESC desc);

    int      (DescSURF::*computeDesc)(KeyPoint *crnt_kpt);
    int      getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size);
    int      getAoDDesc(KeyPoint   *crnt_kpt);
    int      getSURFDesc(KeyPoint  *crnt_kpt);
    int      getESURFDesc(KeyPoint *crnt_kpt);
    int      getSURFDescriptor(KeyPoint *crnt_kpt, bool upright);
    float    gaussian(float x, float y, float sig);
    int      buildDescriptor(const int kpnum, vector<KeyPoint*> &kps,
                             const char *descfn, const float resize_rate);
    int      buildPatchView(const int kpnum, vector<KeyPoint*> &kps,
                             const char *dvfn, const float resize_rate);
protected:
    float    getDx(const int x0, const int y0, const int rd0, const int ts);
    float    getDy(const int x0, const int y0, const int rd0, const int ts);
    float    getHaarX(const int x0, const int y0, const int rd0, float rMat[2]);
    float    getHaarY(const int x0, const int y0, const int rd0, float rMat[2]);
    void     buildIntImage();
    virtual ~DescSURF();
};

#endif
