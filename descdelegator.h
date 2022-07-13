#ifndef DESCDELEGATOR_H
#define DESCDELEGATOR_H

#include "idetector.h"

#include <vector>

/************************************************************
@author: Wanlei Zhao
@date:   27-Mar-2012

Delegator to many local descriptors:
RIFT
SPIN
SIFT

We assume the input patch is normalized intensity levels
***********************************************************/

using namespace std;

class DescDelegator
{
private:
    ///for SIFT
    static const int GRID;
    static const int DSize;
    static const int DEGREE;
    static const int PSIFTLen;
    static const int PatchMag;
    static const int NumOrient;

    ///for RIFT
    static const int   Ori;
    static const int   DistF;
    static const float Theta_per_bin;

    ///for SPIN
    static const int Ints;
    static const int Dist;
    static const float alpha;
    static const float belta;
    static const int   CLR_DEPTH;

    DESC descoption;
    int featLen;
    int PatchSize;

protected:
    vector<vector<float> > GaussianWeight2D(const int winSize);
    vector<vector<float> > gmat;
    int (DescDelegator::*extractFeat)(float *feats, const float *myWin);
    int getSIFTDescriptor(float *feats, const float *myWin);
    int getSPINDescriptor(float *feats, const float *myWin);
    int getRIFTDescriptor(float *feats, const float *myWin);
    int getGBLURDescriptor(float *feats, const float *myWin);

public:
    DescDelegator(DESC desc0, const int PatchSize);
    int getFeatDim();
    int getDescriptor(float *feats, const float *myWin);

    virtual ~DescDelegator();
};

#endif
