#ifndef DESCPCASIFT_H
#define DESCPCASIFT_H

#include "abstractdescriptor.h"
#include "keypoint.h"

class DescPCASIFT: public AbstractDescriptor
{
private:
    ///for descriptor
    static const unsigned int PatchMag;
    static const unsigned int PCALen;
    static const int PatchSize = 41;
    static const int GPLEN = ((PatchSize-2)*(PatchSize-2)*2);
    static const int EPCALEN = 36;

    float  avgs[GPLEN]; ///< average local descriptor value.  Used in PCA.
    float  eigs[EPCALEN][GPLEN]; ///< PCA basis vectors
    float  *gradients;
    char   pcafn[512];

public:
    bool   paramsCheck();
    void   loadPCAMatrix(const char *pcafn);
    int    getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size);
    int    buildDescriptor(const int kpnum,vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
    int    buildPatchView(const int kpnum,vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
    static int getWinSize(KeyPoint *key, float &sizeratio);
    int    (DescPCASIFT::*computeDesc)(const float *win);
    int    getPCASIFTDesc(const float *myWin);
    int    makePatches(const float *myWin);

    static void test();
public:

    DescPCASIFT(DESC desc, const char *pcamapfn);
    virtual ~DescPCASIFT()
    {
        delete [] pcaFeat;
        delete [] gradients;
    }
};

#endif
