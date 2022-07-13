#ifndef DESCASIFT_H
#define DESCASIFT_H

#include "abstractdescriptor.h"
#include "keypoint.h"
#include "filter.h"
#include "image.h"

#include <vector>

using namespace std;

class DescASIFT: public AbstractDescriptor
{
   private:
        ///for descriptor
        static const int PatchSize;
        static const int NumOrient;
        static const int GRID;
        static const int DSize;
        static const int DEGREE;
        static const int PSIFTLen;
        static const int PatchMag;
        vector<vector<float> > gMat;

    public:
        DescASIFT(DESC desc);
        int        buildDescriptor(const int kpnum,vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
        int        buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
        int        calcDescPatch(KeyPoint *keyp,float *myWin,const int Size, const bool _NORM_);
        int        getSIFTDescriptor(const float *myWin);
        static vector<vector<float> > GaussianWeight2D(const int winSize);

        virtual ~DescASIFT();
};

#endif
