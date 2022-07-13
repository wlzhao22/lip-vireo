#ifndef DESCPVIEW_H
#define DESCPVIEW_H

#include "abstractdescriptor.h"
#include "idetector.h"

class DescPView: public AbstractDescriptor
{
    private:
        static const int PatchSize;
        static const int PatchMag;
        static const int FIFT_MAG;
        vector<vector<float> > cgmat;
    public:
        DescPView(DESC desc);
        vector<vector<float> > GaussianWeight2D(const int winSize);
        int calcDescPatch(KeyPoint *keyp, float *myWin, const int Size, const int id);
        int calcPatchCurl(KeyPoint *keyp, float *myWin, const int Size);
        int buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate);
        virtual ~DescPView();
};

#endif

