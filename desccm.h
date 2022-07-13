#ifndef DESCCM_H
#define DESCCM_H

#include "abstractdescriptor.h"
#include "cimage.h"

/******************************
* Color Moment is used as local descriptor for keypoint
* After several runs of experiment, it simply does not work
* I just keep it here. Should not be used, it is highly indistinctive
*
******************************/

class DescCM : public AbstractDescriptor
{

private:
    static const int CM_NUM;
    static const int PatchSize;
    float *CMfeats;
    float *block;
    float *descPatch; //3 channels
    int bl_size;
    int desc_patch_size;

public:
    DescCM(DESC desc);
    virtual ~DescCM();

public:

    //bool setupImage(CImage *srcImg);
    bool reShape();
    bool reShape(const int channel,const int width,const int height);

    int buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate);
    int buildPatchView(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate)
    {
        return 1;
    }

    void getLocalColorMoment(const float *patch);
    void getLocalColorMoment(CImage *myimg,const int height,const int width);
    void RGB2LUV(const float RGB[],float LUV[]);

    static void test();
};

#endif // DESCCM_H
