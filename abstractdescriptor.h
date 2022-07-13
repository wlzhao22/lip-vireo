#ifndef ABSTRACTDESCRIPTOR_H
#define ABSTRACTDESCRIPTOR_H

#include "idetector.h"
#include "keypoint.h"
#include "image.h"
#include "cimage.h"

#include <cassert>
#include <vector>

using namespace std;

/*****************************
Interface to descriptors

@author: Wanlei Zhao
@date:   Dec-2011


******************************/

class AbstractDescriptor
{

protected:
    float *descWin;
    float *featsBin;
    float *pcaFeat;

    AbstractImage *crntImg;
    vector<vector<Image*> > ImageOctaves;

    bool properties[IDetector::Numb_PROP];
    int  KP_PROP_SZ[IDetector::Numb_PROP];
    bool AFF_OUT;
    DESC_FMT _out_format;

protected:

    DESC descoption;
    int featLen, gDim;
    Detector DETECTOR;

    static const float SIGMA;
    static const int   SCALESPEROCTAVE;
    static const float SmthRatio;

public:
    AbstractDescriptor()
    {
        this->descWin      = NULL;
        this->featsBin     = NULL;
        this->pcaFeat      = NULL;
        this->crntImg      = NULL;
        this->AFF_OUT      = false;
        this->_out_format  = _vireo_fmrt;
        this->DETECTOR     = corner;
        KP_PROP_SZ[0] = 1;
        KP_PROP_SZ[1] = 1;
        KP_PROP_SZ[2] = 1;
        KP_PROP_SZ[3] = 2;
        KP_PROP_SZ[4] = 32;
    }

public:
    bool setupImage(AbstractImage *srcImg);
    bool setupImage(CImage *srcImg);
    bool setOutFormat(bool kp_properties[], bool aff_out, DESC_FMT _out_fmt)
    {
        for(int i = 0; i < IDetector::Numb_PROP; i++)
        {
            properties[i] = kp_properties[i];
        }
        this->AFF_OUT      = aff_out;
        this->_out_format  = _out_fmt;
        return true;
    }

    bool setupOctaves(vector<vector<Image*> > &GOctaves);


    int  getDescPatch1(KeyPoint *keyp, float *myWin, const int Size);
    int  getDescPatch(KeyPoint *keyp,  float *myWin, const int Size);

    int  getDescAffPatch(KeyPoint *keyp, float *myWin, const int Size);
    int  getDescPatch3C(KeyPoint *keyp, float *myWin, const int Size);

    void saveDescVGG(const KeyPoint* crnt_kpt, const float *feature,
                     const int dim,  const float resize, ofstream &outStrm);
    void saveDescVireo(const KeyPoint* crnt_kpt, const float *feature, const int dim,
                       const int idx0,  const float resize, ofstream &outStrm);
    void saveDescTrain(const KeyPoint* crnt_kpt, const float *feature, const int dim,
                       const float resize, ofstream &outStrm);

    DESC getDescOption()
    {
        return descoption;
    }

    virtual int buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate) = 0;
    virtual int buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate) = 0;

    static  int findFlip(const float thetas[], const int dbin, const int nbin, float &ratio);
    int getCurl(KeyPoint *keyp,float *myWin,const int Size);

    virtual ~AbstractDescriptor()
    {
        if(descWin != NULL)
        {
            delete [] descWin;
        }
    }
};
#endif
