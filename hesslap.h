#ifndef HESSLAP_H
#define HESSLAP_H

#include "abstractdetector.h"

/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose


Speed-up version: HessLap+nrrSIFT VLAD*     Holidays        mAP=0.745
Speed-up version: HessLap+nrrSIFT VLAD*     Ox5k            mAP=0.373
Speed-up version: HessLap+nrSIFT   BoW*     Ox5k            mAP=0.354
Speed-up version: HessLap+nrSIFT   BoW*     Holidays        mAP=0.626
Speed-up version: HessLap+nrSIFT   HE+BoW*  Holidays        mAP=0.745
Speed-up version: HessLap+nrSIFT   HE+BoW*  Ox5k            mAP=0.399
**/

class HessLap: public AbstractDetector
{

private:

    int numb;
    float thresh;
    float sigma;

    static const bool  INTERP_KEYS;
    static const int   MaxOctaves;
    static const float INITSIGMA;
    static const int   SCALES;
    static const float _SIGMA;

    ///for detector
private:
    static const float sigma0;
    static const float k0;
    static const int   BORDER;
    static const int   THRESH;
    static const float mag;
    static const int   DEGREE;
    static const float cvtRatio;
    static const float thresh_ratio;
    static const int   maxIter;
    static const int   W0;
    static const float err_c ;

    static const float NwKpThresh;
    static const int   NumOrient;
    static const int   DEGPERBIN;

    vector<Image*>  BuildLoGScales(vector<Image *>  & gxScales,vector<Image *> &gyScales);
    vector<Image* > BuildHessScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales,vector<Image *> &gxyScales,vector<Image*> &LoGScales);

    vector<KeyPoint *> FindPeaksScales(const int octave, vector<Image*> HessImages, vector<Image *> & LoGImages);
    vector<KeyPoint*>  FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves);

    float InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds);
    bool  InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal);

    vector<vector<Image *> > BuildHessOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> > & GxyOctaves, vector<vector<Image *> > &LoGOctaves);
    vector<vector<Image *> > BuildLOGOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves);

    bool adaptScale(int &x, int &y, vector<vector<Image *> > & GxxOctaves, vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves, int &scale, int &octave);
    bool adaptAffine(vector<vector<Image *> > &HessOctaves, vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves);

    bool isSpatialPeak(Image *image,const int x,const int y);
    bool isScalePeak(Image * aimage, Image * bimage, Image * cimage, float &logVal, const int x, const int y);

public:
    HessLap();
    int keypDetect(const float sigma, const int level, const int scale);

public:
    bool paramsCheck();
    bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn);
    void writeKeypoint(const char*fn);

    static void test();
};

#endif
