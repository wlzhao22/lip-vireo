#ifndef HESSIAN_H
#define HESSIAN_H

#include "abstractdetector.h"

/**
*@author Wanlei Zhao
*@version 1.30
*All rights reserved by Wan-Lei Zhao
*
*Anyone receives this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
=================Performance evaluatoin=====================
Gaussian-Norm:    Hessian+nrrSIFT VLAD* Holidays    mAP=0.624
Speed-up version: Hessian+nrrSIFT VLAD* Holidays    mAP=0.698
Gaussian-Norm:    Hessian+nrrSIFT VLAD* ox5k        mAP=0.332
Speed-up version: Hessian+nrrSIFT VLAD* ox5k        mAP=0.355
**/

class Hessian: public AbstractDetector
{

private:

    int numb;
    float thresh;
    float sigma;
    static const int MaxOctaves;
    static const int SCALES;
    static const float _SIGMA;
    static const bool INTERP_KEYS;
    static const float INITSIGMA;

    //for detector
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

    static const int NOrient;
    static const int DEGPERBIN;
    static const float NwKpThresh;

    vector<Image*>  BuildLoGScales(vector<Image *> & gxScales, vector<Image *> &gyScales);
    vector<Image* > BuildHessScales(vector<Image *> &gxxScales, vector<Image *> &gyyScales,
                                    vector<Image *> &gxyScales);
    vector<KeyPoint *> FindPeaksScales(const int octave, vector<Image*> &HessImages, vector<Image*> &GxxImages, vector<Image*> &GyyImages);
    vector<KeyPoint*>  FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves);

    float InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds);
    bool InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal);

    vector<vector<Image *> > BuildHessOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,
                                              vector<vector<Image *> > & GxyOctaves);
    vector<vector<Image *> > BuildLOGOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves);

    bool isSpatialPeak(Image *image,const int x,const int y);
    bool isScalePeak(Image * aimage, Image * bimage, Image * cimage, const int x, const int y);

public:
    Hessian();
    ~Hessian();
    int keypDetect(const float sigma,const int level,const int scale);

public:
    bool paramsCheck();
    bool KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn);
    void writeKeypoint(const char*fn);

    static void test();
};

#endif
