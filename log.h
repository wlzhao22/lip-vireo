#ifndef LOG_H
#define LOG_H
#include "abstractdetector.h"

/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wan-Lei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
*
Speed-up version: Hessian+nrrSIFT VLAD* mAP=0.741
**/

class LoG :public AbstractDetector
{
private:

    //parameters for user

    static const int   MaxOctaves;
    static const int   SCALES;
    static const float _SIGMA;
    static const bool  INTERP_KEYS;
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
    static const int   maxIter;
    static const int   W0;
    static const float err_c ;
    int                numb;
    float              curveRatio, thresh, sigma;

    static const int NumOrient;
    static const int DEGPERBIN;
    static const float NwKpThresh;

    vector<vector<Image *> > BuildLOGOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves);

    vector<Image*>  BuildLoGScales(vector<Image *>  & gxScales,vector<Image *> &gyScales);
    bool isPeak(Image * aimage, Image * bimage,Image * cimage,  short &peak, int x, int y);

    vector<KeyPoint *> FindPeaksScales(int octave, vector<Image *> & LoGImages);
    vector<KeyPoint *> FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves);

    bool  adaptAffine(vector<vector<Image *> > &LoGOctaves, vector<vector<Image *> > & GxxOctaves,
                     vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves);

    float InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds);
    bool  InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal);

public:
    LoG();
    int keypDetect(const float sigma, const int level, const int scale);

public:
    bool paramsCheck();
    bool KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn);
    void writeKeypoint(const char*fn);

    static void test();
};
#endif
