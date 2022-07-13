/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
Speed-up version: DoG+nrrSIFT VLAD* mAP=0.712
**/

#ifndef DOG_H
#define DOG_H

#include <vector>
#include <cmath>

#include "abstractdetector.h"
#include "keypoint.h"
#include "filter.h"
#include "image.h"

class DoG: public AbstractDetector
{

protected:

    static const int MaxOctaves;

    static const int   PatchMag;
    static const float INITSIGMA;
    static const float _SIGMA;
    static const int   SCALESPEROCTAVE;

    static const int   NumOrient;
    static const int   DEGPERBIN;
    static const float NwKpThresh;
    static const float LOW_CONTRAST_THRESH;
    static const float EDGE_PEAK_THRESH;
    static const float EDGE_PEAK_L;
    static const int   INTERP_KEYS;

    static const int   BORDER;

    static const int   featLen;
    static const int   ORIENTATION;
    static const int   GRID;
    static const int   DSize;
    static const int   DEGREE;

    static const int   kpnum0;
    static const int   PatchSize;
    static const float mag;

    float *featsBin;
    float *descWin ;

    float sigma, thresh;
    int   numb, RESIZE;;
       /**
        *@param sigma initial convolution window
        *@param threshold for hessian, it is the minimum value required for Det of hessian matrix
        *@param numb you can set a upbound number of points you need, or set a minus one, means no upperbound
        **/

private:

    bool               isPeak(Image * aimage, Image * bimage,Image * cimage, unsigned char &valley,int x, int y);
    vector<KeyPoint *> FindPeaksScales(int octave, vector<Image *> & DImages,vector<Image *> & GImages);
    bool               isEdgePeak(int x, int y, Image * image);
    float              InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds) ;
    bool               InterpKey(int x, int y, int s, vector<Image *> & DImages, float * fx, float * fy, float * fs,float *dogVal);

protected:

    Image                    *ScaleInitImage(Image * image);
    vector<vector<Image *> > BuildDOGOctaves(vector<vector<Image *> > & GOctaves);
    vector<Image *>          BuildDOGScales(vector<Image *> & LImages);
    void                     BuildOctaves(Image * image,vector<vector<Image *> > &octaves);
    vector<Image *>          BuildGaussianScales(Image * image);
    void                     FindPeaksOctaves(vector<vector<Image *> > & DOctaves,vector<vector<Image*> > &GImage);
    vector<KeyPoint *>       FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves);

public:
    DoG();
    bool paramsCheck();
    ~DoG();
    bool KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn);
    static void test();

};
#endif
