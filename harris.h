/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**/

#ifndef HARRIS_H
#define HARRIS_H
#include "abstractdetector.h"

class Harris:public AbstractDetector
{
private:

    ///parameters for user
    int   numb;
    float thresh;
    float sigma;

    ///for detector
private:
    static const float sigma0;
    static const float _SIGMA;
    static const float k0;
    static const int   BORDER;
    static const int   SCALES;
    static const int   THRESH;
    static const float mag;
    static const float dfactor;
    static const int   DEGREE;
    static const int   NumOrient;
    static const int   DEGPERBIN;
    static const float NwKpThresh;


public:
    Harris();
    int keypDetect(const float dsigma, const int level, const int scale);

    ///inherited from AbstractDetector
public:
    bool paramsCheck();
    bool                     isSpatialPeak(Image *image, const int x, const int y);
    vector<KeyPoint*>        FindPeaksScales(const int octave, vector<Image*> HarImages);
    vector<KeyPoint*>        FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves);
    bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn);
    vector<Image* >          BuildHarrisScales(vector<Image *>  &gxScales, vector<Image *> &gyScales);
    vector<vector<Image *> > BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves,
                                                        vector<vector<Image *> > & GyOctaves);

    static void test();

};
#endif
