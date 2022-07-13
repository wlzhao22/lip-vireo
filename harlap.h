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

#ifndef HARLAP_H
#define HARLAP_H

#include "abstractdetector.h"
#include "scriptparser.h"
#include "imagesymmmat.h"
#include "keypoint.h"
#include "image.h"

#include <vector>
#include <cmath>

using namespace std;

/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
Speed-up version: Harlap+nrrSIFT VLAD* mAP=72.9
**/


class HarLap:public AbstractDetector{

    private:
        int numb;
        float thresh;
        float sigma;

        /**
        *factor for differetial scale comparing with integration scale,fixed to 0.7
        *according to M.K.'s paper
        ***/

        static const float dfactor;
        static const int MaxOctaves;
        static const int SCALES;
        static const float _SIGMA;
        static const bool INTERP_KEYS;
        static const float INITSIGMA;
        static const int maxIter;
        static const int W0;
        static const float err_c ;
        static const float cvtRatio;

        /** check Mickolajczyk Krystian's IJCV'04 paper "Scale & Affine Invariant Interest Point Detectors" **/
        /**for ellipse adaption***/

     private:
        static const float sigma0;
        static const float k0;
        static const int BORDER;
        static const int THRESH;
        static const float mag;
        static const int DEGREE;

        static const int NumOrient;
        static const int DEGPERBIN;
        static const float NwKpThresh;

        vector<Image*>  BuildLoGScales(vector<Image *>  & gxxScales,vector<Image *> &gyyScales);
        vector<Image* > BuildHarrisScales(vector<Image *>  &gxScales,vector<Image *> &gyScales);
        vector<Image* > BuildHarrisScales(vector<Image *>  &gxScales, vector<ImageSymmMat*> &paraMats, vector<Image *> &gyScales);

        vector<KeyPoint *> FindPeaksScales(const int octave, vector<Image*> HessImages, vector<Image *> & LoGImages);
        vector<KeyPoint *> FindPeaksScales(const int octave, vector<Image*> HarImages, vector<ImageSymmMat*> &paraMats, vector<Image *> & LoGImages);
        vector<KeyPoint *> FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves);

    public:
        bool adaptAffine(vector<vector<Image *> > &HarrisOctaves, vector<vector<ImageSymmMat *> > & pmatOctaves);
        vector<vector<Image *> > BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves);
        vector<vector<Image *> > BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves, vector<vector<ImageSymmMat *> > & pmatOctaves, vector<vector<Image *> > & GyOctaves);
        vector<vector<Image *> > BuildLOGOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves);

        bool isSpatialPeak(Image *image,const int x,const int y);
        bool isScalePeak(Image * aimage, Image * bimage,Image * cimage, float &logVal,const int x, const int y);

        HarLap();
        int keypDetect(const float sigma,const int level,const int scale);

     public:
        bool paramsCheck();
        bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn);
        void writeKeypoint(const char *fn);

        static void test();
};
#endif
