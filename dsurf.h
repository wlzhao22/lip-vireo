#ifndef DSURF_H
#define DSURF_H

#include "abstractdetector.h"
#include "intimage.h"

/********************************

This is an implementation of SURF, a few of the codes (dominant orientation estimation)
 are extracted from following project
http://code.google.com/p/opensurf1/
@author: Wan-Lei Zhao

   DSURF+RRSIFT+ VLAD*  63.0
**********************************/

class DSURF:public AbstractDetector
{
public:
    DSURF();

    static const int MaxOctaves;
    static const int PatchMag;
    static const float INITSIGMA;
    static const float SIGMA;
    static const int SCALES;
    static const float w0;
    static const float THRESH;

    static const int INTERP_KEYS;
    static const int PatchSize;
    static const float LWTHRESH;
    static const float mag; // 5
    static const int kpnum0;
    static const float PI0;
    static const float PI1;

private:
    float *featsBin;
    float *descWin;
    float thresh;
    IntImage *intImg;
    int   Octs[3][4];
    vector<vector<float> > gMat;

public:
    static IntImage *buildIntImage(Image *img0);

public:
    bool paramsCheck();
    ~DSURF();
    bool               KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn);
    void               BuildSURFOctaves(vector<vector<Image *> > &octaves, vector<vector<Image *> > &LoGOctaves);

    void               FindPeaksOctaves(vector<vector<Image *> > &HessOctaves, vector<vector<Image *> > &LoGOctaves);
    vector<KeyPoint *> FindPeaksScales(const int octave, vector<Image *> &HessImages, vector<Image *> &LoGImages);

    Image*             ScaleInitImage(Image * image);
    void               FindOrientation(vector<KeyPoint*> &keyps);

    float              InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds);
    bool               InterpKey(int x, int y, int s, vector<Image *> & DImages, float * fx, float * fy, float * fs,float *dogVal);

    bool               isSpatialPeak(Image *image, const int x, const int y);
    bool               isScalePeak(Image * aimage, Image * bimage, Image * cimage, const int x, const int y);

    static float       getAngle(float dx, float dy);

    float getDxy(const int x, const int y, const int rd0, const int ts);
    float getDxx(const int x, const int y, const int rd0, const int ts);
    float getDyy(const int x, const int y, const int rd0, const int ts);
    float getDx(const int x,  const int y, const int rd0, const int ts);
    float getDy(const int x,  const int y, const int rd0, const int ts);

    static void test();
};

#endif
