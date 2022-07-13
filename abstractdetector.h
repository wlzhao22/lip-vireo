#ifndef ABSTRACTDETECTOR_H
#define ABSTRACTDETECTOR_H

#include "abstractdescriptor.h"
#include "idetector.h"
#include "keypoint.h"
#include "intimage.h"
#include "cimage.h"
#include "config.h"
#include "image.h"


/********************************************************************
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**********************************************************************/

#include <vector>
#include <map>

using namespace std;

enum GBlur {_Conv2 = 0, _Dx_, _Dy_, _Dxx_, _Dyy_, _Dxy_};


#define _Dx(img, x, y)  ((img->getPixel(x+1, y) - img->getPixel(x-1, y))/2.0);
#define _Dy(img, x, y)  ((img->getPixel(x, y+1) - img->getPixel(x, y - 1))/2.0);
#define _Dxx(img, x, y) (img->getPixel(x+1, y)  + img->getPixel(x-1, y) - 2 * img->getPixel(x, y));
#define _Dyy(img, x, y) (img->getPixel(x, y+1)  + img->getPixel(x, y-1) - 2 * img->getPixel(x, y));
#define _Dxy(img, x, y) ((img->getPixel(x+1, y+1) + img->getPixel(x-1, y-1)- img->getPixel(x-1, y+1) - img->getPixel(x+1, y-1))/4.0);

class AbstractDetector: public IDetector
{

protected:

    bool INIT;
    bool RESIZE;
    static const int imgSzBound;

    map<string, const char*> paras;
    vector<KeyPoint*> kps;
    vector<KeyPoint*> leveli_kps;
    Image *crntimg;
    SelOPT sel_option;
    DESC mydesc;

    ///for descriptor
    AbstractDescriptor *myDescriptor;

    ///for dominant orientation estimation
    IntImage *intImg;
    vector<vector<float> > gauss25;

    ///out properties, e.g., scale and rotation for each detected points

    ///for orientation detection
    /**/
    static const int winSize;
    static const int MaxOctaves;
    static const int kp_num0;
    int       fix_kp_numb;

    ///for detector
    static const int   SCALESPEROCTAVE;
    static const int   BORDER;
    static const int   KNN;
    static const int   nOrient0;
    static const int   nScale0;
    static const float minDist;
    static const int   EPS;
    static const int   KN0;  //for kernel mapping
    static const int   Ki0;  //for kernel mapping

protected:
    static const float Ec;
    static const float El;
    static const float alpha;
    static const float PI0;
    static const float PI4;
    static const float deltaSc0;
    static const float deltaSc1;
    Detector DETECTOR;
    char detectorsID[15][10];

    char tm_logfn[FNLEN];

protected:

    void           BuildOctaves(Image * image, GBlur blur_opt, vector<vector<Image *> > &Ggoctaves, const float _SIGMA0,const int NUM_SCALES0);
    vector<Image*> BuildScales(Image * image,  GBlur blur_opt, const float _SIGMA0,const int NUM_SCALES0);

    float getDx(const int x, const int y, const int rd0, const int ts);
    float getDy(const int x, const int y, const int rd0, const int ts);

public:
    AbstractDetector();
    bool        Init(const char *scriptfn, const char *descopt);
    void        FindOrientation(vector<KeyPoint*> &keyps);


    bool        (AbstractDetector::*saveKpts)(vector<KeyPoint*> &kplst, const int kpnum, const char *dstfn,
                         const float resize_rate0, bool AFF_OUT);

    bool        saveKpVGG(vector<KeyPoint*> &kplst, const int kpnum, const char *dstfn,
                         const float resize_rate0, bool AFF_OUT);

    bool        saveKpVIREO(vector<KeyPoint*> &kplst, const int kpnum, const char *dstfn,
                         const float resize_rate0, bool AFF_OUT);

    virtual bool paramsCheck() = 0;
    virtual bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn) = 0;

    vector<KeyPoint*> getKeypoints()
    {
        return this->kps;
    }

    static int   findFlip(const float thetas[], const int dbin, const int nbin, float &ratio);
    static bool  topkSelect(vector<KeyPoint*> &kplst, const int topk);
    static bool  topkEqDnSelect(vector<KeyPoint*> &kplst, const int w,const int h,const int topk);
    static bool  mapPCA(const float *feats, const float *pcamat, vector<float> &means,
                        vector<float> &var, float *pcaft, const int dim, const int pdim);
    static bool  proj(const float *feats, const float *prjMat, float *prjFt, const int dim, const int pdim);

    bool         detectKeyPoints(const char *srcdir,  const char *kpdir, const char *drdir,
                                 const char *destdir, const char *dvdir, bool ONE_MISSION);
    void         printAVG();

public:
    static float getAngle(float dx, float dy);
    static IntImage *buildIntImage(Image *img0);
    static bool  releaseScales(vector<Image*> &scales);
    static bool  releaseKpList(vector<KeyPoint *> &keyps);

    ~AbstractDetector();
};
#endif
