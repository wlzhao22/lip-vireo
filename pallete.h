#ifndef PALLETE_H
#define PALLETE_H

#include "abstractimage.h"
#include "keypoint.h"
#include "cimage.h"


enum VIREO_COLOR {_vireo_white = 0, _vireo_red = 1, _vireo_blue, _vireo_purple, _vireo_black, _vireo_yellow};

class Pallete
{
private:
    AbstractImage *image;
    static const unsigned int WHITE;;
    static const unsigned int RED;
    static const unsigned int YELLW;
    static const unsigned int PURPL;
    static const unsigned int BLACK;
    static const float  PI0;
    static const float  PI3;

public:
    static bool   fetchColor(VIREO_COLOR color, float mycolor[3]);
    static void   drawEllipse(CImage *img, KeyPoint *crnt_pt, VIREO_COLOR color);
    static void   drawCorners(CImage *view, vector<HCPoint*> curves, const unsigned int lbound);
    static void   lineto(CImage *img, const unsigned int x0, const unsigned int y0, const unsigned int x1, const unsigned int y1, VIREO_COLOR mycolor);
    static void   drawTriangle(CImage *Img, vector<HCPoint*> &curves, VIREO_COLOR fcolor);
    static void   display(CImage *srcimg, const unsigned int x, const unsigned int y, CornType mytype);
    static void   buildView(const char *srcfn, const char *dstfn, vector<KeyPoint*> &kps, VIREO_COLOR fcolor, KeyPtType flt);
    static void   buildEdgeView(const unsigned char* edgeimg, const int width, const int height, const char *dstfn);
    static void   InitView(CImage *img, VIREO_COLOR bgcolor);
    static void   test();
};

#endif
