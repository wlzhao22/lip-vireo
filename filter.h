#ifndef FILTER_H
#define FILTER_H

#include <vector>
#include "image.h"

#define GAUSSKERN 4
#define PI 3.14159265358979323846
#define PI_2 1.57079633
#define PI_23 4.71238898
#define PI2 6.28318531

using namespace std;
/***************************************************************************
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**************************************************************************/


class Filter
{

public:
    static vector<float>          GaussianKernel1D(const float sigma);
    static void                   GaussianKernel1D(const float sigma,const float ratio, vector<float> &kern);
    static vector<float>          GaussianDxKernel1D(const float sigma);
    static vector<float>          GaussianDxxKernel1D(const float sigma);
    static vector<vector<float> > GaussianDxyKernel2D(const float sigma);
    static vector<vector<float> > GaussianKernel2D(const float sigma);
    static vector<vector<float> > GaussianKernel2D(const float sigma, const int radius);
    static void                   GaussianKernel2D(const float sigma, vector<vector<float> > &G2);

    static void   Convolve1DWidth(vector<float> & kern, Image * src, Image * dst);
    static void   Convolve1DHeight(vector<float> & kern, Image * src, Image * dst);
    static void   BlurImage(Image * src, Image * dst, const float sigma);
    static void   BlurImage(Image * src, const float sigma);
    static bool   multiply(Image *img1,Image *img2,Image *dstimg);
    static float *GaussianKernel2D(const float sigma,int &size);
    static bool   getDxx(const float sigma,Image *srcimg,Image *dstimg);
    static bool   getDxy(const float sigma,Image *srcimg,Image *dstimg);
    static Image* getDxx(const Image *srcimg);
    static Image* getDyy(const Image *srcimg);
    static Image* getDxy(const Image *srcimg);
    static Image* getDx(const Image  *srcimg);
    static Image* getDy(const Image  *srcimg);
    static void   nonMaxSuppr(Image  *srcImg);

    static bool   GetPixOrientation(const float x, const float y, const Image * image, float &m,  float &theta);
    static bool   GetPixCurl(const float x, const float y, const Image * image, float &curl);
    static bool   GetPixOrientation(const int x, const int y, const Image * image, float &m, float &m1, float &theta);
    static void   SmoothHistogram(float *hist,const int n);
    static void   SmoothHistogram(float *hist, const int n,const int start);
    static void   clear2Dvector(vector<vector<float> > &G2);
    static float  Average(float *hist,int n);
    static float  fast_expn(const float x);
    static void   normalizeMat(vector<vector<float> > & mat);
    static void   normalizeVec(vector<float> & vec);
    static void   getLoG();
    static void   printMat(vector<vector<float> > &mat);
    static void   ScaleImage(Image *src);
    static Image *ScaleImage(Image *srcImg,const int zoomout,const float sigma,const float sigma0);

    static void test();

};

#endif
