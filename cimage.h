#ifndef CIMAGE_H
#define CIMAGE_H

#include "abstractimage.h"
#include "image.h"

class CImage :public AbstractImage
{
private:
    float *rgb[3];

public:
    CImage(const int width, const int height);

    CImage(const int width, const int height, const unsigned int InitVal);
    CImage(const int width, const int height, const float *new_data);

    CImage(const char * filename);
    ~CImage();

    void setPixelR(const int x, const int y, const float val);

    void setPixelG(const int x, const int y, const float val);

    void setPixelB(const int x, const int y, const float val);

    void setPixel(const int x,const int y,const float*color);

    void getPixel(const int x,const int y,float* RGB);
    void getPixelBI(float col,float row,float *RGB);

    void get2DConVal(const int x,const int y,float pixel[],vector<float> &kern);

    Image *getGrayImage();
    CImage *getBlueImage();

    bool RGBtoLUV();
    void getPixelLUV(const int x,const int y, float *LUV);

    static CImage *doubleSizeImage(CImage *srcImg);
    static CImage *halfSizeImage(CImage *srcImg);
    CImage *clone();

    static void test();
};
#endif
