#ifndef ABSTRACTIMAGE_H
#define ABSTRACTIMAGE_H

#include <iostream>
#include <cassert>
#include <vector>

using namespace std;

class AbstractImage
{

public:
    AbstractImage()
    {
        pix = NULL; width = 0; height = 0;
    }
    int channel, width, height;
    float *pix;

public:
    virtual void getPixel(const int x,const int y,float *RGB) = 0;
    virtual void setPixel(const int x,const int y, const float *color) = 0;

    float get2DConVal(const int x, const int y, vector<float> &kern);
    void  get2DConVal(const int x, const int y, float pixVals[3], vector<float> &kern);

    float getPixel(const int x, const int y);
    bool  isActive()
    {
       if(this->pix != NULL && this->width > 0 && this->height > 0)
       {
           return true;
       }else{
           return false;
       }
    }
    void  setPixel(const int x, const int y, const float val);
    /**
    *get pixel value which has been bilinearly interpolated
    */
    float getPixelBI(float x,  float y) const;

    void save(const char *filename);
    virtual ~AbstractImage()
    {
        delete [] pix;
        pix = NULL;
    }

};
#endif
