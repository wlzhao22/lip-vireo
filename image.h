/**
*@author Wan-Lei Zhao
*@version 1.30
*All rights reserved by Wanlei Zhao
*
**/

#ifndef IMAGE_H
#define IMAGE_H

#include "abstractimage.h"
#include <vector>

using namespace std;

class Image: public AbstractImage
{

public:
    Image(const int width, const int height);
    Image(const int width, const int height, const float *new_data);
    Image(const char *srcFn, const float norm);
    Image(const char *srcFn);
    virtual ~Image(){}

    /**
    *
    *pix = pix * val;
    *  val multiply with each pixel in that image
    */

    void multiply(const float val);

    /**
    *
    *make each pixel exponentialized with a cerntain value
    *pix = pix^val
    */
    void exp(const int val);
    inline void setPixel(const int x, const int y, const float val)
    {
        if(y >= this->height || y < 0)
            return ;

        if(x < 0||x >= this->width)
            return ;
        pix[y*width+x] = val;
    }

    inline void setPixel(const int x, const int y, const float *gray)
    {
        assert(gray);

        if(y >= this->height||y < 0)
            return ;

        if(x < 0|| x >= this->width)
            return ;

        pix[y*width+x] = *gray;
    }

    inline float getPixel(const int x, const int y) const
    {
        if(y>=this->height||y<0)
        {
            return 0;
        }

        if(x<0||x>=this->width)
        {
            return 0;
        }
        return pix[y*width+x];
    }

    inline void  getPixel(const int x, const int y, float *val)
    {
        if(y>=this->height||y<0)
        {
            *val = 0;
            return ;
        }

        if(x<0||x>=this->width)
        {
            *val = 0;
            return ;
        }
        *val = pix[y*width+x];
    }
    /**/


    /**
    *get halfsized image of the original image
    */

    static Image * halfSizeImage(Image *im);

    /**
    * get double sized image of the orignial image
    */

    static Image *doubleSizeImage(Image *srcImg);
    /**
    * Clones a copy of the image.
    * @return Cloned copy.
    */
    Image * clone();

    /**
    *substract image with im2 and store it to 'dst' image
    **/
    void sub(Image * im2, Image * dst);

    float getConValWidth(const int x,  const int y, vector<float> &dkern);
    float getConValHeight(const int x, const int y, vector<float> &dkern);
    Image *getPatch(const int x, const int y, const int win, float *umatrix);

    static void test();
};


#endif
