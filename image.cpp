/********************************************************************

shared from Yan Ke, but extensive changes have been made by Wanlei Zhao

*********************************************************************/

#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "vstring.h"
#include "ioimage.h"
#include "filter.h"
#include "image.h"

Image::Image(int width, int height)
{
    assert(width > 0 && height > 0);
    this->width  = width;
    this->height = height;
    channel      = 1;

    pix = new float[width*height];
    memset(pix, 0, sizeof(float)*width*height);
}

Image::Image(int width, int height, const float * new_data)
{
    assert(width > 0 && height > 0);

    this->width = width;
    this->height = height;
    channel = 1;

    pix = new float[height*width];

    memcpy(pix, new_data, width * height * sizeof(float));
}

Image::Image(const char *srcFn, const float norm)
{
    assert(srcFn);
    int maxval = 0;
    unsigned char *data = NULL;
    channel = 1;

    if(VString::endWith(srcFn, "pgm")||VString::endWith(srcFn, "PGM"))
    {
        data = (unsigned char*)IOImage::read_pgm(srcFn, this->width, this->height, &maxval);
    }
    else if(VString::endWith(srcFn, "bmp")||VString::endWith(srcFn, "BMP"))
    {
        data = (unsigned char*)IOImage::read_bmp(srcFn, this->width, this->height, &maxval,channel);
    }
    else if(VString::endWith(srcFn, "jpg")||VString::endWith(srcFn, "JPG")||VString::endWith(srcFn, "JPEG")|VString::endWith(srcFn,"jpeg"))
    {
        data = (unsigned char*)IOImage::read_jpg(srcFn, this->width, this->height,&maxval,channel);
    }else if(VString::endWith(srcFn, "ppm")||VString::endWith(srcFn, "PPM"))
    {
        data = (unsigned char*)IOImage::read_ppm(srcFn, this->width, this->height,&maxval,channel);
    }else if(VString::endWith(srcFn, "png")||VString::endWith(srcFn,"PNG"))
    {
        data = (unsigned char*)IOImage::read_png(srcFn, this->width, this->height, &maxval, channel);
    }
    else
    {
        cout<<srcFn<<endl;
        this->width = this->height = 0;
        cout<<"Image type is unrecorganizable or suffixed with wrong type!\n";
    }
    if(data == NULL)
    {
        cout<<"Image file "<<srcFn<<" cannot be loaded!\n"<<endl;
    }else{
      assert(this->width  > 0);
      assert(this->height > 0);

      pix = new float[this->width*this->height];
      for (int y = 0; y < this->height; y++)
      {
         for (int x = 0; x < this->width; x++)
         {
             pix[y*this->width+x] = data[y*this->width+x]/norm;
         }
     }
     delete [] data;
    }
}

Image::Image(const char *srcFn)
{
    assert(srcFn);
    int maxval = 0;
    unsigned char *data = NULL;
    channel = 1;

    if(VString::endWith(srcFn, "pgm")||VString::endWith(srcFn, "PGM"))
    {
        data = (unsigned char*)IOImage::read_pgm(srcFn, this->width, this->height, &maxval);
    }
    else if(VString::endWith(srcFn, "bmp")||VString::endWith(srcFn,"BMP"))
    {
        data = (unsigned char*)IOImage::read_bmp(srcFn, this->width, this->height, &maxval,channel);
    }
    else if(VString::endWith(srcFn, "jpg")||VString::endWith(srcFn, "JPG")||VString::endWith(srcFn, "JPEG")||VString::endWith(srcFn, "jpeg"))
    {
        data = (unsigned char*)IOImage::read_jpg(srcFn, this->width, this->height,&maxval, channel);
    }
    else if(VString::endWith(srcFn, "ppm")||VString::endWith(srcFn, "PPM"))
    {
        data = (unsigned char*)IOImage::read_ppm(srcFn, this->width, this->height, &maxval, channel);
    }
    else if(VString::endWith(srcFn, "png")||VString::endWith(srcFn, "PNG"))
    {
        data = (unsigned char*)IOImage::read_png(srcFn, this->width, this->height, &maxval, channel);
    }
    else
    {
        this->width = this->height = 0;
        cout<<"Image type is unrecorganizable or suffixed with wrong type!\n";
    }
    if(data == NULL)
    {
        this->width = this->height = 0;
        cout<<"Image file "<<srcFn<<" cannot be loaded!\n"<<endl;
    }else{
        assert(this->width  > 0);
        assert(this->height > 0);
        pix = new float[this->width*this->height];
        for (int y = 0; y < this->height; y++)
        {
            for (int x = 0; x < this->width; x++)
            {
                pix[y*this->width+x] = data[y*this->width+x];
            }
        }
        delete [] data;
    }
}

void Image::multiply(const float val)
{
    int len = width*height;
    for(int i = 0; i < len; i++)
    {
        this->pix[i] = this->pix[i]*val;
    }
}

void Image::exp(const int val)
{
    int len = this->width*this->height;
    for(int i = 0; i < len; i++)
    {
        this->pix[i] = this->pix[i]*this->pix[i];
    }
}

Image * Image::halfSizeImage(Image *im)
{
    int w = im->width/2;
    int h = im->height/2;
    Image * nim = new Image(w, h);

    for (int j = 0; j < h; j++)
    {
        for (int i = 0; i < w; i++)
        {
            nim->setPixel(i, j, im->getPixelBI(i*2, j*2));
        }
    }

    return nim;
}

Image * Image::doubleSizeImage(Image *srcImg)
{
    int w = srcImg->width*2;
    int h = srcImg->height*2;
    int i,j;

    Image * dstImg = new Image(w, h);

    for (j = 0; j < h; j++)
    {
        for (i = 0; i < w; i++)
        {
            dstImg->setPixel(i, j, srcImg->getPixel(i/2, j/2));
        }
    }

    // interpolate pixels B and I
    for (j = 0; j < h; j += 2)
        for (i = 1; i < w - 1; i += 2)
        {
            dstImg->setPixel(i, j, (srcImg->getPixel(i/2, j/2) + srcImg->getPixel(i/2 + 1, j/2)) / 2.0);
        }


    // interpolate pixels E and G
    for (j = 1; j < h - 1; j += 2)
        for (i = 0; i < w; i += 2)
        {
            dstImg->setPixel(i, j, (srcImg->getPixel(i/2, j/2) + srcImg->getPixel(i/2, j/2 + 1)) / 2.0);
        }

    // interpolate pixel F
    // interpolate pixels E and G
    for (j = 1; j < h - 1; j += 2)
        for (i = 1; i < w - 1; i += 2)
        {
            dstImg->setPixel(i, j, (srcImg->getPixel(i/2, j/2) + srcImg->getPixel(i/2, j/2 + 1)
                                    + srcImg->getPixel(i/2 + 1, j/2) + srcImg->getPixel(i/2 + 1, j/2 + 1))/4.0);
        }
    return dstImg;
}

void Image::sub(Image * im2, Image * dstImg)
{
    assert(im2);
    assert(dstImg);
    assert(width  == im2->width && width   == dstImg->width);
    assert(height == im2->height && height == dstImg->height);

    for (int j = 0; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            dstImg->setPixel(i, j, getPixel(i, j) - im2->getPixel(i, j));
        }
    }
    return ;
}

Image * Image::clone()
{
    Image *im = new Image(width, height, pix);
    return im;
}

Image *Image::getPatch(const int x,const int y,const int win,float *umatrix)
{
    int wc = win/2;
    Image *winImg = new Image(win,win);

    float dx,dy,xl,yl,xr,yr;
    int xp, yp, row, col, irow, icol;
    float tmp;

    for(row = -wc; row <= wc; row++)
    {
        for(col = -wc; col <= wc; col++)
        {
            xr = umatrix[0]*col + umatrix[1]*row;
            yr = umatrix[2]*col + umatrix[3]*row;
            xl = floor(xr);
            yl = floor(yr);
            dx = xr - xl;
            dy = yr - yl;
            xp = (int)round(xr);
            yp = (int)round(yr);
            xp = x + xp;
            yp = y + yp;
            /**/
            tmp = this->getPixel(xp,yp)*(1-dx)*(1-dy); //xp,yp
            tmp += this->getPixel(xp+1,yp)*(1-dy)*dx; //xp+1,yp
            tmp += this->getPixel(xp,yp+1)*(1-dx)*dy; //xp,yp+1
            tmp += this->getPixel(xp+1,yp+1)*dx*dy;   //xp+1,yp+1
            /**/
            irow = row + wc;
            icol = col + wc;
            winImg->setPixel(icol,irow,tmp);
        }
    }

    return winImg;
}


float Image::getConValWidth(const int x, const int y, vector<float> &dkern)
{
    /**
    *blur in the x direction with dkern, and smooth in y direction with ikern
    *sharp denotes pixel location where the operation made on
    **/

    int dscale = dkern.size();
    int dc = dscale/2;
    float pixel = 0.0f, tmp = 0.0f;
    int j = 0;
    pixel = 0.0;
    for (j = -dc; j <= dc; j++)
    {
        tmp = this->getPixel(x+j,y);
        pixel += dkern[dc+j] * tmp;
    }
    return pixel;
}

float Image::getConValHeight(const int x, const int y, vector<float> &dkern)
{
    /**
    *blur in the x direction with dkern, and smooth in y direction with ikern
    *sharp denotes pixel location where the operation made on
    **/

    int j = 0, dscale = dkern.size();
    int dc = dscale/2;
    float pixel = 0, tmp = 0;
    pixel = 0.0;
    for (j = -dc; j <= dc; j++)
    {
        tmp = this->getPixel(x,y+j);
        pixel += dkern[dc+j] * tmp;
    }
    return pixel;
}

void Image::test()
{
    const char *srcfn2  = "/home/wlzhao/datasets/inria/boat/img0.pgm";
    const char *dstfn2  = "/home/wlzhao/datasets/inria/img0.jpg";

    Image *myimg = new Image(srcfn2);
    cout<<myimg->width<<"\t"<<myimg->height<<endl;
    myimg->save(dstfn2);
    delete myimg;
}

