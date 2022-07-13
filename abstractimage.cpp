#include "abstractimage.h"
#include "vstring.h"
#include "ioimage.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

void AbstractImage::save(const char *imgFn)
{
    assert(imgFn);
    unsigned char maxval = 0;

    if(VString::endWith(imgFn, ".pgm")||VString::endWith(imgFn, ".PGM"))
    {
        IOImage::write_pgm(imgFn, width, height, pix, maxval, NULL, channel);
    }
    else if(VString::endWith(imgFn, ".bmp")||VString::endWith(imgFn, ".BMP"))
    {
        IOImage::write_bmp(imgFn, width, height, pix, channel);
    }
    else if(VString::endWith(imgFn, ".ppm")||VString::endWith(imgFn, ".PPM"))
    {
        IOImage::write_ppm(imgFn, pix, width, height, channel);
    }
    else if(VString::endWith(imgFn, ".jpg")||VString::endWith(imgFn, ".JPG"))
    {
        IOImage::write_jpg(imgFn, pix, width, height, channel, 100);
    }
}

float AbstractImage::get2DConVal(const int x, const int y, vector<float> &kern)
{
    int ic = kern.size()/2;
    int i, j, row, col;

    float pixel = 0, gray = 0;
    gray = 0;
    for(i = -ic; i <= ic; i++)
    {
        row = y + i;
        pixel = 0.0;
        for (j = -ic; j <= ic; j++)
        {
            col = x + j;
            this->getPixel(col,row, &gray); //xp,yp
            pixel += kern[ic+j] * gray;
        }
        gray = gray  + pixel*kern[ic+i];
    }
    return gray;
}

void AbstractImage::get2DConVal(const int x, const int y, float pixVals[3], vector<float> &kern)
{
    int ic = kern.size()/2;
    int i, j, row, col;
    assert(pixVals);
    float RGB[3], abc[3];
    pixVals[0] = pixVals[1] = pixVals[2] = 0.0f;

    for(i = -ic; i <= ic; i++)
    {
        row = y + i;
        abc[0] = abc[1] = abc[2] = 0.0f;
        for (j = -ic; j <= ic; j++)
        {
            col = x + j;
            getPixel(col,row,RGB); //xp,yp
            abc[0] += kern[ic+j] * RGB[0];
            abc[1] += kern[ic+j] * RGB[1];
            abc[2] += kern[ic+j] * RGB[2];
        }
        pixVals[0] = pixVals[0] + abc[0]*kern[ic+i];
        pixVals[1] = pixVals[1] + abc[1]*kern[ic+i];
        pixVals[2] = pixVals[2] + abc[2]*kern[ic+i];
    }
    return ;
}

float AbstractImage::getPixel(const int x, const int y)
{
    float RGB[3];
    float gray = 0;
    if(this->channel == 3)
    {
        getPixel(x, y, RGB); //xp,yp
        gray += RGB[0];
        gray += RGB[1];
        gray += RGB[2];
        gray  = gray/3.0f;
    }
    else
    {

        if(x >= 0 && x < width && y >= 0 && y < height)
        {
            unsigned int loc = y*width+x;
            return pix[loc];
        }
        else
        {
            return 0;
        }
    }
    return 0;
}

void AbstractImage::setPixel(const int x, const int y, const float val)
{
    float RGB[3];
    if(this->channel == 3)
    {
        RGB[0] = RGB[1] = RGB[2] = val;
        setPixel(x, y, RGB); //xp,yp
    }
    else
    {
        if(x >= 0 && x < width && y >= 0 && y < height)
        {
            unsigned int loc = y*width+x;
            pix[loc] = val;
        }
        else
        {
            return ;
        }
    }
}

float AbstractImage::getPixelBI(float col, float row) const
{
    int irow = 0, icol = 0, ploc = 0;
    float dx = 0, dy = 0;
    float val = 0.0f;

    irow = (int)floor(row);
    icol = (int)floor(col);
    int hbound = height - 1;
    int wbound = width - 1;

    if(irow >= hbound || irow < 0)
        return 0;

    if(icol >= wbound || icol < 0)
        return 0;

    dy = row -  irow;
    dx = col -  icol;

    ploc =  irow*this->width + icol;
    val =  (1 - dy)*(1 - dx)*pix[ploc];

    if(icol < wbound)
        val += (1 - dy)* dx*pix[ploc + 1];

    if(irow < hbound)
    {
        ploc  = ploc + this->width;
        val  += (1 - dx)*dy*pix[ploc];
        if(icol < wbound)
        {
            val += dx*dy*pix[ploc + 1];
        }
    }

    return val;
}

