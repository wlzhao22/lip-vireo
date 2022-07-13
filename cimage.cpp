#include "cimage.h"
#include "ioimage.h"
#include "vstring.h"
#include "vmath.h"

#include <cstring>
#include <cmath>

/*************************************
*Represent color image in RGB space
*@author Wan-Lei Zhao
*@date 27.Oct.2008
**************************************/

CImage::CImage(const int width, const int height)
{
    assert(width > 0 && height > 0);
    channel  = 3;

    this->width  = width;
    this->height = height;

    pix = new float[width*height*3];
    memset(pix, 0, sizeof(float)*width*height*3);

    for(int i = 0; i < 3; i++)
    {
        rgb[i] = NULL;
    }
}

CImage::CImage(const int width, const int height, const unsigned int InitVal)
{
    assert(width > 0 && height > 0);
    channel  = 3;

    this->width  = width;
    this->height = height;

    pix = new float[3*width*height];
    memset(pix, InitVal, sizeof(float)*width*height*3);

    for(int i = 0; i < 3; i++)
    {
        rgb[i] = NULL;
    }
}

CImage::CImage(const int width, const int height, const float *new_data)
{
    assert(width > 0 && height > 0);
    channel  = 3;

    this->width  = width;
    this->height = height;

    pix = new float[3*height*width];

    memcpy(pix, new_data, 3*width * height * sizeof(float));

    for(int i = 0; i < 3; i++)
    {
        rgb[i] = NULL;
    }
}

CImage::CImage(const char *srcFn)
{
    assert(srcFn);
    int maxval;
    unsigned char *data = NULL;
    int loc, loc1, color;
    int x, y;
    channel = 3;

    if(VString::endWith(srcFn,"pgm")||VString::endWith(srcFn,"PGM"))
    {
        data = (unsigned char*)IOImage::read_pgm(srcFn, this->width, this->height, &maxval);
        color = 1;
    }
    else if(VString::endWith(srcFn,"bmp")||VString::endWith(srcFn,"BMP"))
    {
        data = (unsigned char*)IOImage::read_bmp(srcFn, this->width, this->height, &maxval, channel);
        color = channel;
    }
    else if(VString::endWith(srcFn,"jpg")||VString::endWith(srcFn,"JPG")||VString::endWith(srcFn,"JPEG")|VString::endWith(srcFn,"jpeg"))
    {
        data = (unsigned char*)IOImage::read_jpg(srcFn, this->width, this->height,&maxval, channel);
        color = channel;
    }
    else if(VString::endWith(srcFn,"ppm")||VString::endWith(srcFn,"PPM"))
    {
        data = (unsigned char*)IOImage::read_ppm(srcFn, this->width, this->height,&maxval,channel);
        color = channel;
    }
    else if(VString::endWith(srcFn,"png")||VString::endWith(srcFn,"PNG"))
    {
        data = (unsigned char*)IOImage::read_png(srcFn, this->width, this->height,&maxval,channel);
        color = channel;
    }
    else
    {
        color = 0;
        this->width = this->height = 0;
        cout<<"CImage type is unrecorganizable or suffixed with wrong type!\n";
    }

    if(data == NULL)
    {
        this->width = this->height = 0;
        cout<<"Image file "<<srcFn<<" cannot be loaded!\n"<<endl;
    }
    else
    {
        assert(this->width > 0);
        assert(this->height > 0);

        pix = new float[3*this->width*this->height];
        if(color == 3)
        {
            for (y = 0; y < this->height; y++)
            {
                for (x = 0; x < this->width; x++)
                {
                    loc = color*(y*this->width+x);
                    pix[loc] = data[loc];
                    pix[loc + 1] = data[loc + 1];
                    pix[loc + 2] = data[loc + 2];
                }
            }
        }
        else
        {
            for (y = 0; y < this->height; y++)
            {
                for (x = 0; x < this->width; x++)
                {
                    loc = color*(y*this->width+x);
                    loc1 = channel*(y*this->width+x);
                    pix[loc1] = data[loc];
                    pix[loc1 + 1] = data[loc];
                    pix[loc1 + 2] = data[loc];
                }
            }
        }
        delete [] data;
    }

    for(int i = 0; i < 3; i++)
    {
        rgb[i] = NULL;
    }
}

Image *CImage::getGrayImage()
{
    Image *grayImg = new Image(this->width,this->height);
    int i = 0, j = 0;
    unsigned int loc = 0;
    float val = 0.0f;

    for (j = 0; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            loc = 3*(j*width+i);
            val = (pix[loc] + pix[loc+1]+pix[loc+2])/3.0;
            grayImg->setPixel(i,j,val);
        }
    }

    return grayImg;
}

CImage *CImage::getBlueImage()
{
    CImage *blueImg = new CImage(this->width,this->height);
    int i = 0, j = 0;
    unsigned int loc = 0;
    float val = 0.0f;

    for (j = 0; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            loc = 3*(j*width+i);
            val = pix[loc+2];
            blueImg->setPixelR(i, j, val);
        }
    }

    return blueImg;
}

void CImage::setPixel(const int x,const int y,const float*color)
{
    assert(color);
    if(y >= this->height ||y < 0)
        return ;

    if(x < 0|| x >= this->width)
        return ;

    pix[3*(y*width+x)]   = color[0];
    pix[3*(y*width+x)+1] = color[1];
    pix[3*(y*width+x)+2] = color[2];
}

void CImage::setPixelR(const int x, const int y, const float val)
{
    if(y >= this->height ||y < 0)
        return ;

    if(x < 0 || x >= this->width)
        return ;

    pix[3*(y*width+x)] = val;
}

void CImage::setPixelG(const int x,const int y, const float val)
{
    if(y >= this->height || y < 0)
        return ;

    if(x < 0 || x >= this->width)
        return ;
    pix[3*(y*width+x)+1] = val;
}

void CImage::setPixelB(const int x,const int y, const float val)
{
    if(y >= this->height||y < 0)
        return ;

    if(x < 0||x >= this->width)
        return ;
    pix[3*(y*width+x)+2] = val;
}

void CImage::getPixel(const int x, const int y,float* RGB)
{
    assert(RGB);

    if(y >= this->height || y < 0)
        return ;

    if(x < 0 || x >= this->width)
        return ;

    int loc = 3*(y*width+x);

    RGB[0] = pix[loc];
    RGB[1] = pix[loc+1];
    RGB[2] = pix[loc+2];

    return ;
}

void CImage::getPixelLUV(const int x,const int y, float *LUV)
{
    assert(LUV);

    if(y >= this->height||y < 0)
        return ;

    if(x < 0 ||x >= this->width)
        return ;

    int loc = y*width+x;

    if(rgb[0] == NULL)
    {
        LUV[0] = 0;
        LUV[1] = 0;
        LUV[2] = 0;
    }
    else
    {
        LUV[0] = rgb[0][loc];
        LUV[1] = rgb[1][loc];
        LUV[2] = rgb[2][loc];
    }
    return ;
}


CImage * CImage::halfSizeImage(CImage *im)
{
    int w = im->width/2;
    int h = im->height/2;

    CImage * nim = new CImage(w, h);
    float RGB[3] = {0};

    for (int j = 0; j < h; j++)
    {
        for (int i = 0; i < w; i++)
        {
            im->getPixelBI(i*2, j*2, RGB);
            nim->setPixelR(i, j, RGB[0]);
            nim->setPixelG(i, j, RGB[1]);
            nim->setPixelB(i, j, RGB[2]);
        }
    }

    return nim;
}

CImage * CImage::doubleSizeImage(CImage *srcImg)
{
    int w = srcImg->width*2;
    int h = srcImg->height*2;
    int i = 0, j = 0;
    float RGB[3] = {0};
    float RGB1[3] = {0};

    CImage * dstImg = new CImage(w, h);

    for (j = 0; j < h; j++)
    {
        for (i = 0; i < w; i++)
        {
            srcImg->getPixel(i/2, j/2, RGB);
            dstImg->setPixelR(i, j, RGB[0]);
            dstImg->setPixelG(i, j, RGB[1]);
            dstImg->setPixelB(i, j, RGB[2]);
        }
    }

    /**
      A B C
      E F G
      H I J

      pixels A C H J are pixels from original image
      pixels B E G I F are interpolated pixels
    **/

    // interpolate pixels B and I
    for (j = 0; j < h; j += 2)
        for (i = 1; i < w - 1; i += 2)
        {
            srcImg->getPixel(i/2, j/2, RGB);
            srcImg->getPixel(i/2 + 1, j/2,RGB1);
            dstImg->setPixelR(i, j,(RGB[0]+RGB1[0])/ 2.0);
            dstImg->setPixelG(i, j,(RGB[1]+RGB1[1])/ 2.0);
            dstImg->setPixelB(i, j,(RGB[2]+RGB1[2])/ 2.0);
        }

    // interpolate pixels E and G
    for (j = 1; j < h - 1; j += 2)
        for (i = 0; i < w; i += 2)
        {
            srcImg->getPixel(i/2, j/2,RGB);
            srcImg->getPixel(i/2, j/2+1,RGB1);
            dstImg->setPixelR(i, j,(RGB[0]+RGB1[0])/ 2.0);
            dstImg->setPixelG(i, j,(RGB[1]+RGB1[1])/ 2.0);
            dstImg->setPixelB(i, j,(RGB[2]+RGB1[2])/ 2.0);
        }

    //interpolate pixel F
    //interpolate pixels E and G
    float RGB2[3] = {0};
    float RGB3[3] = {0};

    for (j = 1; j < h - 1; j += 2)
        for (i = 1; i < w - 1; i += 2)
        {
            srcImg->getPixel(i/2, j/2,RGB);
            srcImg->getPixel(i/2 + 1, j/2,RGB1);
            srcImg->getPixel(i/2, j/2+1,RGB2);
            srcImg->getPixel(i/2 + 1, j/2+1,RGB3);

            dstImg->setPixelR(i, j,(RGB[0]+RGB1[0]+RGB2[0]+RGB3[0])/ 4.0);
            dstImg->setPixelG(i, j,(RGB[1]+RGB1[1]+RGB2[1]+RGB3[1])/ 4.0);
            dstImg->setPixelB(i, j,(RGB[2]+RGB1[2]+RGB2[2]+RGB3[2])/ 4.0);
        }
    return dstImg;
}

CImage * CImage::clone()
{
    CImage * im = new CImage(width, height, pix);
    return im;
}

void CImage::getPixelBI(float col, float row, float *RGB)
{
    int irow = 0, icol = 0;
    float rfrac = 0.0f, cfrac = 0.0f;
    float row1 = 0, row2 = 0;

    irow = (int) row;
    icol = (int) col;

    if (irow < 0 || irow >= height
            || icol < 0 || icol >= width)
        return ;

    if (row > height - 1)
        row = height - 1;

    if (col > width - 1)
        col = width - 1;

    rfrac = 1.0 - (row - (float) irow);
    cfrac = 1.0 - (col - (float) icol);

    int row_loc =  irow*width;

    if (cfrac < 1)
    {
        row1 = cfrac * pix[row_loc+icol] + (1.0 - cfrac) * pix[row_loc+icol + 1];
    }
    else
    {
        row1 = pix[row_loc+icol];
    }

    row_loc = (irow+1)*width;
    if (rfrac < 1)
    {
        if (cfrac < 1)
        {
            row2 = cfrac * pix[row_loc+icol] + (1.0 - cfrac) * pix[row_loc+icol + 1];
        }
        else
        {
            row2 = pix[row_loc+icol];
        }
    }
    RGB[0] = row1+row2;

    return ;
}

bool CImage::RGBtoLUV()
{
    int i = 0, j = 0, loc = 0;
    for(i = 0; i < channel; i++)
    {
        if(rgb[i] == NULL)
        {
            rgb[i] =  new float[this->width*this->height];
        }
    }

    float L, M, S, l, m, s;
    float r, g, b;
    float RGB[3];

    for(i = 0; i < this->height; i++)
    {
        for(j = 0; j < this->width; j++)
        {
            this->getPixel(j, i, RGB);
            r = RGB[0];
            g = RGB[1];
            b = RGB[2];
            loc = width*i+j;

            l = 0.3811f*r + 0.5783f*g + 0.0402f*b;
            m = 0.1967f*r + 0.7244f*g + 0.0782f*b;
            s = 0.0241f*r + 0.1288f*g + 0.8444f*b;

            if(l <= 0)
            {
                l = 2.0f;
            }

            if(m <= 0)
            {
                m = 2.0f;
            }

            if(s <= 0)
            {
                s = 2.0f;
            }

            L = VMath::lgx(l, 10);
            M = VMath::lgx(m, 10);
            S = VMath::lgx(s, 10);

            rgb[0][loc] = 0.57735f*L  + 0.57735f*M   + 0.57735f*S;
            rgb[1][loc] = 0.408248f*L + 0.408248f*M   + -0.816496f*S;
            rgb[2][loc] = 0.707107f*L + -0.707107f*M + 0.0f*S;
        }
    }
    return true;
}

CImage::~CImage()
{
    for(int i = 0; i < 3; i++)
    {
        if(rgb[i] != NULL)
        {
            delete [] rgb[i];
            rgb[i] =  NULL;
        }
    }
}

void CImage::test()
{
    const char *srcfn = "/udd/wazhao/datasets/inria/bricks/ppm/img_6592.ppm";
    const char *dstfn = "/udd/wazhao/datasets/inria/img_6592.jpg";

    CImage *myimg = new CImage(srcfn);
    CImage *blue  = myimg->getBlueImage();
    Image  *gray  = myimg->getGrayImage();
    gray->save(dstfn);
    ///blue->save(dstfn);
    ///myimg->save(dstfn);
    delete myimg;
    delete blue;
    delete gray;
}
