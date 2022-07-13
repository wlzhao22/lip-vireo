#ifndef INTIMAGE_H
#define INTIMAGE_H

#include <iostream>
#include <cassert>
#include <cstring>

using namespace std;

class IntImage
{
    public:
        IntImage()
        {
            pix = NULL;
        }

        IntImage(const int w, const int h)
        {
            assert(w > 0);   assert(h > 0);
            width =  w;      height = h;
            pix = new double[w*h];
            memset(pix, 0, sizeof(double)*w*h);
        }

        double *pix;
        int width, height;
        double getPixel(const int x, const int y);
        void   setPixel(const int x, const int y, const double val);
        float  boxIntegral(const int x0, const int y0, const int cols, const int rows);

        virtual ~IntImage()
        {
            if(pix != NULL)
            delete [] pix;
        }
};

#endif // INTIMAGE_H
