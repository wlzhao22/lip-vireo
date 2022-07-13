#ifndef IMAGESYMMMAT_H
#define IMAGESYMMMAT_H

#include <cassert>
#include <cstring>
#include <vector>

using namespace std;

class ImageSymmMat
{

private:
    int w, h, size;
    float *pmat;

public:
    ImageSymmMat(const int width, const int height)
    {
        w = width;
        h = height;
        size = width*height*3;
        pmat = new float[size];
        memset(pmat, 0, size*sizeof(float));
    }

    void pSet(const int x, const int y, const float *mat)
    {
        assert(mat);
        if(x >= w || y >= h)
            return ;

        int cusor = 3*(y*w + x);
        pmat[cusor] = mat[0] ;
        pmat[cusor+1] = mat[1] ;
        pmat[cusor+2] = mat[2] ;
    }

    void pSet(const int x, const int y, const float a, const float b, const float c)
    {
        if(x >= w || y >= h)
            return ;

        int cusor = 3*(y*w + x);
        pmat[cusor] = a;
        pmat[cusor+1] = b;
        pmat[cusor+2] = c;
    }

    void pGet(const int x, const int y, float &a, float &b, float &c)
    {
        if(x >= w || y >= h)
        {
            a = c = 0;
            b = 0;
            return ;
        }

        int cusor = 3*(y*w + x);
        assert(cusor < size);
        a = pmat[cusor];
        b = pmat[cusor+1];
        c = pmat[cusor+2];
    }

    ~ImageSymmMat()
    {
        delete [] pmat;
    }
};

#endif
