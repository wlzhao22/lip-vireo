#include "thinner.h"
#include <cstring>
#include <cstdlib>
#include <cmath>


using namespace std;

/**************************************************************
The following two tables are generated according to the paper.
They are kept as constant all the time.

The names of functions have been defined in accordance with the ones
in the paper.

The following 4 templates have been additionally considered, which have been
missed in original paper.

T01-1
------
0 1 *
0 p 1
0 0 0
5,

T01-4
0 1 0
0 p 1
0 1 0
21

T01-3
0 1 0
0 p 1
1 1 0
53

T01-2
* 1 *
1 p 0
0 0 0
65

***************************************************************/

///#define DEBUG_THINNER

const unsigned int Thinner::Pta2t1[] =
{
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,  0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
    1, 1, 1, 0, 1, 0, 0, 0, 1, 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 1, 1, 1, 0, 1, 0, 0, 0,
    1, 0, 1, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 3, 0, 4, 1, 0, 1, 0, 1, 3,  0, 4, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1, 3, 0, 4, 1, 0, 1, 0,
    1, 3, 0, 4, 1, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 2, 1, 0, 1, 0, 1, 1, 0, 2,
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 2, 1, 0, 1, 0, 1, 1,  0, 2, 1, 0, 1, 0
};

//keep it as original

const unsigned int Thinner::Pta2t2[] =
{
    0, 0, 0, 1, 0, 1, 0, 1, 0, 0,  0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
    0, 3, 0, 3, 0, 0, 0, 0, 0, 1,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 4, 0, 4, 0, 0, 0, 0,
    0, 2, 0, 2, 0, 1, 0, 1, 0, 0,  0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 1, 0, 1, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 1, 1, 0, 1,
    0, 1, 0, 1, 0, 0, 0, 0, 0, 1,  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 1, 0, 1, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 0, 0, 0, 0, 3,  0, 3, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 4, 0, 4, 0, 0, 0, 0, 0, 2,  0, 2, 1, 1, 0, 1, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 1, 0, 1,  0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0
};

Thinner::Thinner(const bool _INVERT_MAP0)
{
    powerk[0] = 1;
    for(unsigned int k = 1; k < 8; k++)
    {
        powerk[k] = powerk[k-1]*2;
    }
    _INVERT_MAP = _INVERT_MAP0;
    if(_INVERT_MAP0)
    {
        sign_func = &Thinner::IVNormSign;
    }
    else
    {
        sign_func = &Thinner::NormSign;
    }
}

unsigned int Thinner::WIN1(const unsigned int x, const unsigned int y)
{
    unsigned int p[8];
    unsigned int k, val = 0;
    unsigned int w = width;
    unsigned int h = height;
    unsigned int loc = w*y+x;

    if(x == 0 || y == 0)
    {
        return 0;
    }

    if(x >= (w-1) || y >= (h-1))
    {
        return 0;
    }

    if((this->*sign_func)(this->imgpix[loc]) == 0)
    {
        return 0;
    }

    for(k = 0; k < 8; k++)
    {
        p[k] = loc;
    }

    p[7] = p[7] - w - 1;
    p[0] = p[0] - w;
    p[1] = p[1] - w + 1;

    p[2] = p[2] + 1;
    p[3] = p[3] + w + 1;
    p[4] = p[4] + w;

    p[5] = p[5] + w - 1;
    p[6] = p[6] - 1;

    for(k = 0; k < 8; k++)
    {
        if((this->*sign_func)(this->imgpix[p[k]]) == 1)
            val += powerk[k];
    }
    return val;
}

unsigned int Thinner::WIN2(const unsigned int x, const unsigned int y)
{
    unsigned int p[8];
    unsigned int k, val = 0;
    unsigned int w = width;
    unsigned int h = height;
    unsigned int loc = w*y+x;

    if(x == 0 || y == 0)
    {
        return 0;
    }

    if(x >= (w-1) || y >= (h-1))
    {
        return 0;
    }

    if((this->*sign_func)(this->imgpix[loc]) == 0)
    {
        return 0;
    }

    for(k = 0; k < 8; k++)
    {
        p[k] = loc;
    }

    p[7] = p[7] + w + 1;
    p[0] = p[0] + w;
    p[1] = p[1] + w - 1;

    p[2] = p[2] - 1;
    p[3] = p[3] - w - 1;
    p[4] = p[4] - w;

    p[5] = p[5] - w + 1;
    p[6] = p[6] + 1;

    for(k = 0; k < 8; k++)
    {
        if((this->*sign_func)(this->imgpix[p[k]]) == 1)
            val += powerk[k];
    }

    return val;
}


void Thinner::WinIdx(const unsigned int x, const unsigned int y, unsigned int *p)
{
    unsigned int k;
    unsigned int w = width;
    unsigned int h = height;

    assert(p);
    for(k = 0; k < 8; k++)
    {
        p[k] = 0;
    }

    if(x == 0 || y == 0)
    {
        return ;
    }

    if(x >= (w-1) || y >= (h-1))
    {
        return ;
    }

    for(k = 0; k < 8; k++)
    {
        p[k] = w*y+x;
    }

    p[7] = p[7] - w - 1;
    p[0] = p[0] - w;
    p[1] = p[1] - w + 1;

    p[2] = p[2] + 1;
    p[3] = p[3] + w + 1;
    p[4] = p[4] + w;

    p[5] = p[5] + w - 1;
    p[6] = p[6] - 1;

    for(k = 0; k < 8; k++)
    {
        p[k] = map_WIN1[p[k]];
    }
    return ;
}


unsigned int Thinner::NormSign(const float val)
{
    if(val > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

unsigned int Thinner::IVNormSign(const float val)
{
    if(floor(val) == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

bool Thinner::calc_WINs()
{
    unsigned int w = width -1;
    unsigned int h = height -1;
    unsigned int i, j, loc;

    memset(map_WIN1, 0, sizeof(unsigned int)*(w+1)*(h+1));

    for(i = 1; i < h; i++)
    {
        loc =  i*(w+1);
        for(j = 1; j < w; j++)
        {
            map_WIN1[loc+j] = Thinner::WIN1(j, i);
            #ifdef DEBUG_THINNER
            if(srcimg->pix[loc+j] == 0)
                counter[0]++;
            else
            {
                counter[1]++;
            }
            #endif
        }
    }
    #ifdef DEBUG_THINNER
    cout<<"display: "<<counter[0]<<"\t"<<counter[1]<<endl;
    #endif

    return true;
}


unsigned int Thinner::NN(const unsigned int x, const unsigned int y)
{
    unsigned int p[8];
    unsigned int k, val = 0;
    unsigned int w = width;
    unsigned int h = height;

    if(x == 0 || y == 0)
    {
        return 0;
    }

    if(x >= (w-1) || y >= (h-1))
    {
        return 0;
    }

    for(k = 0; k < 8; k++)
    {
        p[k] = w*y+x;
    }
    p[7] = p[7] + w + 1;
    p[0] = p[0] + w;
    p[1] = p[1] + w - 1;

    p[2] = p[2] - 1;
    p[3] = p[3] - w - 1;
    p[4] = p[4] - w;

    p[5] = p[5] - w + 1;
    p[6] = p[6] + 1;

    for(k = 0; k < 8; k++)
    {
        val += (this->*sign_func)(this->imgpix[p[k]]);
    }

    return val;
}

bool Thinner::PTA2T(const Image *Img)
{
    width                 = Img->width;
    height                = Img->height;
    unsigned char *tmppix = new unsigned char[width*height];
    unsigned int    x, y, loc;
    for(y = 0; y < height; y++)
    {
         loc = y*Img->width;
         for(x = 0; x < width; x++)
         {
             tmppix[loc+x] = (unsigned char)floor(Img->pix[loc+x]);
         }
    }

    this->PTA2T(tmppix, width, height);

    for(y = 0; y < height; y++)
    {
         loc = y*Img->width;
         for(x = 0; x < width; x++)
         {
             Img->pix[loc+x] = tmppix[loc+x];
         }
    }

    delete [] tmppix;

    return true;
}


bool Thinner::PTA2T(unsigned char *pix, const unsigned int width0, const unsigned int height0)
{
    assert(pix != NULL);
    this->imgpix = pix;
    width        = width0;
    height       = height0;

    bool THIN_NOT_END = true, STABLE = false;
    unsigned int   k, x, y, loc, idx, map_idx;
    unsigned int   win[8] = {0};
    unsigned int   w = width -  1;
    unsigned int   h = height - 1;
    unsigned int   *map_pp = NULL;
    const unsigned int *table_pp = NULL;

    map_WIN1 = new unsigned int[(w+1)*(h+1)];

    while(THIN_NOT_END)
    {
        THIN_NOT_END = false;
        map_pp   = map_WIN1;
        table_pp = Pta2t1;

        for(k = 1; k <= 2; k++)
        {
            this->calc_WINs();
            for(y = 1; y < h; y++)
            {
                loc = y*this->width;
                for(x = 1; x < w; x++)
                {
                    idx = loc + x;
                    if((this->*sign_func)(this->imgpix[idx]) == 0)
                    {
                        continue;
                    }

                    this->WinIdx(x, y, win);

                    STABLE = false;
                    map_idx = map_pp[idx];
                    switch(table_pp[map_idx])
                    {
                    case 0:
                    {
                        break;
                    }
                    case 1:
                    {
                        STABLE = true;
                        break;
                    }
                    case 2:
                    {
                        if(table_pp[win[0]] == 0)
                        {
                            STABLE = true;
                        }
                        break;
                    }
                    case 3:
                    {
                        if(table_pp[win[6]] == 0)
                        {
                            STABLE = true;
                        }
                        break;
                    }
                    case 4:
                    {
                        if((table_pp[win[0]] == 0) && (table_pp[win[6]] == 0))
                        {
                            STABLE = true;
                        }
                        break;
                    }
                    }//switch-of
                    if(STABLE)
                    {
                        this->imgpix[idx] = _INVERT_MAP?255:0;
                        THIN_NOT_END     = true;
                    }
                }
            }//[row, col]
        }//for-k-loop
    }//while-loop

    delete [] map_WIN1;
    this->imgpix = NULL;
    return true;
}

void Thinner::test()
{
    const char *srcfn1 = "/home/wlzhao/datasets/query1_edge.pgm";
    const char *dstfn  = "/home/wlzhao/datasets/query_edge.pgm";
    Image *srcimg = new Image(srcfn1);
    Thinner *mythin = new Thinner(true);

    mythin->PTA2T(srcimg);
    srcimg->save(dstfn);
    delete mythin;
    delete srcimg;
    return ;
}
