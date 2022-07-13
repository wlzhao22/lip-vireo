#ifndef THINNER_H
#define THINNER_H

#include "image.h"

/**********************************************************************************
Y.Y. Zhang, P.S.P. Wang, "A Parallel Thinning Algorithm with Two-Subiteration that
Generates One-Pixel-Wide Skeletons", ICPR'96


Implemented by Wan-Lei Zhao

@author Wan-Lei Zhao
@date   25-11-2011

***********************************************************************************/

class Thinner
{
    private:
        unsigned char *imgpix;
        unsigned int *map_WIN1;
        unsigned int powerk[8], width, height;
        static const unsigned int Pta2t1[256];
        static const unsigned int Pta2t2[256];
        bool _INVERT_MAP;

    protected:

        unsigned int NormSign(const float val);
        unsigned int IVNormSign(const float val);

        unsigned int (Thinner::*sign_func)(const float val);

        unsigned int WIN1(const unsigned int x, const unsigned int y);
        unsigned int WIN2(const unsigned int x, const unsigned int y);

        void WinIdx(const unsigned int x, const unsigned int y, unsigned int *p);

        unsigned int NN(const unsigned int x, const unsigned int y);
        bool         calc_WINs();

    public:
        Thinner(const bool IS_INVERT_MAP);
        bool         PTA2T(unsigned char *pix, const unsigned int width, const unsigned int height);
        bool         PTA2T(const Image *Img);

        virtual ~Thinner(){};
        static void test();
};

#endif
