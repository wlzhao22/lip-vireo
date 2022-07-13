#ifndef VIEWBOARD_H
#define VIEWBOARD_H

#include "image.h"

class ViewBoard
{
    public:
        ViewBoard(const int kpnum, const int patchSize);
        bool     addPatch(const float *patch, const int c0, const int r0);
        int      saveView(const char *dstfn);
        virtual ~ViewBoard();

    private:
        static const int rmargin, cmargin, col0;
        int w0, h0;
        int row, col, sx, sy;

        Image *img_view;
};

#endif
