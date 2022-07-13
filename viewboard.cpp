#include "viewboard.h"
#include <cmath>


const int ViewBoard::rmargin = 3;
const int ViewBoard::cmargin = 3;
const int ViewBoard::col0    = 13;

ViewBoard::ViewBoard(const int kpnum, const int patchSize)
{
    this->row = ceil((float)kpnum/(float)col0);
    w0 = patchSize*col0 + (col0-1)*ViewBoard::cmargin + 2;
    h0 = patchSize*row  + (row-1)*ViewBoard::rmargin  + 2;

    img_view = new Image(w0, h0);

    sx = sy = 0;
    for(int i = 0; i < h0; i++)
    {
        for(int j = 0; j < w0; j++)
        {
            img_view->setPixel(j, i, 255);
        }
    }
}

bool ViewBoard::addPatch(const float *patch, const int c0, const int r0)
{
    assert(patch);
    int i, j, x, y, idx = 0;

    if((sx+c0) > this->img_view->width)
    {
        sy += (ViewBoard::rmargin + r0);
        sx  = 0;
    }

    if((sy+r0) >= this->img_view->height)
    {
        cout<<"Alert: view board is full! (in ViewBoard::addPatch)\n";
        return false;
    }

    for(i = 0; i < c0; i++)
    {
        y = i + sy;
        for(j = 0; j < r0; j++)
        {
            //copy patch
            x = j + sx;
            img_view->setPixel(x, y, patch[idx]);
            idx++;
        }
    }

    sx += (ViewBoard::cmargin + c0);
    if((sx+c0) > this->img_view->width)
    {
        sy += (ViewBoard::rmargin + r0);
        sx  = 0;
    }

    return true;
}

int ViewBoard::saveView(const char *dstfn)
{
    assert(dstfn);
    img_view->save(dstfn);
    return 0;
}

ViewBoard::~ViewBoard()
{
    delete img_view;
}
