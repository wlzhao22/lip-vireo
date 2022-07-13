#include "intimage.h"
#include "vmath.h"

#include <cstdlib>

using namespace std;

void IntImage::setPixel(const int x, const int y, const double val)
{
    if(x >= 0 && y >= 0)
    {
        int loc = y*this->width + x;
        this->pix[loc] = val;
    }
    return ;
}


double IntImage::getPixel(const int x, const int y)
{
    if(x >= 0 && y >= 0)
    {
        int loc = y*this->width + x;
        return this->pix[loc];
    }
    return 0;
}


float IntImage::boxIntegral(const int x0, const int y0, const int cols, const int rows)
{
    int step = this->width;

    int r1 = y0 > this->height?(this->height-1):(y0 - 1);
    int c1 = x0 > this->width ?(this->width- 1):(x0 - 1);
    int r2 = (y0 + rows) > this->height?(this->height-1):(y0 + rows- 1);
    int c2 = (x0 + cols) > this->width ?(this->width- 1):(x0 + cols- 1);

    float A = 0.0f;
    float B = 0.0f;
    float C = 0.0f;
    float D = 0.0f;

    if (r1 >= 0 && c1 >= 0) A = this->pix[r1 * step + c1];
    if (r1 >= 0 && c2 >= 0) B = this->pix[r1 * step + c2];
    if (r2 >= 0 && c1 >= 0) C = this->pix[r2 * step + c1];
    if (r2 >= 0 && c2 >= 0) D = this->pix[r2 * step + c2];

    return VMath::max(0.0f, D - B - C + A);
}
