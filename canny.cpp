#include "canny.h"
#include "image.h"

#include <cstring>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

const float Canny::STEPWIDTH               = 0.01;
const float Canny::PI0                     = 3.14159265;
const unsigned int Canny::EDGE_NUM_LIMIT   = 10;
const unsigned int Canny::MAX_FOLLOW_EDGES = 200;
const float Canny::BOOSTBLURFACTOR         = 90.0f;
const float Canny::NOISE                   = 100.0f;

void Canny::Gkernel(const float sigma, float *kernel, const int windowsize)
{
    assert(kernel);
    int i, center;
    float x, fx, sum = 0.0f;

    center = windowsize / 2;
    for(i=0; i < windowsize; i++)
    {
        x = (float)(i - center);
        fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
        kernel[i] = fx;
        sum += fx;
    }
    for(i = 0; i < windowsize; i++)
    {
        kernel[i] =  kernel[i]/sum;
    }
    return ;
}

float *Canny::GB1Dkernel(const float sigma, int & winSize)
{
    winSize = (int)floor(3.0*sigma + 0.5f);
    winSize = winSize < 1?1:winSize;
    winSize = 2*winSize + 1;

      int c = winSize/2, i = 0;
    float *kernel  = new float[winSize];
    float s2 = sigma*sigma, sum = 0, v, r;

    for(i = 0; i < winSize; i++)
    {
        r = i - c;
        v = (1.0 / (sqrt(2 * PI0) * sigma) )* exp(-(r*r) / (2 * s2));
        kernel[i] = v;
        sum += v;
    }

    for(i = 0; i < winSize; i++)
    {
        kernel[i] = kernel[i] / sum;
    }
    return kernel;
}

void Canny::radian_direct(const float *delta_x, const float *delta_y, const int rows,
                             const int cols, float **dir_radians, int xdirtag, int ydirtag)
{
    int r, c, pos;
    double dx, dy;
    assert(delta_x);
    assert(delta_y);

    float *dirim = new float[rows*cols];
    memset(dirim, 0, rows*cols*sizeof(float));

    *dir_radians = dirim;

    for(r = 0,pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++,pos++)
        {
            dx = (double)delta_x[pos];
            dy = (double)delta_y[pos];

            if(xdirtag == 1) dx = -dx;
            if(ydirtag == -1) dy = -dy;

            dirim[pos] = (float)angle_radians(dx, dy);
        }
    }
    delete [] dirim;
}

double Canny::angle_radians(const double x, const double y)
{
    double xu, yu, ang, retval;

    xu = fabs(x);
    yu = fabs(y);

    if((xu == 0))
    {
        retval=(0);
    }
    else
    {
        if ((yu == 0))
        {
            retval=(0);
        }
        else
        {
            ang = atan(yu/xu);
            if(x >= 0)
            {
                if(y >= 0) retval=(ang);
                else retval=(2*M_PI - ang);
            }
            else
            {
                if(y >= 0)
                    retval = (M_PI - ang);
                else
                    retval = (M_PI + ang);
            }
        }
    }

    return retval;
}

float Canny::mag_x_y(const float *delta_x, const float *delta_y, const int rows, const int cols,
                          float **magnitude)
{
    int r, c, pos;
    float sq1, sq2, max = 0;

    /** Allocate an image to store the magnitude of the gradient **/
    for(r = 0, pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++, pos++)
        {
            sq1 = delta_x[pos] * delta_x[pos];
            sq2 = delta_y[pos] * delta_y[pos];
            (*magnitude)[pos] = sqrt(sq1 + sq2);
            if((*magnitude)[pos] > max)
            {
                max = (*magnitude)[pos];
            }
        }
    }
    return max;
}


void Canny::get_dxdy(const float *blurimg, const int rows0, const int cols0, float **delta_x, float **delta_y)
{
    int r, c, pos;
    /** Compute the x-derivative. Adjust the derivative at the borders to avoid losing pixels **/

    for(r = 0; r < rows0; r++)
    {
        pos = r * cols0;
        (*delta_x)[pos] = blurimg[pos+1] - blurimg[pos];
        pos++;
        for(c = 1; c < (cols0-1); c++, pos++)
        {
            (*delta_x)[pos] = blurimg[pos+1] - blurimg[pos-1];
        }

        (*delta_x)[pos] = blurimg[pos] - blurimg[pos-1];
    }

    /** Compute the y-derivative. Adjust the derivative at the borders to avoid losing pixels **/

    for(c = 0; c < cols0; c++)
    {
        pos = c;
        (*delta_y)[pos] = blurimg[pos+cols0] - blurimg[pos];
        pos += cols0;
        for(r = 1; r < (rows0-1); r++, pos += cols0)
        {
            (*delta_y)[pos] = blurimg[pos+cols0] - blurimg[pos-cols0];
        }
        (*delta_y)[pos] = blurimg[pos] - blurimg[pos-cols0];
    }
}


void Canny::GBlur(const float *image, const int rows, const int cols, const float sigma, float **blurimg)
{
    int   r, c, rr, cc;
    int   dim, center;
    float *tempim, dot, sum;
    float *kernel = NULL;

    /*** Create a 1-dimensional gaussian smoothing kernel **/
    kernel = GB1Dkernel(sigma, dim);
    center = dim / 2;

    float *dirim = new float[rows*cols];
    memset(dirim, 0, rows*cols*sizeof(float));

    /*** Allocate a temporary buffer image and the smoothed image**/
    tempim = dirim;

    /*** Blur in the x - direction **/
    for(r = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++)
        {
            dot = sum = 0.0;
            for(cc = (-center); cc <= center; cc++)
            {
                if(((c+cc) >= 0) && ((c+cc) < cols))
                {
                    dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
                    sum += kernel[center+cc];
                }
            }
            tempim[r*cols+c] = dot/sum;
        }
    }
    /*** Blur in the y - direction.**/
    for(c = 0; c < cols; c++)
    {
        for(r = 0; r < rows; r++)
        {
            sum = 0.0;   dot = 0.0;
            for(rr=(-center); rr<=center; rr++)
            {
                if(((r+rr) >= 0) && ((r+rr) < rows))
                {
                    dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                    sum += kernel[center+rr];
                }
            }
            (*blurimg)[r*cols+c] = dot*BOOSTBLURFACTOR/sum;
        }
    }
    delete [] dirim;
    delete [] kernel;
    return ;
}


void Canny::non_max_supp(const float *mag, const float *gradx, const float *grady,
                         const int nrows, const int ncols, unsigned char *result)
{
    int   rowcount, colcount;
    const float *magrowptr, *magptr;
    const float *gxrowptr, *gxptr;
    const float *gyrowptr, *gyptr;
    float z1, z2;
    float m00, gx = 0, gy = 0;
    float mag1, mag2, xperp = 0, yperp = 0;
    unsigned char *resultrowptr, *resultptr;

    assert(mag);
    assert(gradx);
    assert(grady);
    assert(result);

    /*** Zero the edges of the result image **/
    memset(result, 0, sizeof(unsigned char)*nrows*ncols);
    /*** Suppress non-maximum points **/
    magrowptr    = mag + ncols + 1;
    gxrowptr     = gradx + ncols + 1;
    gyrowptr     = grady + ncols + 1;
    resultrowptr = result + ncols + 1;

    for(rowcount = 1; rowcount < (nrows-2);  rowcount++)
    {
        magptr    = magrowptr;
        gxptr     = gxrowptr;
        gyptr     = gyrowptr;
        resultptr = resultrowptr;
        for(colcount = 1; colcount < (ncols-2);  colcount++)
        {
            m00 = *magptr;
            if(m00 == 0)
            {
                *resultptr = (unsigned char) NOEDGE;
                xperp  = yperp = 0.0f;
            }else
            {
                xperp = -(gx = *gxptr)/((float)m00);
                yperp = (gy = *gyptr)/((float)m00);
            }

            if(gx >= 0)
            {
                if(gy >= 0)
                {
                    if (gx >= gy)
                    {
                        /** 111, Left point **/
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                        /**111, Right point **/
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /**110, Left point **/
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /**110, Right point **/
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if(gx >= -gy)
                    {
                        /**101, Left point **/
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                        /**101, Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }else
                    {
                        /**100, Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /**100, Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                    }
                }
            }else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {
                        /**011, Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /**011, Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }else
                    {
                        /**010, Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /**010, Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }else
                {
                    if (-gx > -gy)
                    {
                        /**001, Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /**001, Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /**000, Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /**000, Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            }

            /** Now determine if the current point is a maximum point **/
            if ((mag1 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }else
            {

                if ((mag2 > 0.0))
                {
                    *resultptr = (unsigned char) NOEDGE;
                }
                else
                {
                    if (mag2 == 0.0)
                        *resultptr = (unsigned char) NOEDGE;
                    else
                        {
                            *resultptr = (unsigned char) POSSIBLE_EDGE;
                        }
                }
            }
            /**visit next column**/
            magptr++;
            gxptr++;
            gyptr++;
            resultptr++;
        } /**inner-for loop */

        /**visit next row**/
        magrowptr   += ncols;
        gyrowptr    += ncols;
        gxrowptr    += ncols;
        resultrowptr+=ncols;
    }/**outer-for loop */
    return ;
}



void Canny::non_max_supp(const float *mag, const float *gradx, const float *grady,
                         const int nrows, const int ncols, unsigned char *result, CannyDg dg_ornt)
{
    int   rowcount, colcount;
    const float *magrowptr, *magptr;
    const float *gxrowptr, *gxptr;
    const float *gyrowptr, *gyptr;
    float z1, z2;
    float m00, gx = 0, gy = 0;
    float mag1, mag2, xperp = 0.0f, yperp = 0.0f;
    unsigned char *resultrowptr, *resultptr;

    assert(mag);
    assert(gradx);
    assert(grady);
    assert(result);

    /*** Zero the edges of the result image **/
    memset(result, 0, sizeof(unsigned char)*nrows*ncols);
    /*** Suppress non-maximum points **/
    magrowptr    = mag + ncols + 1;
    gxrowptr     = gradx + ncols + 1;
    gyrowptr     = grady + ncols + 1;
    resultrowptr = result + ncols + 1;

    for(rowcount = 1; rowcount < (nrows-2);  rowcount++)
    {
        magptr    = magrowptr;
        gxptr     = gxrowptr;
        gyptr     = gyrowptr;
        resultptr = resultrowptr;
        for(colcount = 1; colcount < (ncols-2);  colcount++)
        {
            m00 = *magptr;
            if(m00 == 0)
            {
                *resultptr = (unsigned char) NOEDGE;
                xperp = yperp = 0.0f;
            }else
            {
                xperp = -(gx = *gxptr)/((float)m00);
                yperp = (gy = *gyptr)/((float)m00);
            }

            if(gx >= 0 && dg_ornt == _pdx)
            {
                if(gy >= 0/** && dg_ornt == _pdy**/)
                {
                    if (gx >= gy)
                    {
                        /** 111, Left point **/
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                        /**111, Right point **/
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /**110, Left point **/
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /**110, Right point **/
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                    }
                }
                else if(gy < 0 /** && dg_ornt == _ndy**/)
                {
                    if(gx >= -gy)
                    {
                        /**101, Left point **/
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                        /**101, Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }else
                    {
                        /**100, Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /**100, Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                    }
                }
            }else if(gx < 0 && dg_ornt == _ndx)
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {
                        /**011, Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /**011, Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }else
                    {
                        /**010, Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /**010, Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }else
                {
                    if (-gx > -gy)
                    {
                        /**001, Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /**001, Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /**000, Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /**000, Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            }else{
                mag1 = mag2 = 0.0f;
            }

            /** Now determine if the current point is a maximum point **/
            if ((mag1 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }else
            {

                if ((mag2 > 0.0))
                {
                    *resultptr = (unsigned char) NOEDGE;
                }
                else
                {
                    if (mag2 == 0.0)
                        *resultptr = (unsigned char) NOEDGE;
                    else
                        {
                            *resultptr = (unsigned char) POSSIBLE_EDGE;
                        }
                }
            }
            /**visit next column**/
            magptr++;
            gxptr++;
            gyptr++;
            resultptr++;
        } /**inner-for loop */

        /**visit next row**/
        magrowptr   += ncols;
        gyrowptr    += ncols;
        gxrowptr    += ncols;
        resultrowptr+=ncols;
    }/**outer-for loop */
    return ;
}

void Canny::track_edges(unsigned char *edgemapptr, const float *edgemagptr, const short lowval, const int cols)
{
    const float *centermagptr, *tempmagptr;
    unsigned char *centermapptr;
    unsigned char *tempmapptr;
    unsigned int array_ctr = 0;
    unsigned char* mapptr[MAX_FOLLOW_EDGES+1];
    const float* magptr[MAX_FOLLOW_EDGES+1];
    unsigned int i, follow_ctr = 0;
    int x[8] = {1,1,0,-1,-1,-1,0,1};
    int y[8] = {0,1,1,1,0,-1,-1,-1};

    centermapptr = edgemapptr;
    centermagptr = edgemagptr;
    do  /* MAXITER = MAX_FOLLOW_EDGES */
    {
        for(i = 0; i < 8; i++) /* MAXITER = 8 */
        {
            tempmapptr = centermapptr - y[i]*cols + x[i];
            tempmagptr = centermagptr - y[i]*cols + x[i];
            if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
            {
                *tempmapptr = (unsigned char) EDGE;
                if (array_ctr < MAX_FOLLOW_EDGES)
                {
                    mapptr[array_ctr]= tempmapptr;
                    magptr[array_ctr]= tempmagptr;
                    array_ctr++;
                }
            }
        }
        centermapptr = mapptr[follow_ctr];
        centermagptr = magptr[follow_ctr];
        follow_ctr++;
    }while(follow_ctr <= array_ctr);
}


/****************************************************************************
 * Initialize the edge map to possible edges everywhere the non-maximal
 * suppression suggested, there could be an edge except for the border. At
 * the border we say there can not be an edge because it makes the
 * follow_edges algorithm more efficient to not worry about tracking an
 * edge off the side of the image.
 ****************************************************************************/

void Canny::apply_hysteresis(const float *mag, const float max_mag, unsigned char *nms, const int rows, const int cols,
                             const float tlow, const float thigh, unsigned char *edge)
{
    assert(max_mag > 0);

    int c, r, pos, numedges, highcount, lowthreshold, highthreshold;
    int num_bin = (unsigned int)floor(max_mag);
    int *hist = new int[num_bin+1];
    int edge_ctr = 0;
    int maximum_mag = 0;

    memset(hist, 0, sizeof(int)*num_bin);

    for(r = 0, pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++,pos++)
        {
            if(nms[pos] == POSSIBLE_EDGE)
            {
                edge[pos] = POSSIBLE_EDGE;
            }else{
                edge[pos] = NOEDGE;
            }
        }
    }

    for(r = 0, pos = 0; r < rows; r++,pos += cols)
    {
        edge[pos] = NOEDGE;
        edge[pos+cols-1] = NOEDGE;
    }

    pos = (rows-1) * cols;

    for(c = 0; c < cols; c++,pos++)
    {
        edge[c]   = NOEDGE;
        edge[pos] = NOEDGE;
    }

    /// Compute the histogram of the magnitude image. use histogram to compute hysteresis thresholds

    for(r = 0; r < num_bin; r++)
    {
        hist[r] = 0;
    }
    unsigned int idx = 0;
    for(r = 0, pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++,pos++)
        {
            if(edge[pos] == POSSIBLE_EDGE)
            {
                idx = (unsigned int)floor(mag[pos]); //wanlei zhao, potential bugs;
                hist[idx]++;
            }
        }
    }

    /// Compute the number of pixels that passed the nonmaximal suppression

    for(r = 1, numedges = 0; r < num_bin; r++)
    {
        if(hist[r] != 0)
        {
            maximum_mag = r;
        }
        numedges += hist[r];
    }

    highcount = (int)floor(numedges * thigh + 0.5);
    r = 1;
    numedges = hist[1];
    {
        int tmp = 0;
        while((tmp = (r < (maximum_mag-1)) && (numedges < highcount)))
        {
            r++;
            numedges += hist[r];
        }
    }

    highthreshold = r;
    lowthreshold = (int)floor(highthreshold * tlow + 0.5);

    /// looks for pixels above the highthreshold to locate edges, calls follow_edges to continue the edge

    for(r = 0, pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++, pos++)
        {
            if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)) // && edge_ctr < EDGE_NUM_LIMIT
            {
                edge[pos] = EDGE;
                edge_ctr++;
                track_edges((edge+pos), (mag+pos), lowthreshold, cols);
            }
        }
    }

    /// Set all the remaining possible edges to non-edges

    for(r = 0, pos = 0; r < rows; r++)
    {
        for(c = 0; c < cols; c++,pos++)
        {
            if(edge[pos] != EDGE)
            {
                edge[pos] = NOEDGE;
            }
        }
    }

    delete [] hist;
    return ;
}

unsigned char *Canny::canny(const float *image, const int rows, const int cols, const float sigma,
                  const float tlow, const float thigh)
{
    unsigned char *nms;        /** Points that are local maximal magnitude **/
    float *magnitude;
    float *dir_radians;
    float max_mag;

    unsigned char *uctemp = new unsigned char[rows*cols];
    unsigned char *edge   = new unsigned char [rows*cols];

    float *blurimg        = new float[rows*cols];
    float *delta_x        = new float[rows*cols];
    float *delta_y        = new float[rows*cols];

    memset(uctemp,  0, rows*cols*sizeof(unsigned char));
    memset(blurimg, 0, rows*cols*sizeof(float));
    memset(delta_x, 0, rows*cols*sizeof(float));
    memset(delta_y, 0, rows*cols*sizeof(float));
    memset(edge,    0, rows*cols*sizeof(unsigned char));

    /** Perform gaussian smoothing on the image using the input standard deviation **/
    GBlur(image, rows, cols, sigma, &blurimg);
    get_dxdy(blurimg, rows, cols, &delta_x, &delta_y);

    /** Compute the direction up the gradient, in radians that are
    * specified counteclockwise from the positive x-axis.   **/

    radian_direct(delta_x, delta_y, rows, cols, &dir_radians, -1, -1);
    magnitude = blurimg;
    max_mag = mag_x_y(delta_x, delta_y, rows, cols, &magnitude);

    /** Perform non-maximal suppression. **/
    nms = uctemp;
    non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

    /** Use hysteresis to mark the edge pixels.**/
    apply_hysteresis(magnitude, max_mag, nms, rows, cols, tlow, thigh, edge);

    delete [] uctemp;
    delete [] delta_x;
    delete [] delta_y;
    delete [] blurimg;

    return edge;
}

void Canny::test()
{
    int count, i, j;
    float sigma, tlow, thigh;

    sigma = 6.8;    tlow  = 0.15;  thigh = 0.75;

    const char *srcimg   = "/home/wlzhao/datasets/vgg/graf/graf1.jpg";
    const char *dstimg1   = "/home/wlzhao/lectures/mir/lec6/pics/graf1_edge3.jpg";

    Image *src  = new Image(srcimg);
    Image *dst  = new Image(src->width, src->height);
    unsigned char *tmpdst;

    tmpdst = Canny::canny(src->pix, src->height, src->width, sigma, tlow, thigh);
    count = 0;
    for(i = 0; i < src->height; i++)
    {
        for(j = 0; j < src->width; j++)
        {
            dst->pix[count] = tmpdst[count];
            count++;
        }
    }

    dst->save(dstimg1);

    delete [] tmpdst;
    return ;
}
