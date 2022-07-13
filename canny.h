#ifndef CANNY_H
#define CANNY_H

/**********************************************************************************************
This code is collected from Mike Heath, who was once with University of South Florida. The code is
 modified by Wan-lei Zhao on Nov. 2011. Potential bugs have been removed in this updated version.

Website for original code: http://marathon.csee.usf.edu/edge/edge_detection.html

Step:
*   1) Convolve the image with a separable gaussian filter.
*   2) Take the dx and dy the first derivatives using [-1,0,1] and [1,0,-1].
*   3) Compute the magnitude: sqrt(dx*dx+dy*dy).
*   4) Perform non-maximal suppression.
*   5) Perform hysteresis.

@author: Mike Heath, Wan-Lei Zhao
@date:   2/15/96; 11/7/2011
@parameters:
@sigma [>=1, 10]: The standard deviation of the gaussian smoothing filter
@tlow  [0, 1]:    Specifies the low value to use in hysteresis.
@thigh [0, 1]:    Specifies the high value to use in hysteresis. This fraction (0-1)
*                 specifies the percentage point in a histogram of the gradient of
*                 the magnitude. Magnitude values of zero are not counted in the
*                 histogram.

@tips: use as a starting point are: sigma 0.60-2.40, tlow 0.20-0.50 and thigh 0.60-0.90

***********************************************************************************************/

enum CannyDg {_pdx = 0, _ndx, _pdy, _ndy};

class Canny
{
    static const float STEPWIDTH;
    static const int   NOEDGE = 255;
    static const int   POSSIBLE_EDGE = 128;
    static const int   EDGE   = 0;
    static const float BOOSTBLURFACTOR;
    static const unsigned char NOTHING = 50;
    static const unsigned char CURVE   = 120;
    static const float NOISE;
    static const float PI0;

    static const unsigned int EDGE_NUM_LIMIT;
    static const unsigned int MAX_FOLLOW_EDGES;

protected:

    static void  GBlur(const float *image, const int rows, const int cols, const float sigma, float **blurimg);
    static float *GB1Dkernel(const float sigma, int & winSize);
    static void  Gkernel(const float sigma, float *kernel, const int windowsize);
    static void  get_dxdy(const float *smoothedim, const int rows, const int cols, float **delta_x, float **delta_y);
    static float mag_x_y(const float *delta_x, const float *delta_y, const int rows, const int cols, float **magnitude);
    static void  track_edges(unsigned char *edgemapptr, const float *edgemagptr, const short lowval, const int cols);

    static void  apply_hysteresis(const float *mag, const float max_mag, unsigned char *nms, const int rows, const int cols,
                                 const float tlow, const float thigh, unsigned char *edge);
    static void  radian_direct(const float *delta_x, const float *delta_y, const int rows, const int cols,
                                 float **dir_radians, int xdirtag, int ydirtag);

    static double angle_radians(const double x, const double y);
    static void   non_max_supp(const float *mag, const float *gradx, const float *grady,
                         const int nrows, const int ncols, unsigned char *result);
    static void   non_max_supp(const float *mag, const float *gradx, const float *grady,
                         const int nrows, const int ncols, unsigned char *result, CannyDg dg_ornt);
public:
    static unsigned char *canny(const float *image, const int rows, const int cols,
                                 const float sigma, const float tlow, const float thigh);
    Canny(){}
    virtual ~Canny(){}
    static void test();
};

#endif
