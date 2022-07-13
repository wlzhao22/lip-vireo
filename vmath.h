/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wan-Lei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**/

#ifndef VMATH_H
#define VMATH_H

#include <vector>
using namespace std;

class VMath {

public:
    static void  l2norm(float *vect,  const int dim);
    static int   l2norm(float *vectr, const int dim, const int nseg);
    static void  l1norm(float *vect,  const int dim);
    static void  l1norm(float *vect,  const int dim, const float factor);

    static bool  SIFTNorm(float *vect,     const int dim);
    static bool  sqrtSIFTNorm(float *vect, const int dim);

    static bool  powerLaw(float *vectr, const unsigned int d0, const float p0);
    static float dst_cos(const float vect1[], const float vect2[], const unsigned int d0);

    static void  smoothHist(float *vect,      const int n);
    static void  SmoothHistogram(float *vect, const int n);
    static void  smoothonePass(float *vect,   const int n);
    static void  smooth(float vect[], const unsigned int idxS, const unsigned int idxE);

    static float lgx(const float a, const float b);
    static float absx(const float val);

    static void  mInv33(float a[3][3], float b[3][3]);
    static float Det(const float *mat, const unsigned int row);

    static bool  invSymMat(const float a,      const float b, const float c, float mat[2][2]);
    static bool  sqrtSymMat(const float a, const float b, const float c, float U[2][2], float eigs[2]);
    static bool  eigvlSymMat(const float a,    const float b, const float c, float &eigen1, float &eign2);
    static bool  eigvtSymMat(const float a,    const float b, const float c, float eignv[2][2], float eigs[2]);
    static bool  inv_sqrtSymMat(const float a, const float b, const float c, float mat[2][2], const float thr0);
    static void  inv_sqrtSymMat(const float a, const float b, const float c, float mat[2][2], float ei[2], bool _norm_);

    static bool  normMat(const float m[2][2], float nmat[2][2]);
    static bool  normMat(const float m[2][2], float nmat[2][2], const float thr0);
    static bool  normMat(const float m[2][2], float nmat[2][2],  float &e1, float &e2, const float  thr0);
    static void  normMat(vector<vector<float> > & mat);

    static float *project(float *vect, const int dim, const float *prjMat,
                          const int row, const int col, const int dstdim);

    static void  printMat(const float mat[2][2]);
    static void  printVect(const float *vect, const int dim);
    static float maxVec(const float *v, const int len);
    static int   mod(const int x, const int b);
    static int   Sign(const float val);
    static float getMaxParabola(float x1, float y1, float x2, float y2,float x3, float y3);

    static void  dec2hex(int number, int* hex_array);
    static int   min(const int x, const int y);
    static int   max(const int x, const int y);

    /**modified Bessel function of the first kind**/
    static double besseli(const unsigned int nu, const float x);
    static double factorial(const unsigned int n);


    static void test();
};
#endif
