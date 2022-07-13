#ifndef HAAR_H
#define HAAR_H
#include "image.h"


class Haar
{

public:
    static int upbound;
    static void haar1d(float *vec, int n);
    static void haar1(float *vec, const int n, const int w);
    static void haar2(float *matrix, const int rows, const int cols);
    static void haar2(Image *srcimg);
    static int test();
    static int test1();
};

#endif
