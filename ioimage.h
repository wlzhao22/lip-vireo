/**
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
**/

#ifndef IOIMAGE_H
#define IOIMAGE_H

#include <iostream>
#include <fstream>

extern "C"
{
    #include "./include/jpeglib.h"
    #include <png.h>
}

using namespace std;

class IOImage
{
   static const int PNG_BYTES_TO_CHECK = 8;
public:
    IOImage();
    static char *read_pgm(const char *fn, int &width, int &height,int *MaxVal);
    static int  write_pgm(const char *fn, const int w, const int h, const float* data, const unsigned char maxval,
                            const char* comment_string, const int channel);
    static int  write_pgm(const char *fn, const int w, const int h, const unsigned char*data, const unsigned char maxval,
                            const char* comment_string, const int channel);

    static char *read_bmp(const char *fn, int &width, int &height,int *maxVal,const int channel);
    static void write_bmp(const char *fn, const int w,const int h, const float*data,const int channel);
    static void write_bmp(const char *fn, const int w, const int h, const unsigned char *body,const int channel);
    static char *decmp_bmp(char *pbuffer, const int w, const int h, int *maxVal);

    static char *read_jpg(const char *srcfn,int &w, int &h, int *MaxVal, const int channel);
    static void write_jpg(const char *srcfn, const unsigned char *data, const int w,
                          const int h,const int ch, const int quality);
    static void write_jpg(const char *srcfn, const float *data,const int w,
                          const int h, const int ch, const int quality);

    static char *read_ppm(const char *srcfn, int &w, int &h, int *MaxVal, const int channel);
    static void write_ppm(const char *srcfn, const unsigned char *data,const int w, const int h, const int ch);
    static void write_ppm(const char *srcfn, const float *data,const int w, const int h, const int ch);

    static char *read_png(const char *srcfn, int &w, int &h, int *maxVal, const int channel);
    static void write_png(const char* srcfn, const int width, const int height, const char *data, const int ch);
    static void write_png(const char* srcfn, const int width, const int height, const float *data, const int ch);
    static void  test();
    static void testbmp();
    virtual ~IOImage();
};
#endif

