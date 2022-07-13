#ifndef NONDETECTOR_H
#define NONDETECTOR_H

#include "abstractdetector.h"

class NonDetector:public AbstractDetector
{

public:
    NonDetector();
public:
    bool paramsCheck();
    bool KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn);
    void writeKeypoint(const char*fn);
    ~NonDetector(){};

protected:
    bool keypDetect(const char *dstfn);
public:
    static void test();
    static void convert(const char *srcfn, const char* dstfn);
    static void convert2vgg(const char *srcfn, const char* dstfn);

};

#endif
