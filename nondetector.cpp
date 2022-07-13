#include "nondetector.h"
#include "keypoint.h"
#include "iotool.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

NonDetector::NonDetector()
{
    cout<<"Detector ................................ non\n";
    this->DETECTOR = non;
}

bool NonDetector::paramsCheck()
{
    return true;
}

bool NonDetector::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    bool success = this->keypDetect(dstfn);
    if(!success)
    {
        cout<<"Loading keypoint from '"<<dstfn<<"' failed!\n";
        NonDetector::releaseKpList(this->kps);
        return false;
    }
    this->crntimg = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    this->leveli_kps = kps;

    if(strcmp(descfn, "")&&this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, 1.0);
    }

    delete this->crntimg;
    return true;
}

bool NonDetector::keypDetect(const char *kpfn)
{
    bool success = IOTool::getKeypoint(kpfn, this->kps);
    return success;
}

void NonDetector::convert(const char *srcfn, const char* dstfn)
{
    ofstream outStrm(dstfn);
    ifstream inStrm(srcfn);
    assert(outStrm.is_open());
    assert(inStrm.is_open());
    int dim, i;
    float numb, a, b, c, laplace, x, y;
    inStrm>>dim;
    dim = dim - 1;
    float *feat = new float[dim];
    inStrm>>numb;
    outStrm<<(int)round(numb)<<" 4 "<<dim<<endl;
    while(!inStrm.eof())
    {
        inStrm>>x;
        inStrm>>y;
        inStrm>>a;
        inStrm>>b;
        inStrm>>c;
        inStrm>>laplace;
        //cout<<x<<" "<<y<<" "<<a<<" "<<b<<" "<<c<<" "<<laplace<<endl;
        for(i = 0; i < dim; i++)
        {
            inStrm>>feat[i];
            //cout<<feat[i]<<" ";
        }
        //cout<<endl;
        //exit(0);
        outStrm<<(int)round(x)<<" "<<(int)round(y)<<endl;
        for(i = 1; i <= dim; i++)
        {
            if(i%12)
            {
                outStrm<<(feat[i-1])<<" ";
            }
            else
            {
                outStrm<<(feat[i-1])<<endl;
            }
        }
        if(dim%12)
        {
            outStrm<<endl;
        }

    }
    inStrm.close();
    outStrm.close();

}


void NonDetector::convert2vgg(const char *srcfn, const char* dstfn)
{
    ofstream outStrm(dstfn);
    ifstream inStrm(srcfn);
    assert(outStrm.is_open());
    assert(inStrm.is_open());
    int dim, i;
    float numb, a, b, c, x, y, lap;
    inStrm>>dim;
    dim = dim - 1;
    float *feat = new float[dim];
    inStrm>>numb;
    outStrm<<dim<<endl<<(int)round(numb)<<endl;
    while(!inStrm.eof())
    {
        inStrm>>x;
        inStrm>>y;
        inStrm>>a;
        inStrm>>b;
        inStrm>>c;
        inStrm>>lap;
        for(i = 0; i < dim; i++)
        {
            inStrm>>feat[i];
        }
        outStrm<<x<<" "<<y<<" ";
        outStrm<<a<<" "<<b<<" "<<c<<" ";
        for(i = 0; i < dim-1; i++)
        {
            outStrm<<feat[i]<<" ";
        }
        outStrm<<feat[i]<<endl;
    }
    inStrm.close();
    outStrm.close();

}


void NonDetector::test()
{
    const char *keypf0 = "/home/wallace/src/cpp/lip-vireo/data/desc/img_659";
    const char *desc0  = "/home/wallace/src/cpp/lip-vireo/data/desc/img_659";
    /**
    AbstractDetector *mydetector = new NonDetector();
    mydetector->Init(config, "SIFT");
    mydetector->KeypointBuild(pgmf1,keypf1,desc1);
    **/
    char keyp[1204], desc1[1024];
    for(int i = 2; i <= 7; i++)
    {
        sprintf(keyp, "%s%d.desc", keypf0, i);
        sprintf(desc1, "%s%d.pkeys", desc0, i);
        //NonDetector::convert(keyp, desc1);
        NonDetector::convert2vgg(keyp, desc1);
    }
}
