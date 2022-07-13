#ifndef DESCLJET_H
#define DESCLJET_H

#include "abstractdescriptor.h"
#include "filter.h"

#define LjSize 31

/************************************************************************
*@author Wanlei Zhao
*@version 1.0
*All rights reserved by Wanlei Zhao
*
*Anyone receive this code should not redistribute it to other people
*without permission of the author
*
*This code should only be used for non-commercial purpose
************************************************************************/

class DescLJet:public AbstractDescriptor
{

private:
    static const float MyPI2 ;
    static const float MyPI ;
    static const float lj_sigma;
    static const float factor;
    static const float factor1;
    static const int LjFeatLen;
private:
    float *lj_win;
    float *localjet;

    //kernels for all 14 jets

    float kernel_dx[LjSize];
    //3:dxx
    float kernel_dxx[LjSize][LjSize];
    //4:dxy
    float kernel_dxy[LjSize][LjSize];
    //5:dy2
    float kernel_dyy[LjSize][LjSize];
    //6:dxxx
    float kernel_dxxx[LjSize][LjSize];
    //7:dxxy
    float kernel_dxxy[LjSize][LjSize];
    //8:dxyy
    float kernel_dxyy[LjSize][LjSize];
    //9:dyyy
    float kernel_dyyy[LjSize][LjSize];
    //10:dxxxx
    float kernel_dxxxx[LjSize][LjSize];
    //11:dxxxy
    float kernel_dxxxy[LjSize][LjSize];
    //12:dxxyy
    float kernel_dxxyy[LjSize][LjSize];
    //13:dxyyy
    float kernel_dxyyy[LjSize][LjSize];
    //14:dyyyy
    float kernel_dyyyy[LjSize][LjSize];

public:
    DescLJet(DESC desc);
    ~DescLJet();

    int buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate);
    int buildPatchView(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate);

    int getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size);

    void getLocalJet(const float *patch);

    float inv_sigma(const float sigma,const int n);
    float Conv2D(float kernel[][LjSize],const float *patch);
    float Conv1x(float kernel[],const float *patch);
    float Conv1y(float kernel[],const float *patch);

    //static void cleanMat(float mat);

    void GaussianDx(const float sigma, float kernel[]);
    void GaussianDy(const float sigma,float kernel[]);
    void GaussianDx2_2D(const float sigma,float mat[][LjSize]);
    void GaussianDy2_2D(const float sigma,float mat[][LjSize]);
    void GaussianDxy_2D(const float sigma,float mat[][LjSize]);
    void GaussianDx3_2D(const float sigma,float mat[][LjSize]);
    void GaussianDx2y_2D(const float sigma,float mat[][LjSize]);
    void GaussianDxy2_2D(const float sigma,float mat[][LjSize]);
    void GaussianDy3_2D(const float sigma,float mat[][LjSize]);
    void GaussianDx4_2D(const float sigma,float mat[][LjSize]);
    void GaussianDx3y_2D(const float sigma,float mat[][LjSize]);
    void GaussianDx2y2_2D(const float sigma,float mat[][LjSize]);
    void GaussianDxy3_2D(const float sigma,float mat[][LjSize]);
    void GaussianDy4_2D(const float sigma,float mat[][LjSize]);
};
#endif
