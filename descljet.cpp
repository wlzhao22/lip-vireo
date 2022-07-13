#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "viewboard.h"
#include "descljet.h"
#include "filter.h"
#include "vmath.h"

const float DescLJet::MyPI2    = 6.28318531f;
const float DescLJet::MyPI     = 3.14159265358979323846f;
const float DescLJet::lj_sigma = 5.0f;
const int LJFeatLen            = 15;
const float DescLJet::factor   = 1;//400000;
const float DescLJet::factor1  = 1;//2000;

DescLJet::DescLJet(DESC desc)
{
    this->lj_win = new float[LjSize*LjSize];
    this->localjet = new float[LJFeatLen];
    this->featLen = LJFeatLen;
    this->descoption = desc;

    switch(desc)
    {
    case LJET:
    {
        cout<<"Descriptor .............................. unnormalized LJET\n";
        this->featLen = LJFeatLen;
        break;
    }
    case NLJET:
    {
        cout<<"Descriptor .............................. normalized LJET\n";
        this->featLen = LJFeatLen;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = LJFeatLen;
    }
    }

    DescLJet::GaussianDx(DescLJet::lj_sigma,       kernel_dx);
    DescLJet::GaussianDx2_2D(DescLJet::lj_sigma,   kernel_dxx);
    DescLJet::GaussianDxy_2D(DescLJet::lj_sigma,   kernel_dxy);
    DescLJet::GaussianDy2_2D(DescLJet::lj_sigma,   kernel_dyy);
    DescLJet::GaussianDx3_2D(DescLJet::lj_sigma,   kernel_dxxx);
    DescLJet::GaussianDx2y_2D(DescLJet::lj_sigma,  kernel_dxxy);
    DescLJet::GaussianDxy2_2D(DescLJet::lj_sigma,  kernel_dxyy);
    DescLJet::GaussianDy3_2D(DescLJet::lj_sigma,   kernel_dyyy);
    DescLJet::GaussianDx4_2D(DescLJet::lj_sigma,   kernel_dxxxx);
    DescLJet::GaussianDx3y_2D(DescLJet::lj_sigma,  kernel_dxxxy);
    DescLJet::GaussianDx2y2_2D(DescLJet::lj_sigma, kernel_dxxyy);
    DescLJet::GaussianDxy3_2D(DescLJet::lj_sigma,  kernel_dxyyy);
    DescLJet::GaussianDy4_2D(DescLJet::lj_sigma,   kernel_dyyyy);
}

int DescLJet::buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    ofstream outStrm(descfn);

    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(0);
    }

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
         if(properties[i])
         nproperty += KP_PROP_SZ[i];
    }

    int kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }
        kpnumb++;
    }

    if(this->_out_format == _vgg_fmrt)
        outStrm<<this->featLen<<endl<<kpnumb<<endl;
    else if(this->_out_format == _vireo_fmrt)
        outStrm<<kpnumb<<" "<<nproperty<<" "<<this->featLen<<endl;

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }

        this->getNormDescPatch(crntKpt, this->lj_win, LjSize);
        this->getLocalJet(this->lj_win);

        if(this->descoption == NLJET)
        {
            VMath::l2norm(this->lj_win, this->featLen);
        }

        if(this->_out_format == _vgg_fmrt)
        {
            this->saveDescVGG(crntKpt, this->localjet, this->featLen, resize_rate, outStrm);
        }
        else if(this->_out_format == _vireo_fmrt)
        {
            this->saveDescVireo(crntKpt, this->localjet, kpnumb, this->featLen, resize_rate, outStrm);
        }
        kpnumb++;
    }

    outStrm.close();
    return kpnumb;
}

int DescLJet::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }
        kpnumb++;
    }

    ViewBoard *myview = new ViewBoard(kpnumb, LjSize);

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }
        this->getNormDescPatch(crntKpt, this->lj_win, LjSize);
        myview->addPatch(this->lj_win, LjSize, LjSize);
    }
    myview->saveView(dvfn);
    delete myview;
    return 0;
}

int DescLJet::getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sine = 0, cosine = 1;
    float rPtx = 0, rPty = 0;
    float gray = 0;


    int radius = (int)floor(Size/2.0);
    float sc = (2*keyp->iscale)/Size;
    int i, j, irow = 0;


    sine = sin(keyp->ori);
    cosine = cos(keyp->ori);

    sine   = sine*sc;
    cosine = cosine*sc;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            rPtx = j*cosine + i*sine;
            rPty = i*cosine - j*sine;

            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }
    return 1;
}


void DescLJet::getLocalJet(const float *patch)
{
    //0: pixel value
    int center = LjSize*(LjSize/2)+LjSize/2;
    this->localjet[0] = patch[center]*DescLJet::factor1;

    //1:dx
    this->localjet[1] = this->Conv1x(kernel_dx, patch);
    ///cout<<this->localjet[1]<<endl;
    //cout<<"dx\n";

    //2:dy
    this->localjet[2] = this->Conv1y(kernel_dx, patch);
    //cout<<"dy\n";

    //3:dxx
    this->localjet[3] = this->Conv2D(kernel_dxx,patch);
    //cout<<"dxx\n";

    //4:dxy
    this->localjet[4] = this->Conv2D(kernel_dxy,patch);
    //cout<<"dxy\n";

    //5:dy2
    this->localjet[5] = this->Conv2D(kernel_dyy,patch);
    //cout<<"dy2\n";

    //6:dxxx
    this->localjet[6] = this->Conv2D(kernel_dxxx,patch);
    //cout<<"dxxx\n";

    //7:dxxy
    this->localjet[7] = this->Conv2D(kernel_dxxy,patch);
    //cout<<"dxxy\n";

    //8:dxyy
    this->localjet[8] = this->Conv2D(kernel_dxyy,patch);
    //cout<<"dxyy\n";

    //9:dyyy
    this->localjet[9] = this->Conv2D(kernel_dyyy, patch);
    //cout<<"dyyy\n";

    //10:dxxxx
    this->localjet[10] = this->Conv2D(kernel_dxxxx, patch);
    //cout<<"dxxxx\n";

    //11:dxxxy
    this->localjet[11] = this->Conv2D(kernel_dxxxy, patch);
    //cout<<"dxxxy\n";

    //12:dxxyy
    this->localjet[12] = this->Conv2D(kernel_dxxyy,patch);
    //cout<<"dxxyy\n";

    //13:dxyyy
    this->localjet[13] = this->Conv2D(kernel_dxyyy,patch);
    //cout<<"dxyyy\n";

    //14:dyyyy
    this->localjet[14] = this->Conv2D(kernel_dyyyy, patch);
    //cout<<"dyyyy\n";


}

float DescLJet::Conv1x(float kernel[],const float *patch)
{
    int ix,loc;
    //assert(kernel.size() == DescLJet::LjSize);
    float sumVal = 0,pixel;
    int ic = LjSize/2;
    loc = ic*LjSize;

    for(ix = 0; ix < LjSize; ix++)
    {
        loc = loc+ ix;
        pixel = patch[loc];
        sumVal = sumVal+kernel[ix]*pixel;
    }

    return sumVal;
}

float DescLJet::Conv1y(float kernel[],const float *patch)
{
    int iy,loc = 0;
    //assert(kernel.size() == DescLJet::LjSize);
    float sumVal = 0,pixel;
    int ic = LjSize/2;

    loc = ic;
    for(iy = 0; iy < LjSize; iy++)
    {
        pixel = patch[loc];
        sumVal = sumVal+kernel[iy]*pixel;
        loc += LjSize;
    }

    return sumVal;
}

float DescLJet::Conv2D(float kernel[][LjSize],const float *patch)
{
    int ix,iy,loc,col;
    //assert(kernel.size() == DescLJet::LjSize);
    float sumVal = 0,pixel;
    int ic = LjSize/2;

    loc = 0;
    for(iy = -ic; iy <= ic; iy++)
    {
        for(ix = -ic; ix <= ic; ix++)
        {
            col = ic + ix;
            pixel = patch[loc+col];
            sumVal = sumVal+kernel[iy+ic][ix+ic]*pixel;
        }
        loc = loc + LjSize;
    }

    return sumVal;
}

float DescLJet::inv_sigma(const float sigma,const int n)
{
    float v = pow(sigma,n);
    v = 1.0f/v;
    return v;
}

void DescLJet::GaussianDx(const float sigma, float kern[])
{
    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    assert (sigma > 0);
    //ranges from [-3.0*sigma,+3.0*sigma];

    int dim = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));
    dim = dim*2+1;
    float s2 = sigma * sigma;
    int c = dim / 2;
    float v = 0;

    for (int ix = -c; ix <= c; ix++)
    {
        v = DescLJet::inv_sigma(sigma, 3);
        v = v*inv_sqrt_pi2 * exp(-(ix*ix) / (2 * s2));
        kern[c+ix] = v;
        kern[c-ix] = v;

    }
    VMath::l1norm(kern, LjSize);
    return ;
}

void DescLJet::GaussianDx2_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);
    //ranges from [-3.0*sigma,+3.0*sigma];

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    //inv_sqrt_pi2 = 1.0;
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,3);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = v*((x*x)/s2 - 1.0);
            mat[c+y][c+x] = v*factor;
            //cout<<v<<" ";
        }
        //cout<<endl;
    }
    return ;

}

void DescLJet::GaussianDy2_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);
    //ranges from [-3.0*sigma,+3.0*sigma];

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,3);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = v*((y*y)/s2 - 1.0);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDxy_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);
    //ranges from [-3.0*sigma,+3.0*sigma];

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = v*x*y;
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDx3_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = x*v*(3-(x*x)/s2);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDx2y_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    //inv_sqrt_pi2 = 1;
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = y*v*(1-(x*x)/s2);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDxy2_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = x*v*(1-(y*y)/s2);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDy3_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    int y = 0, x = 0;
    float v = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = y*v*(3-(y*y)/s2);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDx4_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    float v = 0;
    int y = 0, x = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = v*(3.0-(6*x*x)/s2+(x*x*x*x)/(s2*s2));
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDx3y_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    float v = 0;
    int y = 0, x = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,7);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = x*y*v*((x*x)/s2-3.0);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDx2y2_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    float v = 0;
    int y = 0, x = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = x*y*v*(1-(y*y)/s2-(x*x)/s2+(x*x*y*y)/(s2*s2));
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDxy3_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    float v = 0;
    int y = 0, x = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,7);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = x*y*v*((y*y)/s2-3.0);
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

void DescLJet::GaussianDy4_2D(const float sigma,float mat[][LjSize])
{
    assert (sigma > 0);

    float inv_sqrt_pi2 = 1.0/sqrt(MyPI2);
    int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

    float s2 = sigma * sigma;
    float v = 0;
    int y = 0, x = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = DescLJet::inv_sigma(sigma,5);
            v =  v*inv_sqrt_pi2 * exp(-(x*x + y*y) / (2 * s2));
            v = v*(3.0-(6*y*y)/s2+(y*y*y*y)/(s2*s2));
            mat[c+y][c+x] = v*factor;
        }
    }
    return ;
}

DescLJet::~DescLJet()
{
    if(this->lj_win != NULL)
    {
        delete [] this->lj_win;
        this->lj_win = NULL;
    }

    if(this->localjet != NULL)
    {
        delete [] this->localjet;
        this->localjet = NULL;
    }
}

