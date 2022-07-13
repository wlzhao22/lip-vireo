#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "viewboard.h"
#include "descsurf.h"
#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

const int DescSURF::PatchSize  = 21;
const int DescSURF::GRID       = 4;
const int DescSURF::CHANNEL    = 4;
const int DescSURF::Radius     = 10;
const int DescSURF::iScale0    = 20;
const float DescSURF::dScale0  = 3.0; // in the paper, this has been set to 2.0f*gscale

/**
@author:    Wan-Lei Zhao
@date:      06/Jan/2012
@copyright: All rights are reserved by the author and patent holder.


Without interpolation, the performance is already acceptable. However,
interpolation will boost the performance.

**/

DescSURF::DescSURF(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];
    ///gmat = DescSURF::GaussianWeight2D(PatchSize);
    this->intImg   = NULL;

    switch(desc)
    {
    case SURF:
    {
        cout<<"Descriptor .............................. un-normalized SURF\n";
        computeDesc = & DescSURF::getSURFDesc;
        this->rChnl    = DescSURF::CHANNEL;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;

        break;
    }
    case NSURF:
    {
        cout<<"Descriptor .............................. normalized SURF\n";
        computeDesc = &DescSURF::getSURFDesc;
        this->rChnl    = DescSURF::CHANNEL;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;
        break;
    }
    case NAOD:
    {
        cout<<"Descriptor .............................. normalized AoD\n";
        computeDesc = &DescSURF::getAoDDesc;
        this->rChnl    = DescSURF::CHANNEL;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;
        break;
    }
    case AOD:
    {
        computeDesc = &DescSURF::getAoDDesc;
        this->rChnl    = DescSURF::CHANNEL;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;
        cout<<"Descriptor .............................. un-normalized AoD\n";
        break;
    }
    case NESURF:
    {
        cout<<"Descriptor .............................. normalized Extended SURF\n";
        computeDesc = &DescSURF::getESURFDesc;
        this->rChnl    = DescSURF::CHANNEL*2;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;
        break;
    }
    case ESURF:
    {
        cout<<"Descriptor .............................. Extended SURF\n";
        computeDesc = &DescSURF::getESURFDesc;
        this->rChnl    = DescSURF::CHANNEL*2;
        this->featsBin = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->pcaFeat  = new float[DescSURF::GRID*DescSURF::GRID*this->rChnl];
        this->featLen  = DescSURF::GRID*DescSURF::GRID*this->rChnl;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
    }
    }
    this->descoption = desc;
}

float DescSURF::getDx(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral(x0, (y0 - rd0), rd0, ts)
                - this->intImg->boxIntegral((x0-rd0), (y0-rd0), rd0, ts);
    return val;
}

float DescSURF::getDy(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral((x0-rd0), y0, ts, rd0) - this->intImg->boxIntegral((x0-rd0), (y0-rd0), ts, rd0);
    return val;
}

float DescSURF::getHaarX(const int x0, const int y0, const int rd0, float rMat[2])
{
     int pSize = 2*rd0, i = 0, j = 0;
     float fx = 0, fy = 0, val = 0, x = 0, y = 0;
     float lSum = 0, rSum = 0, haarX = 0;
     for(i = 0; i < pSize; i++)
     {
        y = i - rd0;
        for(j = 0; j < pSize; j++)
        {
              x = j - rd0;
              fx = x0 + rMat[0]*x - rMat[1]*y;
              fy = y0 + rMat[1]*x + rMat[0]*y;
             val = this->crntImg->getPixelBI(fx, fy);
             if(j < rd0)
             {
                lSum += val;
             }else{
                rSum += val;
             }
        }
     }
     haarX = rSum - lSum;
     return haarX;
}

float DescSURF::getHaarY(const int x0, const int y0, const int rd0, float rMat[2])
{
     int pSize = 2*rd0, i = 0, j = 0;
     float fx = 0, fy = 0, val = 0, x = 0, y = 0;
     float uSum = 0, bSum = 0, haarY = 0;
     for(i = 0; i < pSize; i++)
     {
        for(j = 0; j < pSize; j++)
        {
               x = j - rd0;
              fx = x0 + rMat[0]*x - rMat[1]*y;
              fy = y0 + rMat[1]*x + rMat[0]*y;
             val = this->crntImg->getPixelBI(fx, fy);
             if(i < rd0)
             {
                uSum += val;
             }else{
                bSum += val;
             }
        }
     }
     haarY = bSum - uSum;
     return haarY;
}

float DescSURF::gaussian(float x, float y, float sig)
{
    return 1.0f/(2.0f*PI*sig*sig) * exp( -(x*x+y*y)/(2.0f*sig*sig));
}

int DescSURF::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }

    this->buildIntImage();

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
         if(properties[i])
         nproperty += KP_PROP_SZ[i];
    }

    kpnumb = 0;
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

        (*this.*computeDesc)(crntKpt);
        if(this->descoption == NSURF || this->descoption == NAOD)
        {
            VMath::l2norm(this->featsBin, this->featLen);
        }

        if(this->_out_format == _vgg_fmrt)
        {
            this->saveDescVGG(crntKpt, this->featsBin, this->featLen, resize_rate, outStrm);
        }
        else if(this->_out_format == _vireo_fmrt)
        {
            this->saveDescVireo(crntKpt, this->featsBin, kpnumb, this->featLen, resize_rate, outStrm);
        }
        kpnumb++;
    }

    outStrm.close();
    delete this->intImg;
    this->intImg = NULL;
    return kpnumb;
}


int DescSURF::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt = NULL;
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

    ViewBoard *myview = new ViewBoard(kpnumb, PatchSize);

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }
        this->getNormDescPatch(crntKpt, this->descWin, PatchSize);
        myview->addPatch(this->descWin, PatchSize, PatchSize);
    }
    myview->saveView(dvfn);
    delete myview;
    return 0;
}

int DescSURF::getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);

    int i = 0, j = 0;
    int radius = (int)floor(Size/2.0);
    float sine = sin(keyp->ori);
    float cosine = cos(keyp->ori);
    float sc = (DescSURF::iScale0*keyp->gscale)/radius;
    float gray = 0;
    sine   = sc*sine;
    cosine = sc*cosine;
    float rPtx = 0, rPty = 0;

    int irow = 0;
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

void DescSURF::buildIntImage()
{
    unsigned int h = this->crntImg->height, w = this->crntImg->width;
    unsigned int i, j, loc;
    if(this->intImg != NULL)
    {
        delete this->intImg;
        this->intImg = NULL;
    }
    this->intImg = new IntImage(w, h);
    double val = 0;
    int x, y;

    for(i = 0; i < h; i++)
    {
        y = i;
        for(j = 0; j < w; j++)
        {
            x = j;
            loc = i*w + j;     ///[x, y]
            val = this->intImg->getPixel(x, y-1) + this->intImg->getPixel(x-1, y)
                  + this->crntImg->pix[loc] - this->intImg->getPixel(x-1, y-1);
            this->intImg->setPixel(x, y, val);
        }
    }

    return ;
}

int DescSURF::getAoDDesc(KeyPoint *crnt_kpt)
{
    float rbin = 0, cbin = 0, sc = 0, co = 0, si = 0, fx = 0, fy = 0;
    float rfactor = 0, cfactor = 0, dr = 0, dc = 0, weight = 0;
    float haar[4], rhaar[4], chaar[4], rMat[2];
    float dxs[4],  dys[4], Dx = 0, Dy = 0;
    int x = 0, y = 0, rd = 0, ts, x0, y0, rbin0, cbin0;
    int i, j, k, idxr, idxc, idx, r, c, _idx, m = 0;

    float bWdth = PatchSize/(GRID + 0.0f);
    memset(this->featsBin, 0, this->featLen*sizeof(float));

    sc = (DescSURF::iScale0*crnt_kpt->gscale)/(Radius + 0.0);
    co = cos(crnt_kpt->ori);
    si = sin(crnt_kpt->ori);
    rMat[0] = co; rMat[1] = si;
    ts = (int)round(2.0f*crnt_kpt->gscale); //2.0f*s is suggested in CVIU SURF paper
    rd = ts/2;

    for(i = 1; i < PatchSize - 1; i++)
    {
        rbin = (i + 0.0)/bWdth;
        rbin = rbin - 0.5f;
        rbin0 = (int)floor(rbin);
        rbin0 = (rbin0 < 0)?0:rbin0;
        dr = rbin - rbin0;
        y0 = i - Radius;

        for(j = 1; j < PatchSize - 1; j++)
        {
            cbin = (j + 0.0)/bWdth;
            cbin = cbin - 0.5f;
            cbin0 = (int)floor(cbin);
            cbin0 = (cbin0 < 0)?0:cbin0;
            dc = cbin - cbin0;
            x0 = j - Radius;

            fx = crnt_kpt->x + (x0*co - y0*si)*sc;
            fy = crnt_kpt->y + (x0*si + y0*co)*sc;

            x = (int)floor(fx);
            y = (int)floor(fy);
            dxs[0] = this->getHaarX(x, y, rd, rMat);
            dys[0] = this->getHaarY(x, y, rd, rMat);
            dxs[1] = dys[1] = dxs[2] = dys[2] = dxs[3] = dys[3] = 0;
            weight  = gmat[i][j];
            Dx = Dy = 0.0f;
            for(k = 0; k < 4; k++)
            {
                Dx += dxs[k];
                Dy += dys[k];
            }
            haar[0] = weight*Dx;
            haar[1] = weight*Dy;
            haar[2] = weight*fabs(Dx);
            haar[3] = weight*fabs(Dy);

            for(r = 0; r <= 1; r++)
            {
                idxr  = rbin0 + r;
                if(idxr >= 0 && idxr < DescSURF::GRID)
                {
                    rfactor = ((r == 0)?(1 - dr):dr);
                    for(m = 0; m < 4; m++)
                    {
                        rhaar[m] = haar[m]*rfactor;
                    }
                    _idx = idxr*DescSURF::GRID;
                    for(c = 0; c <= 1; c++)
                    {
                        idxc = cbin0 + c;
                        if(idxc >= 0 && idxc < DescSURF::GRID)
                        {
                            cfactor = ((c == 0)?(1 - dc):dc);
                            idx = (_idx + idxc)*this->rChnl;
                            for(m = 0; m < 4; m++)
                            {
                                chaar[m] = rhaar[m]*cfactor;
                                this->featsBin[idx+m] += chaar[m];
                            }
                        }
                    }///for(c)
                }
            }///for(r)
        }///for-loop column
    }///for-loop row
    return 1;
}


int DescSURF::getSURFDesc(KeyPoint *crnt_kpt)
{
    float rbin = 0, cbin = 0, sc = 0, co = 0, si = 0, fx = 0, fy = 0;
    float rfactor = 0, cfactor = 0, dr = 0, dc = 0, weight = 0;
    float haarX,  haarY, Dx = 0, Dy = 0;
    float haar[4], rhaar[4], chaar[4];
    int x = 0, y = 0, rd = 0, ts, x0, y0, rbin0, cbin0;
    int i, j, idxr, idxc, idx, r, c, _idx, m = 0;

    float bWdth = PatchSize/(GRID + 0.0f);
    memset(this->featsBin, 0, this->featLen*sizeof(float));

    sc = (DescSURF::iScale0*crnt_kpt->gscale)/(Radius + 0.0);
    co = cos(crnt_kpt->ori);
    si = sin(crnt_kpt->ori);
    ts = (int)round(DescSURF::dScale0*crnt_kpt->gscale);
    rd = ts/2;
    gmat = Filter::GaussianKernel2D(3.3*crnt_kpt->dscale, Radius);

    for(i = 1; i < PatchSize - 1; i++)
    {
        rbin = (i + 0.0)/bWdth;
        rbin = rbin - 0.5f;
        rbin0 = (int)floor(rbin);
        rbin0 = (rbin0 < 0)?0:rbin0;
        dr = rbin - rbin0;
        y0 = i - Radius;

        for(j = 1; j < PatchSize - 1; j++)
        {
            cbin = (j + 0.0)/bWdth;
            cbin = cbin - 0.5f;
            cbin0 = (int)floor(cbin);
            cbin0 = (cbin0 < 0)?0:cbin0;
            dc = cbin - cbin0;
            x0 = j - Radius;

            fx = crnt_kpt->x + (x0*co - y0*si)*sc;
            fy = crnt_kpt->y + (x0*si + y0*co)*sc;

            x  = (int)floor(fx);     y = (int)floor(fy);
            haarX = this->getDx(x, y, rd, ts);
            haarY = this->getDy(x, y, rd, ts);
            weight  = gmat[i][j];
            Dx = haarX*co  + haarY*si;
            Dy = haarX*-si + haarY*co;
            haar[0] = weight*Dx;
            haar[1] = weight*Dy;
            haar[2] = weight*fabs(Dx);
            haar[3] = weight*fabs(Dy);

            for(r = 0; r <= 1; r++)
            {
                idxr  = rbin0 + r;
                if(idxr >= 0 && idxr < DescSURF::GRID)
                {
                    rfactor = ((r == 0)?(1 - dr):dr);
                    for(m = 0; m < 4; m++)
                    {
                        rhaar[m] = haar[m]*rfactor;
                    }
                    _idx = idxr*DescSURF::GRID;
                    for(c = 0; c <= 1; c++)
                    {
                        idxc = cbin0 + c;
                        if(idxc >= 0 && idxc < DescSURF::GRID)
                        {
                            cfactor = ((c == 0)?(1 - dc):dc);
                            idx = (_idx + idxc)*this->rChnl;
                            for(m = 0; m < 4; m++)
                            {
                                chaar[m] = rhaar[m]*cfactor;
                                this->featsBin[idx+m] += chaar[m];
                            }
                        }
                    }///for(c)
                }
            }///for(r)
        }///for-loop column
    }///for-loop row

    Cleaner::clear2DArray(gmat);

    return 1;
}

int DescSURF::getESURFDesc(KeyPoint *crnt_kpt)
{
    float rbin = 0, cbin = 0, sc = 0, co = 0, si = 0, fx = 0, fy = 0;
    float rfactor = 0, cfactor = 0, dr = 0, dc = 0, weight = 0;
    float haarX = 0,  haarY = 0, Dx = 0, Dy = 0; ///intpol[4],
    float haar[8] = {0}, rhaar[8] = {0}, chaar[8] = {0};
    int x = 0, y = 0, rd = 0, ts, x0, y0, rbin0, cbin0;
    int i, j, idxr, idxc, idx, r, c, _idx, m = 0;

    float bWdth = PatchSize/(GRID + 0.0f);
    memset(this->featsBin, 0, this->featLen*sizeof(float));

    sc = (DescSURF::iScale0*crnt_kpt->gscale)/(Radius + 0.0);
    co = cos(crnt_kpt->ori);
    si = sin(crnt_kpt->ori);
    ts = (int)round(DescSURF::dScale0*crnt_kpt->gscale);
    rd = ts/2;

    gmat = Filter::GaussianKernel2D(3.3*crnt_kpt->dscale, Radius);

    for(i = 1; i < PatchSize-1; i++)
    {
        rbin = (i + 0.0)/bWdth;
        rbin = rbin - 0.5f;
        rbin0 = (int)floor(rbin);
        rbin0 = (rbin0 < 0)?0:rbin0;
        dr = rbin - rbin0;
        y0 = i - Radius;

        for(j = 1; j < PatchSize - 1; j++)
        {
            cbin = (j + 0.0)/bWdth;
            cbin = cbin - 0.5f;
            cbin0 = (int)floor(cbin);
            cbin0 = (cbin0 < 0)?0:cbin0;
            dc = cbin - cbin0;
            x0 = j - Radius;

            fx = crnt_kpt->x + (x0*co - y0*si)*sc;
            fy = crnt_kpt->y + (x0*si + y0*co)*sc;
             x = (int)floor(fx);
             y = (int)floor(fy);

            /**
            deltx = fx - x;
            delty = fy - y;
            intpol[0] = (1 - deltx)*(1 - delty);
            intpol[1] = deltx*(1 - delty);
            intpol[2] = (1 - deltx)*delty;
            intpol[3] = deltx*delty;
            **/

            haarX  = this->getDx(x, y, rd, ts);
            haarY  = this->getDy(x, y, rd, ts);

            weight  = gmat[i][j];
            Dx = haarX*co  + haarY*si;
            Dy = haarX*-si + haarY*co;

            if(Dy < 0)
            {
                haar[0] = weight*Dx;
                haar[2] = weight*fabs(Dx);
            }
            else
            {
                haar[1] = weight*Dx;
                haar[3] = weight*fabs(Dx);
            }
            if(Dx < 0)
            {
                haar[4] = weight*Dy;
                haar[6] = weight*fabs(Dy);
            }
            else
            {
                haar[5] = weight*Dy;
                haar[7] = weight*fabs(Dy);
            }

            for(r = 0; r <= 1; r++)
            {
                idxr  = rbin0 + r;
                if(idxr >= 0 && idxr < DescSURF::GRID)
                {
                    rfactor = ((r == 0)?(1 - dr):dr);
                    for(m = 0; m < 8; m++)
                    {
                        rhaar[m] = haar[m]*rfactor;
                    }

                    _idx = idxr*DescSURF::GRID;
                    for(c = 0; c <= 1; c++)
                    {
                        idxc = cbin0 + c;
                        if(idxc >= 0 && idxc < DescSURF::GRID)
                        {
                            cfactor = ((c == 0)?(1 - dc):dc);
                            idx = (_idx + idxc)*this->rChnl;
                            for(m = 0; m < 8; m++)
                            {
                                chaar[m] = rhaar[m]*cfactor;
                                this->featsBin[idx+m] += chaar[m];
                            }
                        }
                    }///for(c)
                }
            }///for(r)
        }///for-loop column
    }///for-loop row

    Cleaner::clear2DArray(gmat);

    return 1;
}

DescSURF::~DescSURF()
{
    vector<vector<float> >::iterator iter;
    for(iter = gmat.begin(); iter != gmat.end(); iter++)
    {
        iter->clear();
    }
    gmat.clear();
    if(intImg != NULL)
    {
        delete intImg;
        intImg = NULL;
    }
    if(this->featsBin != NULL)
    {
        delete [] featsBin;
        featsBin = NULL;
    }

    if(this->pcaFeat != NULL)
    {
        delete [] pcaFeat;
        pcaFeat = NULL;
    }
}
