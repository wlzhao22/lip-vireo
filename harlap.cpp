#include "harlap.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "vmath.h"
#include "filter.h"
#include "cleaner.h"

///#define HARLAP_DEBUG

/****Constants Initialization**/

const int HarLap::MaxOctaves  = 4;
const int HarLap::SCALES      = 7;
const float HarLap::_SIGMA    = 1.4;
const float HarLap::INITSIGMA = 0.5;
const bool  HarLap::INTERP_KEYS = false;
const float HarLap::mag      = 7.0f;

const int HarLap::BORDER     = 5;
const int HarLap::THRESH     = 60;//13
const float HarLap::dfactor  = 0.9f;
const int HarLap::maxIter    = 10;
const int HarLap::W0         = 10;
const float HarLap::err_c    = 0.90f;
const float HarLap::cvtRatio = 10.0f;

const int HarLap::DEGREE     = 10;
const int HarLap::NumOrient  = 36;
const int HarLap::DEGPERBIN  = (360 / NumOrient);
const float HarLap::NwKpThresh= 0.8;
const float HarLap::k0       = 0.06f;

/****************************************************************
*            lap = _sigma*_sigma*(dxx+dyy);
*            lap = _sigma*(dxx+dyy); //this one is better
*            The modified harris function is from [xuli CUHK],
*            the website which contains details
*            about his proposal has been removed about 1 year ago

*Xu Li's method, is not really stable
*               delta = Det*Det - 4*trace;
                delta = sqrt(delta);
                eign1 = (trace + delta)/2;
                eign2 = (trace - delta)/2;
                if(eign1 >= eign2)
                {
                    eratio = eign1/eign2;
                }
                else
                {
                    eratio = eign2 / eign1;
                }

****************************************************************/

/**implementation log************************************

I try to add following code to FindOrientation function. My purpose
is to make SIFT invariant to flip operation, however, I failed.
Try to add flip oriented mag is wong. So what's the reason
//added for flip invariant

					degree = -1*theta / PI * 180.0 + 180.0;
					index = ((int) (degree / DEGREE));
					thetas[index] += m * gmat[y + c][x + c];
17/01/2010
******************************************************/

HarLap::HarLap()
{
    cout<<"Detector ................................ HarLap\n";
    this->DETECTOR   = harlap;
    this->sel_option = THRSH;
}

bool HarLap::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    if(this->paras.find(argv[1]) != this->paras.end())
    {
        this->thresh = atof(paras["thresh"]);
    }else{
        this->thresh = THRESH;
    }
    cout<<"Thresh .................................. "<<this->thresh<<endl;

    if(paras.find(argv[2]) != paras.end())
    {
        this->fix_kp_numb = atoi(paras["topk"]);
        this->sel_option  = TOPK;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    if(paras.find(argv[3]) != paras.end())
    {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option  = DENS;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    return true;
}

vector<vector<Image *> > HarLap::BuildLOGOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves)
{
    vector<vector<Image *> > LoGOctaves;
    vector<Image *> gxxScales;
    vector<Image *> gyyScales;
    vector<Image *> LoGScales;

    for (unsigned int i = 0; i < GxxOctaves.size(); i++)
    {
        gxxScales = GxxOctaves[i];
        gyyScales = GyyOctaves[i];
        LoGScales = BuildLoGScales(gxxScales,gyyScales);
        LoGOctaves.push_back(LoGScales);
    }

    return LoGOctaves;
}

vector<Image* > HarLap::BuildLoGScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales)
{
    vector<Image* > LoGScales;
    Image *crntDxxImage,*crntDyyImage,*crntLoGImage;
    float dxx, dyy, _sigma, lap;
    int x, y, width, height;
    unsigned int nscale;

    for(nscale = 0; nscale < gxxScales.size(); nscale++)
    {
        crntDxxImage = gxxScales[nscale];
        crntDyyImage = gyyScales[nscale];
        width        = crntDxxImage->width;
        height       = crntDxxImage->height;
        crntLoGImage = new Image(width, height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dxx    = crntDxxImage->getPixel(x, y);
                dyy    = crntDyyImage->getPixel(x, y);
                _sigma = HarLap::_SIGMA * pow(2.0, nscale / (float) SCALES);
                lap    = _sigma*(dxx+dyy);
                ///lap = _sigma*_sigma*(dxx+dyy);
                ///lap = lap<0?(-1*lap):lap;

                crntLoGImage->setPixel(x, y, lap);
            }
        }
        LoGScales.push_back(crntLoGImage);
    }

    return LoGScales;
}

vector<vector<Image *> > HarLap::BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves)
{
    vector<vector<Image *> > HarrisOctaves;
    vector<Image *> gxScales, gyScales;
    vector<Image *> HarrisScales;

    for (unsigned int i = 0; i < GxOctaves.size(); i++)
    {
        gxScales = GxOctaves[i];
        gyScales = GyOctaves[i];
        HarrisScales = BuildHarrisScales(gxScales,gyScales);
        HarrisOctaves.push_back(HarrisScales);
    }

    return HarrisOctaves;
}

vector<vector<Image *> > HarLap::BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves, vector<vector<ImageSymmMat *> > & pmatOctaves, vector<vector<Image *> > & GyOctaves)
{
    vector<vector<Image *> > HarrisOctaves;
    vector<Image *> gxScales, gyScales;
    vector<Image *> HarrisScales;

    for (unsigned int i = 0; i < GxOctaves.size(); i++)
    {
        gxScales = GxOctaves[i];
        gyScales = GyOctaves[i];
        vector<ImageSymmMat*> paraMats;
        HarrisScales = BuildHarrisScales(gxScales, paraMats, gyScales);
        HarrisOctaves.push_back(HarrisScales);
        pmatOctaves.push_back(paraMats);
    }

    return HarrisOctaves;
}

vector<Image* > HarLap::BuildHarrisScales(vector<Image *>  &gxScales,vector<Image *> &gyScales)
{
    vector<Image* > HarrisScales;
    Image *crntDxImg = NULL, *crntDyImg = NULL, *crntDxDyImg = NULL, *crntHarrisImg = NULL;
    float dx2 = 0, dy2 = 0, Det = 0, trace = 0, harris = 0;
    float _sigma = 0, dsigma = 0, dsigma2 = 0;
    int x = 0, y = 0, width = 0, height = 0;
    float dx = 0, dy = 0, dxdy = 0;
    unsigned int nscale = 0;
    vector<float> Gkern;

    if(gxScales.size() == 0 || gyScales.size() == 0)
    {
        return HarrisScales;
    }

    Image *tmpimg = new Image(gxScales[0]->width, gxScales[0]->height);

    for(nscale = 0; nscale < gxScales.size(); nscale++)
    {
        crntDxImg = gxScales[nscale];
        crntDyImg = gyScales[nscale];

        width   = crntDxImg->width;
        height  = crntDxImg->height;
        crntDxDyImg = new Image(crntDxImg->width,crntDxImg->height);
        _sigma  = HarLap::_SIGMA * pow(2.0,  (nscale+1.0) / (float) SCALES);
        dsigma  = _sigma*HarLap::dfactor;
        dsigma2 = dsigma*dsigma;

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx   = crntDxImg->getPixel(x,y);
                dy   = crntDyImg->getPixel(x,y);
                dxdy = dx*dy*dsigma2;
                crntDxDyImg->setPixel(x,y,dxdy);
            }
        }

        crntDxImg->exp(2);
        crntDyImg->exp(2);
        crntDxImg->multiply(dsigma2);
        crntDyImg->multiply(dsigma2);
        Gkern = Filter::GaussianKernel1D(_sigma);
        Filter::Convolve1DWidth(Gkern,  crntDxImg,  tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,       crntDxImg);
        Filter::Convolve1DWidth(Gkern,  crntDyImg,  tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,       crntDyImg);
        Filter::Convolve1DWidth(Gkern,  crntDxDyImg,tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,       crntDxDyImg);
        crntHarrisImg = new Image(crntDxImg->width, crntDxImg->height);
        Gkern.clear();

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx2   = crntDxImg->getPixel(x, y);
                dy2   = crntDyImg->getPixel(x, y);
                dxdy  = crntDxDyImg->getPixel(x, y);
                Det   = dx2*dy2 - dxdy*dxdy;
                trace = dx2+dy2;

                trace = trace*trace;
                harris = Det - HarLap::k0*trace;
                crntHarrisImg->setPixel(x, y, harris);
            }
        }

        HarrisScales.push_back(crntHarrisImg);
        delete crntDxDyImg;
    }

    delete tmpimg;
    return HarrisScales;
}

vector<Image* > HarLap::BuildHarrisScales(vector<Image *>  &gxScales, vector<ImageSymmMat*> &paraMats, vector<Image *> &gyScales)
{
    vector<Image* > HarrisScales;
    Image *crntDxDyImage = NULL, *crntHarrisImage = NULL;
    Image *crntDxImage = NULL, *crntDyImage = NULL;
    float dx2 = 0, dy2 = 0, Det = 0, trace = 0;
    float dx = 0, dy = 0, dxdy = 0, harris = 0;
    float _sigma = 0, dsigma = 0, dsigma2 = 0;
    int   x = 0, y = 0, width = 0, height = 0;
    unsigned int nscale = 0;
    vector<float> Gkern;

    if(gxScales.size() == 0 || gyScales.size() == 0)
    {
        return HarrisScales;
    }

    Image *tmpimg = new Image(gxScales[0]->width,gxScales[0]->height);
    ImageSymmMat *crnt_paraMat;

    for(nscale = 0; nscale < gxScales.size(); nscale++)
    {
        crntDxImage = gxScales[nscale];
        crntDyImage = gyScales[nscale];

        width   = crntDxImage->width;
        height  = crntDxImage->height;
        crntDxDyImage = new Image(crntDxImage->width,crntDxImage->height);
        _sigma  = HarLap::_SIGMA * pow(2.0,  (nscale+1.0f) / (float) SCALES);
        dsigma  = _sigma*HarLap::dfactor;
        dsigma2 = dsigma*dsigma;

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx = crntDxImage->getPixel(x, y);
                dy = crntDyImage->getPixel(x, y);
                dxdy = dx*dy*dsigma2;
                crntDxDyImage->setPixel(x, y, dxdy);
            }
        }

        crntDxImage->exp(2);
        crntDyImage->exp(2);
        crntDxImage->multiply(dsigma2);
        crntDyImage->multiply(dsigma2);
        Gkern = Filter::GaussianKernel1D(_sigma);
        Filter::Convolve1DWidth(Gkern,  crntDxImage, tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg, crntDxImage);
        Filter::Convolve1DWidth(Gkern,  crntDyImage, tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,  crntDyImage);
        Filter::Convolve1DWidth(Gkern,  crntDxDyImage, tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg, crntDxDyImage);
        Gkern.clear();

        crntHarrisImage = new Image(crntDxImage->width, crntDxImage->height);
        crnt_paraMat    = new ImageSymmMat(crntDxImage->width, crntDxImage->height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx2  = crntDxImage->getPixel(x,y);
                dy2  = crntDyImage->getPixel(x,y);
                dxdy = crntDxDyImage->getPixel(x,y);
                Det  = dx2*dy2 - dxdy*dxdy;
                trace = dx2 + dy2;
                crnt_paraMat->pSet(x, y, dx2, dxdy, dy2);
                harris = (dsigma2*Det)/trace;
                crntHarrisImage->setPixel(x, y, harris);
            }
        }

        HarrisScales.push_back(crntHarrisImage);
        paraMats.push_back(crnt_paraMat);
        delete crntDxDyImage;
    }

    delete tmpimg;

    return HarrisScales;
}

bool HarLap::isSpatialPeak(Image *image, const int x, const int y)
{
    float center =image->getPixel(x, y);
    int xi = 0, yi = 0;
    float p2 = 0;
    if(center > 0)
    for (yi = y - 1; yi <= y + 1; yi++)
    {
        for (xi = x - 1; xi <= x + 1; xi++)
        {
            p2 = image->getPixel(xi, yi);
            if (center < p2)
                return false;
        }
    }else{
        return false;
    }
    return true;
}

bool HarLap::isScalePeak(Image * aimage, Image *bimage, Image *cimage, float &logVal, const int x, const int y)
{
    assert(aimage);
    assert(bimage);
    assert(cimage);

    Image *ims[3] = {NULL};

    ims[0] = aimage;
    ims[1] = bimage;
    ims[2] = cimage;

    float center = bimage->getPixel(x, y);
    float p2 = 0;
    logVal = center;

    if (center > 0.0)
    {
        for (int si = 0; si < 3; si++)
        {
            p2 = ims[si]->getPixel(x, y);
            if (center < p2)
                return false;
        }
        return true;
    }
    else
    {
        /**/
        for (int si = 0; si < 3; si++)
        {
            p2 = ims[si]->getPixel(x, y);
            if (center > p2)
                return false;
        }
        return true;
        /**/
        //return false;
    }
}

bool HarLap::adaptAffine(vector<vector<Image *> > &HarrisOctaves, vector<vector<ImageSymmMat *> > & pmatOctaves)
{
    int ioctave, iter, is, x, y, xi, yi, x_flr, y_flr, x_m, y_m;
    float a, b, c, xHarris, gray, dx, dy, xn, yn;
    float m[2][2], u[2][2], u_tmp[2][2];
    bool CNVG;
    unsigned int i;

    vector<ImageSymmMat *> pMats;
    ImageSymmMat *crnt_pMat;
    KeyPoint *crntKp;
    vector<Image *> pHarris;
    Image *crntImg;

    for (i = 0; i < this->kps.size(); i++)
    {
        crntKp = this->kps[i];

        is = crntKp->scale;
        x = xn = (int)round(crntKp->sx);
        y = yn = (int)round(crntKp->sy);
        ioctave = crntKp->octave;
        vector<Image *> &pHarris = HarrisOctaves[ioctave];
        crntImg = pHarris[is];
        vector<ImageSymmMat *> &pMats = pmatOctaves[ioctave];
        crnt_pMat = pMats[is];
        iter = 0;
        CNVG = false;
        crnt_pMat->pGet(x, y, a, b, c);
        CNVG = VMath::inv_sqrtSymMat(a, b, c, u, HarLap::err_c);
        VMath::normMat(u,u);
        CNVG = true;

        while(iter < HarLap::maxIter && !CNVG)
        {
            xHarris = 0;
            x_m = 0;
            y_m = 0;

            for(yi = -HarLap::W0; yi <= HarLap::W0; yi++)
            {
                for(xi = -HarLap::W0; xi <= HarLap::W0; xi++)
                {
                    xn    = x + u[0][0]*xi + u[0][1]*yi;
                    yn    = y + u[1][0]*xi + u[1][1]*yi;
                    x_flr = floor(xn);
                    y_flr = floor(yn);
                    dx    = xn - x_flr;
                    dy    = yn - y_flr;
                    gray  = crntImg->getPixel(x_flr, y_flr)*(1-dx)*(1-dy);
                    gray += crntImg->getPixel(x_flr+1,y_flr)*dx*(1-dy);
                    gray += crntImg->getPixel(x_flr, y_flr+1)*(1-dx)*dy;
                    gray += crntImg->getPixel(x_flr+1,y_flr+1)*dx*dy;

                    if(gray > xHarris)
                    {
                        xHarris = gray;
                        x_m = xi;
                        y_m = yi;
                    }
                }///end for
            }///end for

            if(x_m != 0 || y_m != 0 )
            {
                xn = x + u[0][0]*x_m + u[0][1]*y_m;
                yn = y + u[1][0]*x_m + u[1][1]*y_m;
                xn = xn < 0?0:xn;
                yn = yn < 0 ? 0 : yn;
                x_flr = floor(xn);
                y_flr = floor(yn);
                dx = xn - x_flr;
                dy = yn - y_flr;

                crnt_pMat->pGet(x_flr, y_flr, a, b, c);
                m[0][0] = a*(1-dx)*(1-dy);
                m[0][1] = b*(1-dx)*(1-dy);
                m[1][1] = c*(1-dx)*(1-dy);
                crnt_pMat->pGet(x_flr+1, y_flr, a, b, c);
                m[0][0] += a*dx*(1-dy);
                m[0][1] += b*dx*(1-dy);
                m[1][1] += c*dx*(1-dy);
                crnt_pMat->pGet(x_flr, y_flr+1,  a, b, c);
                m[0][0] += a*dy*(1-dx);
                m[0][1] += b*dy*(1-dx);
                m[1][1] += c*dy*(1-dx);
                crnt_pMat->pGet(x_flr+1, y_flr+1,  a, b, c);
                m[0][0] += a*dx*dy;
                m[0][1] += b*dx*dy;
                m[1][1] += c*dx*dy;
                m[1][0] = m[0][1];

                CNVG = VMath::inv_sqrtSymMat(m[0][0], m[0][1], m[1][1], m, HarLap::err_c);
                if(m[0][0] == 1 && m[1][1] == 1 && m[0][1] == 0 && m[1][0] == 0)
                {
                    xn = x;
                    yn = y;
                    CNVG = true;
                }
                else
                {
                    u_tmp[0][0] = m[0][0]*u[0][0] + m[0][1]*u[1][0];
                    u_tmp[0][1] = m[0][0]*u[0][1] + m[0][1]*u[1][1];
                    u_tmp[1][0] = m[1][0]*u[0][0] + m[1][1]*u[1][0];
                    u_tmp[1][1] = m[1][0]*u[0][1] + m[1][1]*u[1][1];
                    VMath::normMat(u_tmp,u);
                    iter++;
                }
                x = (int)round(xn);
                y = (int)round(yn);
            }
            else
            {
                break;
            }

        }

        crntKp->sx = xn;
        crntKp->sy = yn;
        crntKp->x  = (int)round(xn * pow(2.0, ioctave));
        crntKp->y  = (int)round(yn * pow(2.0, ioctave));

        if(!VMath::normMat(u, u, crntKp->e1, crntKp->e2, HarLap::cvtRatio))
        {
            crntKp->KP = false;
        }

        if(crntKp->x > (this->crntimg->width - HarLap::BORDER) || crntKp->y > (this->crntimg->height- HarLap::BORDER))
        {
            crntKp->KP = false;
        }

        crntKp->a = u[0][0];
        crntKp->b = u[0][1];
        crntKp->c = u[1][1];
    }
    return true;
}

vector<KeyPoint *> HarLap::FindPeaksScales(const int octave, vector<Image*> HarImages, vector<Image *> & LoGImages)
{
    vector<KeyPoint *> peaks;

    Image * kpfound = new Image(LoGImages[0]->width, LoGImages[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0, logVal = 0;
    unsigned int s = 0;
    int x = 0, y = 0;

    float contrast_thresh = (float) THRESH;

    for (s = 1; s < HarImages.size()-1; s++)
    {
        for (y = BORDER; y < (HarImages[0]->height - BORDER); y++)
        {
            for (x = BORDER; x < (HarImages[0]->width - BORDER); x++)
            {
                funcVal = HarImages[s]->getPixel(x, y);

                if (funcVal <= contrast_thresh)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                if (!isSpatialPeak(HarImages[s], x, y))
                    continue;

                if(!isScalePeak(LoGImages[s-1], LoGImages[s], LoGImages[s+1], logVal, x, y))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                fs  = fs < 0?0:fs;

                KeyPoint * peak = new KeyPoint();

                peak->x = (int)round(fx * pow(2.0, octave));
                peak->y = (int)round(fy * pow(2.0, octave));

                peak->dscale    = HarLap::_SIGMA * pow(2.0, octave + fs / (float) SCALES);
                peak->iscale    = peak->dscale*HarLap::mag;
                peak->funcVal   = funcVal;
                peak->img_width = this->crntimg->width;
                peak->ori       = 0;
                peak->scale     = s;
                peak->fscale    = fs;
                peak->gscale    = octave + fs/HarLap::SCALES;
                peak->sx        = fx;
                peak->sy        = fy;
                peak->octave    = octave;
                peak->KP        = true;
                peak->div       = logVal > 0?1:-1;

                leveli_kps.push_back(peak);
                kpfound->setPixel((int)(fx + 0.5), (int)(fy + 0.5), 1);
            }
        }
        peaks.insert(peaks.begin(),leveli_kps.begin(),leveli_kps.end());
        leveli_kps.clear();
    }
    delete kpfound;
    return peaks;
}

vector<KeyPoint *> HarLap::FindPeaksScales(const int octave, vector<Image*> HarImages, vector<ImageSymmMat*> &paraMats, vector<Image *> & LoGImages)
{
    vector<KeyPoint *> peaks;

    Image * kpfound = new Image(LoGImages[0]->width, LoGImages[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0, logVal = 0;
    int x = 0, y = 0, ux = 0, uy = 0;
    unsigned int s = 0;
    float contrast_thresh = this->thresh;

    for (s = 1; s < HarImages.size()-1; s++)
    {
        for (y = BORDER; y < (HarImages[0]->height - BORDER); y++)
        {
            for (x = BORDER; x < (HarImages[0]->width - BORDER); x++)
            {
                funcVal = HarImages[s]->getPixel(x, y);

                if (funcVal <= contrast_thresh)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                if (!isSpatialPeak(HarImages[s], x, y))
                    continue;

                if(!isScalePeak(LoGImages[s-1],LoGImages[s],LoGImages[s+1],logVal,x,y))
                    continue;

                /**useful
                if (isEdgePeak(x, y, LoGImages[s]))
                	continue;
                **/

                fx = x;
                fy = y;
                fs = s;

                fs  = fs < 0?0:fs;
                cout<<funcVal<<endl;

                KeyPoint * peak = new KeyPoint();

                peak->x = (int)round(fx * pow(2.0, octave));
                peak->y = (int)round(fy * pow(2.0, octave));

                peak->dscale = HarLap::_SIGMA * pow(2.0, octave + fs / (float) SCALES);
                peak->iscale = peak->dscale*HarLap::mag;
                peak->funcVal = funcVal;

                peak->scale  = s;
                peak->fscale = fs;
                peak->gscale = octave + fs/HarLap::SCALES;
                peak->sx = fx;
                peak->sy = fy;
                peak->octave = octave;
                peak->KP = true;
                peak->ori = 0;
                ux = floor(fx+0.5);
                uy = floor(fy+0.5);
                //cout<<peak->dscale<<endl;
                /**
                paraMats[us]->pGet(ux, uy, m[0][0], m[0][1], m[1][1]);
                m[1][0] = m[0][1];
                VMath::inv_sqrtSymMat(m[0][0], m[0][1], m[1][1],m,0);
                peak->a =m[0][0];
                peak->b = m[0][1];
                peak->c = m[1][1];
                **/

                leveli_kps.push_back(peak);
                kpfound->setPixel(ux, uy, 1);
            }
        }
        peaks.insert(peaks.begin(),leveli_kps.begin(),leveli_kps.end());
        leveli_kps.clear();
    }
    delete kpfound;
    //cout<<"end finding\n";
    return peaks;
}

vector<KeyPoint *> HarLap::FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves)
{
    vector<KeyPoint * > newkps;
    vector<vector<float> > gmat;
    KeyPoint *newkp = NULL;

    float sigma = 0.0, m = 0.0F, theta = 0, degree = 0;
    int   c = 0, index = 0, j = 0;
    int   indexa, indexb, indexc;
    float thetaa, thetab, thetac;
    unsigned int i = 0;
    float maxval = 0, maxp = 0;

    bool valid = true;
    float thetas[NumOrient] = {0};
    float weight = 0.0f;

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *HarLap::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size()/2;

        for(j = 0; j < NumOrient; j++)
        {
            thetas[j] = 0;
        }
        for (int y = -c; y <= c; y++)
        {
            for (int x = -c; x <= c; x++)
            {
                if (sqrt((float) x*x + y*y) > 3.0 * sigma)
                    continue;

                valid = Filter::GetPixOrientation((int) (kps[i]->sx + x + 0.5), (int) (kps[i]->sy + y + 0.5),
                                                  GOctaves[kps[i]->octave][kps[i]->scale], m, theta);

                if(valid)
                {
                    degree = theta / PI * 180.0 + 180.0;
                    index = ((int) (degree / DEGREE));
                    index = index%NumOrient;
                    weight = gmat[y + c][x + c]*m;
                    thetas[index] += weight;
                }
            }
        }
        vector<vector<float> >::iterator it;
        vector<float> crntvect;
        for(it = gmat.begin(); it != gmat.end(); it++)
        {
            crntvect = *it;
            crntvect.clear();
        }
        gmat.erase(gmat.begin(),gmat.end());

        for (int j = 0; j < 6; j++)
            VMath::SmoothHistogram(thetas, NumOrient);

        maxval = VMath::maxVec(thetas, NumOrient);

        for (int j = 0; j < NumOrient; j++)
        {
            if (thetas[j] < maxval * NwKpThresh)
                continue;

            indexa = VMath::mod(j - 1, NumOrient);
            indexb = j;
            indexc = VMath::mod(j + 1, NumOrient);
            thetaa = thetas[indexa];
            thetab = thetas[indexb];
            thetac = thetas[indexc];
            if (!(thetab > thetaa && thetab > thetac))
                continue;

            maxp = VMath::getMaxParabola(-1, thetaa, 0, thetab, 1, thetac);

            if(thetas[j] == maxval)
            {
                kps[i]->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NumOrient - PI;
                ///kps[i]->flip = findFlip(thetas, j, NumOrient, rate);
            }
            else
            {

                newkp = new KeyPoint();
                memcpy(newkp, kps[i], sizeof(KeyPoint));
                if(/**/this->mydesc == ERIFT || this->mydesc == NERIFT || /**/
                        this->mydesc == NSPIN || this->mydesc == SPIN)
                {
                    newkp->KP = false;
                }
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NumOrient - PI;
                ///newkp->flip = findFlip(thetas, j, NumOrient, rate);
                newkps.push_back(newkp);
            }

        }
    }
    return newkps;
}

bool HarLap::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn)
{
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    AbstractDetector::releaseKpList(this->kps);
    Image *oriImg = NULL, *tmpImg = NULL;

    vector<vector<Image*> > GxOctaves, GyOctaves;
    vector<vector<Image*> > GxxOctaves, GyyOctaves;
    vector<vector<Image *> > LoGOctaves;
    vector<vector<Image *> > HarrisOctaves;
    vector<vector<Image *> > GaussOctaves;
    vector<vector<ImageSymmMat *> > pmatOctaves;
    vector<Image *> Scales;
    vector<KeyPoint *> peaks;
    unsigned int ioctave;

    int minsz = (this->crntimg->width>this->crntimg->height)?this->crntimg->height:this->crntimg->width;
    if(minsz < 50)
    {
        tmpImg = Image::doubleSizeImage(this->crntimg);
        delete this->crntimg;
        this->crntimg = tmpImg;
        this->resize_rate = 0.5;
        this->RESIZE = true;
    }
    else if(minsz > imgSzBound)
    {
        this->resize_rate = 1.0f;
        while(minsz > imgSzBound)
        {
            tmpImg = Image::halfSizeImage(this->crntimg);
            delete this->crntimg;
            this->crntimg = tmpImg;
            minsz =  (this->crntimg->width>this->crntimg->height)?this->crntimg->height:this->crntimg->width;
            this->resize_rate = this->resize_rate*2;
            this->RESIZE = true;
        }
    }
    else
    {
        this->resize_rate = 1.0f;
        this->RESIZE = false;
    }

    oriImg = this->crntimg;
    tmpImg = Filter::ScaleImage(this->crntimg,0,_SIGMA,INITSIGMA);
    this->crntimg = tmpImg;

    this->intImg  = AbstractDetector::buildIntImage(this->crntimg);

    this->BuildOctaves(this->crntimg, _Dx_,   GxOctaves,    HarLap::_SIGMA, HarLap::SCALES);
    this->BuildOctaves(this->crntimg, _Dy_,   GyOctaves,    HarLap::_SIGMA, HarLap::SCALES);
    this->BuildOctaves(this->crntimg, _Dxx_,  GxxOctaves,   HarLap::_SIGMA, HarLap::SCALES);
    this->BuildOctaves(this->crntimg, _Dyy_,  GyyOctaves,   HarLap::_SIGMA, HarLap::SCALES);
    this->BuildOctaves(this->crntimg, _Conv2,GaussOctaves, HarLap::_SIGMA, HarLap::SCALES);

    LoGOctaves = this->BuildLOGOctaves(GxxOctaves, GyyOctaves);
    HarrisOctaves = this->BuildHarrisOctaves(GxOctaves, pmatOctaves, GyOctaves);

    Cleaner::releaseOctaves(GxxOctaves);
    Cleaner::releaseOctaves(GyyOctaves);
    Cleaner::releaseOctaves(GyOctaves);
    Cleaner::releaseOctaves(GxOctaves);

    for(ioctave = 0; ioctave < LoGOctaves.size(); ioctave++)
    {
        peaks = this->FindPeaksScales(ioctave,HarrisOctaves[ioctave],LoGOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    Cleaner::releaseOctaves(LoGOctaves);

    /**/
    peaks = this->FindOrientByGrad(kps, GaussOctaves);
    if(peaks.size() > 0)
    {
        this->kps.insert(kps.begin(),peaks.begin(),peaks.end());
        peaks.clear();
    }

    stable_sort(kps.begin(),kps.end(),KeyPoint::keypCompF);

    /**/
    switch(this->sel_option)
    {
    case 0:
    {
        AbstractDetector::topkSelect(kps, this->fix_kp_numb);
        break;
    }
    case 1:
    {
        AbstractDetector::topkEqDnSelect(kps, this->crntimg->width, this->crntimg->height, this->fix_kp_numb);
        break;
    }
    default:
    {
        break;
    }
    }

    if(this->AFF_OUT)
    {
        this->adaptAffine(HarrisOctaves, pmatOctaves);
    }

    Cleaner::releaseOctaves(pmatOctaves);
    Cleaner::releaseOctaves(HarrisOctaves);
    Cleaner::releaseOctaves(GaussOctaves);

    if(strcmp(dstfn,""))
    {
        (this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
    }

    delete this->crntimg;
    this->crntimg = oriImg;

    if(strcmp(descfn, "")&& this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }

    if(strcmp(dvfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }

    delete this->crntimg;
    delete this->intImg;
    this->intImg  = NULL;
    this->crntimg = NULL;
    return true;
}

void HarLap::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new HarLap();

    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 4; i++)
    {
        sprintf(vsrc,  "/home/wlzhao/datasets/trec03/img%d.jpg",   i);
        sprintf(vdst,  "/home/wlzhao/datasets/trec03/img%d.keys",  i);
        sprintf(vdesc, "/home/wlzhao/datasets/trec03/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }

}
