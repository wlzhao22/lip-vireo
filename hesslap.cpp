#include "hesslap.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <cmath>

#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

/**constant initialization**/

const int HessLap::MaxOctaves    = 4;
const int HessLap::SCALES        = 6;
const float HessLap::_SIGMA      = 2.0;
const float HessLap::INITSIGMA   = 0.5;
const bool  HessLap::INTERP_KEYS = true;

const float HessLap::mag         = 4.0;
const int HessLap::BORDER        = 5;
const int HessLap::THRESH        = 350;

const int HessLap::DEGREE        = 10;
const int HessLap::NumOrient     = 36;
const int HessLap::DEGPERBIN     = (360 / NumOrient);
const float HessLap::NwKpThresh  = 0.8;

const float HessLap::cvtRatio    = 10.0f;

const int HessLap::maxIter       = 10;
const int HessLap::W0            = 10;
const float HessLap::err_c       = 0.80f;

/****************************************
* I try use a dynamic threshold, however, no significant improvement I can get,
even, it actually results in a few bigger number of interest points

@date Mar.2nd.2009

const float HessLap::thresh_ratio = 2.5f;
**********************************************/

/**
*            lap = _sigma*_sigma*(dxx+dyy);
*            lap = _sigma*(dxx+dyy); //this one is better
*
**/

HessLap::HessLap()
{
    cout<<"Detector ................................ HessLap\n";
    this->DETECTOR   = hesslap;
    this->sel_option = THRSH;
}

bool HessLap::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    if(this->paras.find(argv[1]) != this->paras.end())
    {
        this->thresh = atof(paras["thresh"]);
    }else{
        this->thresh = THRESH;
    }
    cout<<"Thresh .................................. "<<this->thresh<<endl;
    if(this->paras.find(argv[2]) != this->paras.end())
    {
        this->fix_kp_numb = atoi(paras["topk"]);
        this->sel_option  = TOPK;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    if(this->paras.find(argv[3]) != this->paras.end())
    {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option  = DENS;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    return true;
}

float HessLap::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
{

    float Dp[3] = {0}; // first derivative of D with respect to x, y, s

    float Dpp[3][3] = {{0}}; // Hessian of D

    Dp[0] = (DI[s]->getPixel(x+1, y) - DI[s]->getPixel(x-1, y)) / 2.0; // Dx
    Dp[1] = (DI[s]->getPixel(x, y+1) - DI[s]->getPixel(x, y-1)) / 2.0; // Dy
    Dp[2] = (DI[s+1]->getPixel(x, y) - DI[s-1]->getPixel(x, y)) / 2.0; // Ds

    // Dxx
    Dpp[0][0] = (DI[s]->getPixel(x+1, y) + DI[s]->getPixel(x-1, y)
                 - 2.0 * DI[s]->getPixel(x, y));

    // Dyy
    Dpp[1][1] = (DI[s]->getPixel(x, y+1) + DI[s]->getPixel(x, y-1)
                 - 2.0 * DI[s]->getPixel(x, y));

    // Dzz
    Dpp[2][2] = (DI[s+1]->getPixel(x, y) + DI[s-1]->getPixel(x, y)
                 - 2.0 * DI[s]->getPixel(x, y));


    // Dxy = Dyx
    Dpp[0][1] = Dpp[1][0] = (DI[s]->getPixel(x+1, y+1) - DI[s]->getPixel(x-1, y+1)
                             - DI[s]->getPixel(x+1, y-1) + DI[s]->getPixel(x-1, y-1)) / 4.0;

    // Dxs = Dsx
    Dpp[0][2] = Dpp[2][0] = (DI[s+1]->getPixel(x+1, y) - DI[s+1]->getPixel(x-1, y)
                             - DI[s-1]->getPixel(x+1, y) + DI[s-1]->getPixel(x-1, y)) / 4.0;

    // Dys = Dsy
    Dpp[1][2] = Dpp[2][1] = (DI[s+1]->getPixel(x, y+1) - DI[s+1]->getPixel(x, y-1)
                             - DI[s-1]->getPixel(x, y+1) + DI[s-1]->getPixel(x, y-1)) / 4.0;


    float invDpp[3][3];

    VMath::mInv33(Dpp, invDpp);

    // Solve for delta positions
    *dx = 0;
    for (int i = 0; i < 3; i++)
        *dx -= invDpp[0][i] * Dp[i];

    *dy = 0;
    for (int i = 0; i < 3; i++)
        *dy -= invDpp[1][i] * Dp[i];

    *ds = 0;
    for (int i = 0; i < 3; i++)
        *ds -= invDpp[2][i] * Dp[i];

    //printf("Interp: %f %f %f\n", *dx, *dy, *ds);

    float val = DI[s]->getPixel(x, y);
    val += 0.5 * (Dp[0] * *ds + Dp[1] * *dy + Dp[2] * *ds);

    return fabs(val);
}

bool HessLap::InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal)
{
    bool addkey = true;
    int moves_left = 5;
    int tx = x;
    int ty = y;
    int ts = s;

    float dx, dy, ds,val;
    bool updated;
    float contrast_thresh = 0.8*this->thresh;//(float)NUM_SCALES;

    do
    {
        moves_left--;
        updated = false;

        val = InterpKeyStep(tx, ty, ts, LoGImages, &dx, &dy, &ds);
        if (val < contrast_thresh)
        {
            addkey = false;
            continue;
        }


        if (dx > 0.6 && tx < LoGImages[0]->width - 3)
        {
            tx++;
            updated = true;
        }
        else if (dx < -0.6 && tx > 3)
        {
            tx--;
            updated = true;
        }


        if (dy > 0.6 && ty < LoGImages[0]->height - 3)
        {
            ty++;
            updated = true;
        }
        else if (dy < -0.6 && ty > 3)
        {
            ty--;
            updated = true;
        }

    }
    while (moves_left > 0 && updated);
    *dogVal = val;

    if (addkey && fabs(dx) < 1.5 && fabs(dy) < 1.5 && fabs(ds) < 1.5)
    {
        *fx = tx + dx;
        *fy = ty + dy;
        *fs = ts + ds;
        return true;
    }

    return false;
}

vector<vector<Image *> > HessLap::BuildLOGOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves)
{
    vector<vector<Image *> > LoGOctaves;
    vector<Image *> gxxScales;
    vector<Image *> gyyScales;

    for (unsigned int i = 0; i < GxxOctaves.size(); i++)
    {
        gxxScales = GxxOctaves[i];
        gyyScales = GyyOctaves[i];
        vector<Image *> LoGScales = BuildLoGScales(gxxScales, gyyScales);
        LoGOctaves.push_back(LoGScales);
    }

    return LoGOctaves;
}

vector<Image* > HessLap::BuildLoGScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales)
{
    vector<Image* > LoGScales;
    Image *crntDxxImage,*crntDyyImage,*crntLoGImage;
    float dxx,dyy;
    float _sigma,lap;
    int x, y, width, height;
    unsigned int nscale;

    for(nscale = 0; nscale < gxxScales.size(); nscale++)
    {
        crntDxxImage = gxxScales[nscale];
        crntDyyImage = gyyScales[nscale];
        width = crntDxxImage->width;
        height = crntDxxImage->height;
        crntLoGImage = new Image(crntDxxImage->width,crntDxxImage->height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dxx = crntDxxImage->getPixel(x,y);
                dyy = crntDyyImage->getPixel(x,y);
                _sigma = HessLap::_SIGMA * pow(2.0,  nscale / (float) SCALES);
                lap = _sigma*_sigma*(dxx+dyy);
                //lap = _sigma*(dxx+dyy);
                //lap = lap<0?(-1*lap):lap;
                crntLoGImage->setPixel(x,y,lap);
            }
        }
        LoGScales.push_back(crntLoGImage);
    }

    return LoGScales;
}

vector<vector<Image *> > HessLap::BuildHessOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves, vector<vector<Image *> > &LoGOctaves)
{
    vector<Image *> gxxScales, gyyScales, gxyScales;
    vector<vector<Image *> > HessOctaves;

    for (unsigned int i = 0; i < GxxOctaves.size(); i++)
    {
        gxxScales = GxxOctaves[i];
        gyyScales = GyyOctaves[i];
        gxyScales = GxyOctaves[i];
        vector<Image*> LoGScales;
        vector<Image *> HessScales;
        HessScales = BuildHessScales(gxxScales,gyyScales,gxyScales,LoGScales);
        HessOctaves.push_back(HessScales);
        LoGOctaves.push_back(LoGScales);
    }

    return HessOctaves;
}

vector<Image* > HessLap::BuildHessScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales,vector<Image *> &gxyScales,vector<Image*> &LoGScales)
{
    vector<Image* > HessScales;
    Image *crntDxxImg = NULL, *crntDyyImg = NULL, *crntDxyImg = NULL, *crntHessImg = NULL, *crntLoGImg = NULL;
    float _sigma = 0, det = 0, trace2 = 0, lap = 0;
    int x = 0, y = 0, width = 0, height = 0;
    float dxx = 0, dyy = 0 ,dxy = 0;
    unsigned int nscale = 0;

    for(nscale = 0; nscale < gxxScales.size(); nscale++)
    {
        crntDxxImg = gxxScales[nscale];
        crntDyyImg = gyyScales[nscale];
        crntDxyImg = gxyScales[nscale];

        width  = crntDxxImg->width;
        height = crntDxxImg->height;
        crntHessImg = new Image(crntDxxImg->width, crntDxxImg->height);
        crntLoGImg  = new Image(crntDxxImg->width, crntDxxImg->height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dxx = crntDxxImg->getPixel(x, y);
                dyy = crntDyyImg->getPixel(x, y);
                dxy = crntDxyImg->getPixel(x, y);
                _sigma = HessLap::_SIGMA * pow(2.0,  nscale / (float) SCALES);
                det = _sigma*_sigma*(dxx*dyy - dxy*dxy);
                trace2 = dxx+dyy;

                lap = _sigma*trace2;
                crntHessImg->setPixel(x, y, det);
                crntLoGImg->setPixel(x,  y, lap);
            }
        }

        HessScales.push_back(crntHessImg);
        LoGScales.push_back(crntLoGImg);
    }

    return HessScales;
}

bool HessLap::isSpatialPeak(Image *image,const int x,const int y)
{
    float center =image->getPixel(x, y);
    int xi = 0, yi = 0;
    float p2 = 0.0f;

    if(center > 0.0)
    {
        for (yi = y - 1; yi <= y + 1; yi++)
        {
            for (xi = x - 1; xi <= x + 1; xi++)
            {
                p2 = image->getPixel(xi, yi);
                if (center < p2)
                    return false;
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool HessLap::adaptScale(int &x, int &y, vector<vector<Image *> > & GxxOctaves, vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves, int &scale, int &octave)
{
    scale = 1;
    octave  = 0;
    float xlamda = 0;
    int locx, locy;
    unsigned int is, ioctave;
    Image *crnt_Gyy, *crnt_Gxx, *crnt_Gxy;
    float a, b, c, e1, e2, ratio = 1.0f;
    locx = x;
    locy = y;

    for(ioctave = 0; ioctave < GyyOctaves.size(); ioctave++)
    {
        vector<Image*> &crnt_GyyScales = GyyOctaves[ioctave];
        vector<Image*> &crnt_GxxScales = GxxOctaves[ioctave];
        vector<Image*> &crnt_GxyScales = GxyOctaves[ioctave];
        for (is = 1; is < crnt_GyyScales.size()-1; is++)
        {
            crnt_Gxx = crnt_GxxScales[is];
            crnt_Gxy = crnt_GxyScales[is];
            crnt_Gyy = crnt_GyyScales[is];
            a = crnt_Gxx->getPixel(locx, locy);
            b = crnt_Gxy->getPixel(locx, locy);
            c = crnt_Gyy->getPixel(locx, locy);
            VMath::eigvlSymMat(a, b, c, e1, e2);
            if(e2 != 0)
            {
                ratio = e1/e2;
                if(ratio > xlamda)
                {
                    xlamda = ratio;
                    scale = is;
                    octave = ioctave;
                    x = locx;
                    y = locy;
                    cout<<x<<"\t"<<y<<"\t"<<scale<<"\t"<<octave<<"\t"<<xlamda<<endl;
                }
            }
        }
        locx = (int)round(locx/2);
        locy = (int)round(locy/2);
    }

    if(ratio < 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool HessLap::adaptAffine(vector<vector<Image *> > &HessOctaves, vector<vector<Image *> > &GxxOctaves,vector<vector<Image *> > &GyyOctaves, vector<vector<Image *> >& GxyOctaves)
{
    vector<KeyPoint*>::iterator it;
    KeyPoint *crntKp;
    Image *crnt_Gyy, *crnt_Gxx, *crnt_Gxy;
    int  is, x, y;
    float a, b, c, xn, yn;
    float u[2][2], ei[2], mi[2][2];
    unsigned int i = 0;
    int ioctave = 0;

    for (i = 0; i < this->kps.size(); i++)
    {
        crntKp = this->kps[i];
        xn = crntKp->sx;
        yn = crntKp->sy;
        x = (int)round(xn);
        y = (int)round(yn);

        ioctave = crntKp->octave;
        is = crntKp->scale;

        vector<Image*> &crnt_GyyScales = GyyOctaves[ioctave];
        vector<Image*> &crnt_GxxScales = GxxOctaves[ioctave];
        vector<Image*> &crnt_GxyScales = GxyOctaves[ioctave];
        crnt_Gxx  = crnt_GxxScales[is];
        crnt_Gxy  = crnt_GxyScales[is];
        crnt_Gyy  = crnt_GyyScales[is];

        a = crnt_Gxx->getPixel(x, y);
        b = crnt_Gxy->getPixel(x, y);
        c = crnt_Gyy->getPixel(x, y);

        crntKp->sx = xn;
        crntKp->sy = yn;
        crntKp->x = (int)round(xn * pow(2.0, ioctave));
        crntKp->y = (int)round(yn * pow(2.0, ioctave));

        VMath::inv_sqrtSymMat(a, b, c, u, ei, true);
        VMath::eigvtSymMat(u[0][0], u[0][1], u[1][1], mi, ei);
        crntKp->sori = atan2(mi[1][0], mi[0][0]);
        crntKp->e1 = ei[0];
        crntKp->e2 = ei[1];
        crntKp->a  = u[0][0];
        crntKp->b  = u[0][1];
        crntKp->c  = u[1][1];
    }//for_loop

    return true;
}

bool HessLap::isScalePeak(Image * aimage, Image * bimage,Image * cimage, float &logVal, const int x, const int y)
{
    assert(aimage);
    assert(bimage);
    assert(cimage);

    Image * ims[3];
    ims[0] = aimage;
    ims[1] = bimage;
    ims[2] = cimage;

    float center = bimage->getPixel(x, y);
    float p2;
    logVal       = center;

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
        for (int si = 0; si < 3; si++)
        {
            p2 = ims[si]->getPixel(x, y);
            if (center > p2)
                return false;
        }
        return true;
    }
}

vector<KeyPoint *> HessLap::FindPeaksScales(const int octave, vector<Image*> HessImages, vector<Image *> & LoGImages)
{
    vector<KeyPoint *> peaks;
    Image * kpfound = new Image(LoGImages[0]->width, LoGImages[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0, logVal = 0;
    int x = 0, y = 0;
    unsigned int s;

    float contrast_thresh = this->thresh;// / (float) NUM_SCALES;

    for (s = 1; s < HessImages.size()-1; s++)
    {

        for (y = BORDER; y < (HessImages[0]->height - BORDER); y++)
        {
            for (x = BORDER; x < (HessImages[0]->width - BORDER); x++)
            {
                funcVal = HessImages[s]->getPixel(x, y);
                if (funcVal <= contrast_thresh)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                if (!isSpatialPeak(HessImages[s], x, y))
                    continue;

                if(!isScalePeak(LoGImages[s-1],LoGImages[s],LoGImages[s+1], logVal, x , y))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                if (INTERP_KEYS)
                {
                    if (!InterpKey(x, y, s, HessImages, &fx, &fy, &fs,&funcVal))
                        continue;
                }

                fs  = fs < 0?0:fs;

                KeyPoint * peak = new KeyPoint();

                peak->fx = fx * pow(2.0, octave);
                peak->fy = fy * pow(2.0, octave);
                peak->x  = (int)round(peak->fx);
                peak->y  = (int)round(peak->fy);

                peak->dscale  = HessLap::_SIGMA * pow(2.0, octave + fs / (float) HessLap::SCALES);
                peak->iscale  = peak->dscale*HessLap::mag;
                peak->funcVal = funcVal;
                peak->ori     = 0;

                peak->scale   = s;
                peak->fscale  = fs;
                peak->gscale  = octave + fs/HessLap::SCALES;
                peak->sx      = fx;
                peak->sy      = fy;
                peak->octave  = octave;
                peak->KP      = true;
                peak->div     = logVal > 0?1:-1;
                peak->img_width = this->crntimg->width;

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

vector<KeyPoint*> HessLap::FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint*> newkps;
    KeyPoint *newkp;

    float sigma, m, theta, degree;
    int c, index;
    int indexa, indexb, indexc, j;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval,maxp;

    bool valid;
    float thetas[NumOrient] = {0};
    float weight;

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *HessLap::_SIGMA;
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
                    index  = ((int) (degree / DEGREE));
                    weight = m*gmat[y + c][x + c];
                    index  = index%NumOrient;
                    thetas[index] +=  weight;

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

                newkps.push_back(newkp);

            }
        }
    }
    return  newkps ;
}

bool HessLap::KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char *dvfn)
{
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *oriImg = NULL, *tmpImg = NULL;

    vector<vector<Image*> > GxxOctaves, GyyOctaves, GxyOctaves;
    vector<vector<Image *> > LoGOctaves;
    vector<vector<Image *> > HessOctaves;
    vector<vector<Image *> > GaussOctaves;
    vector<KeyPoint *>  extra_peaks;

    unsigned int ioctave = 0;

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
    tmpImg = Filter::ScaleImage(this->crntimg, 0, HessLap::_SIGMA, INITSIGMA);
    this->crntimg = tmpImg;

    this->intImg = AbstractDetector::buildIntImage(this->crntimg);

    this->BuildOctaves(this->crntimg, _Dxx_, GxxOctaves, HessLap::_SIGMA, HessLap::SCALES);
    this->BuildOctaves(this->crntimg, _Dyy_, GyyOctaves, HessLap::_SIGMA, HessLap::SCALES);
    this->BuildOctaves(this->crntimg, _Dxy_, GxyOctaves, HessLap::_SIGMA, HessLap::SCALES);
    this->BuildOctaves(this->crntimg, _Conv2,GaussOctaves, HessLap::_SIGMA, HessLap::SCALES);
    HessOctaves = this->BuildHessOctaves(GxxOctaves, GyyOctaves, GxyOctaves, LoGOctaves);

    for(ioctave = 0; ioctave < LoGOctaves.size(); ioctave++)
    {
        vector<KeyPoint *> peaks = this->FindPeaksScales(ioctave,HessOctaves[ioctave],LoGOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    if(this->AFF_OUT)
    {
        this->adaptAffine(HessOctaves, GxxOctaves, GyyOctaves, GxyOctaves);
    }

    Cleaner::releaseOctaves(GxxOctaves);
    Cleaner::releaseOctaves(GyyOctaves);
    Cleaner::releaseOctaves(GxyOctaves);
    Cleaner::releaseOctaves(LoGOctaves);
    Cleaner::releaseOctaves(HessOctaves);

    extra_peaks = this->FindOrientByGrad(kps, GaussOctaves);

    Cleaner::releaseOctaves(GaussOctaves);

    if(extra_peaks.size() > 0)
    {
        this->kps.insert(kps.begin(),extra_peaks.begin(),extra_peaks.end());
        extra_peaks.clear();
    }

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

    if(strcmp(dstfn,""))
    {
        (this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
    }

    delete this->crntimg;
    this->crntimg = oriImg;

    /**/
    if(strcmp(descfn,"") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }
    /**/

    if(strcmp(dvfn, "")&&this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }

    delete this->intImg;
    delete this->crntimg;
    this->crntimg = NULL;
    this->intImg  = NULL;

    return true;
}

void HessLap::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new HessLap();
    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 5; i++)
    {
        sprintf(vsrc, "/home/wlzhao/datasets/vgg/graf/img%d.jpg",    i);
        sprintf(vdst, "/home/wlzhao/datasets/vgg/graf/img%d.keys",   i);
        sprintf(vdesc, "/home/wlzhao/datasets/vgg/graf/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }

    //mydetector->printAVG();
}
