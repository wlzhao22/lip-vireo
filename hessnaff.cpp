#include "hessnaff.h"

#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

/**constant initialization**/

const int HessnAff::MaxOctaves    = 4;
const int HessnAff::SCALES        = 6;
const float HessnAff::_SIGMA      = 1.4;
const float HessnAff::INITSIGMA   = 0.5;
const bool  HessnAff::INTERP_KEYS = true;

const float HessnAff::mag         = 5.1962; //optimal 3.0 // 3*sqrt(3)
const int HessnAff::BORDER        = 5;
const int HessnAff::THRESH        = 600; //100

const int HessnAff::DEGREE        = 10;
const int HessnAff::NOrient       = 36;
const int HessnAff::DEGPERBIN     = (360 / NOrient);
const float HessnAff::NwKpThresh  = 0.8;

const float HessnAff::cvtRatio    = 10.0f;

const int HessnAff::maxIter       = 10;
const int HessnAff::W0            = 10;
const float HessnAff::err_c       = 0.80f;

/**
* I try to use a dynamic threshold, however, no significant improvement I can get,
it actually produces a few bigger number of interest points

@date Mar.2nd.2009

const float HessnAff::thresh_ratio = 2.5f;
**/

/**
I modified the implementation completely. Now it performs as good as HessnAff-Laplacian.
However, it is really HessnAff detector

@date 11-May-2011
@author Wan-Lei Zhao
**/

/**
I modified the implementation a little bit, it performs better and faster

@date 11-May-2015
@author Wan-Lei Zhao


Shuangquan Feng and I worked out a new way to estimate the affine region. The affine region looks reasonably good,
however it does not perform well. So I build Hessian-Affine upon my Hessian detector
**/


HessnAff::HessnAff()
{
    cout<<"Detector ................................ HessnAff\n";
    this->DETECTOR = hesaff;
    this->intImg   = NULL;
    this->sel_option = THRSH;
    this->AFF_OUT    = false;
}

bool HessnAff::paramsCheck()
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
        this->sel_option = TOPK;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    if(this->paras.find(argv[3]) != this->paras.end())
    {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option = DENS;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    return true;
}

float HessnAff::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
{

    float Dp[3] = {0}; // first derivative of D with respect to x, y, s

    float Dpp[3][3] = {{0}}; // HessnAff of D

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

bool HessnAff::InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal)
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

vector<vector<Image *> > HessnAff::BuildHessOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves)
{
    vector<Image *> gxxScales, gyyScales, gxyScales;
    vector<vector<Image *> > HessOctaves;

    for (unsigned int i = 0; i < GxxOctaves.size(); i++)
    {
        gxxScales = GxxOctaves[i];
        gyyScales = GyyOctaves[i];
        gxyScales = GxyOctaves[i];
        vector<Image *> HessScales;
        HessScales = BuildHessScales(gxxScales, gyyScales, gxyScales);
        HessOctaves.push_back(HessScales);
    }

    return HessOctaves;
}

vector<Image* > HessnAff::BuildHessScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales,vector<Image *> &gxyScales)
{
    vector<Image* > HessScales;
    Image *crntDxxImg = NULL, *crntDyyImg = NULL, *crntDxyImg = NULL, *crntHessImg = NULL;
    float dxx = 0.0f, dyy = 0.0f, dxy = 0.0f;
    int x = 0, y = 0, width = 0, height = 0;
    float _sigma = 0, det = 0;
    unsigned int nscale = 0;

    for(nscale = 0; nscale < gxxScales.size(); nscale++)
    {
        crntDxxImg = gxxScales[nscale];
        crntDyyImg = gyyScales[nscale];
        crntDxyImg = gxyScales[nscale];

        width   = crntDxxImg->width;
        height  = crntDxxImg->height;
        crntHessImg = new Image(crntDxxImg->width, crntDxxImg->height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dxx = crntDxxImg->getPixel(x, y);
                dyy = crntDyyImg->getPixel(x, y);
                dxy = crntDxyImg->getPixel(x, y);
                _sigma = HessnAff::_SIGMA * pow(2.0,  nscale / (float) SCALES);
                det = _sigma*_sigma*(dxx*dyy - dxy*dxy);

                crntHessImg->setPixel(x, y, det);
            }
        }
        HessScales.push_back(crntHessImg);
    }

    return HessScales;
}

bool HessnAff::isSpatialPeak(Image *image,const int x,const int y)
{
    float center =image->getPixel(x, y);
    float p2 = 0.0f;
    int xi = 0, yi = 0;

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
        for (yi = y - 1; yi <= y + 1; yi++)
        {
            for (xi = x - 1; xi <= x + 1; xi++)
            {
                p2 = image->getPixel(xi, yi);
                if (center > p2)
                    return false;
            }
        }
        return true;
    }
}

bool HessnAff::isScalePeak(Image * aimage, Image * bimage,Image * cimage, const int x, const int y)
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

    if (center > 0.0)
    {
        for (int si = 0; si < 3; si++)
        {
            p2 = ims[si]->getPixel(x, y);
            if (center < p2)
                return false;
        }
    }
    else
    {
        for (int si = 0; si < 3; si++)
        {
            p2 = ims[si]->getPixel(x, y);
            if (center > p2)
                return false;
        }
        return false;
    }
    return true;
}

vector<KeyPoint *> HessnAff::FindPeaksScales(const int octave, vector<Image*> &HessImages, vector<Image*> &GxxImages, vector<Image*> &GyyImages)
{
    vector<KeyPoint *> peaks;
    Image * kpfound = new Image(HessImages[0]->width, HessImages[0]->height);
    float fx = 0.0f, fy = 0.0f, fs = 0.0f, funcVal = 0.0f, trace = 0.0f;
    int x = 0, y = 0;
    unsigned int s = 0;
    float contrast_thresh = this->thresh;

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

                if(!isScalePeak(HessImages[s-1], HessImages[s], HessImages[s+1], x, y))
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
                trace = GxxImages[s]->getPixel(x, y) + GyyImages[s]->getPixel(x, y);

                KeyPoint * peak = new KeyPoint();

                peak->x       = (int)round(fx * pow(2.0, octave));
                peak->y       = (int)round(fy * pow(2.0, octave));
                peak->dscale  = HessnAff::_SIGMA * pow(2.0, octave + fs / (float) SCALES);
                peak->octSigma = HessnAff::_SIGMA * pow(2.0, fs / (float) SCALES);
                peak->iscale  = peak->dscale*HessnAff::mag;
                peak->funcVal = funcVal;
                peak->ori     = 0;
                peak->scale   = s;
                peak->fscale  = fs;
                peak->gscale  = octave + fs/HessnAff::SCALES;
                peak->sx      = x;
                peak->sy      = y;
                peak->octave  = octave;
                peak->octIndex = octave;
                peak->KP      = true;
                peak->div     = trace > 0 ? 1:-1;
                peak->img_width = this->crntimg->width;
                peak->fx      = x * pow(2.0, octave);
                peak->fy      = y * pow(2.0, octave);

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

bool HessnAff::getdRds(const float block[5][5], const float sigma, float &dxds, float &dyds)
{
    float dxx, dyy, dxy, xn, xp, xy, yn, yp, xpyp, xnyn, xnyp, xpyn;
    float det = 0, v = 0, xnn, xpp, ynn, ypp, dxxx, dxxy, dxyy, dyyy;

    xy = block[2][2];
    xn = block[2][1];    xp = block[2][3];
    yn = block[1][2];    yp = block[3][2];
    xpyp = block[3][3];  xnyn = block[1][1];
    xnyp = block[3][1];  xpyn = block[1][3];

    dxx = xn + xp - 2 * xy;
    dyy = yn + yp - 2 * xy;
    dxy = (xpyp + xnyn - xnyp - xpyn) / 4;
    det = dxx * dyy - dxy * dxy;
    v   = -sigma / det;

    xnn = block[2][0];    xpp = block[2][4];
    ynn = block[0][2];    ypp = block[4][2];

    dxxx = (xpp - 2 * xp + 2 * xn - xnn) / 2;
    dyyy = (ypp - 2 * yp + 2 * yn - ynn) / 2;
    dxyy = (2 * (xn - xp) + xpyn + xpyp - xnyn - xnyp) / 2;
    dxxy = (2 * (yn - yp) + xnyp + xpyp - xnyn - xpyn) / 2;

    dxds = v * (dyy * dxxx + dyy * dxyy - dxy * dxxy - dxy * dyyy);
    dyds = v * (dxx * dyyy + dxx * dxxy - dxy * dxxx - dxy * dxyy);

    return true;
}

float HessnAff::getdHds(const float block[5][5], float sigma)
{
    float t = sigma;
    float dxx, dyy, dxy, dxxxx, dyyyy, dxxyy, dxxxy, dxyyy, dhs = 0;
    float xn = block[2][1]; //window->getPixel(cx - 1, cy); // x negative 1
    float xp = block[2][3];//window->getPixel(cx + 1, cy); // x positive 1
    float xy = block[2][2];//window->getPixel(cx, cy);
    float yn = block[1][2];//window->getPixel(cx, cy - 1);
    float yp = block[3][2];//window->getPixel(cx, cy + 1);
    float xpyp = block[3][3];//window->getPixel(cx + 1, cy + 1);
    float xnyn = block[1][1];//window->getPixel(cx - 1, cy - 1);
    float xnyp = block[3][1];//window->getPixel(cx - 1, cy + 1);
    float xpyn = block[1][3];//window->getPixel(cx + 1, cy - 1);

    dxx   = xn + xp - 2 * xy;
    dyy   = yn + yp - 2 * xy;
    dxy   = (xpyp + xnyn - xnyp - xpyn) / 4;
    dxxxx = block[2][0]  + block[2][4] - 4 * (xn + xp) + 6 * xy;
    dyyyy = block[0][2]  + block[4][2] - 4 * (yn + yp) + 6 * xy;
    dxxyy = xpyp + xpyn  + xnyp + xnyn - 2 * (xn + xp + yn + yp) + 4 * xy;
    dxxxy = (block[1][0] + block[3][4] - block[3][0] - block[1][4])/4 + (xnyp + xpyn - xnyn - xpyp)/2;
    dxyyy = (block[0][1] + block[4][3] - block[0][3] - block[4][1])/4 + (xnyp + xpyn - xnyn - xpyp)/2;
    dhs   = ((dxxxx + dxxyy) * dyy + (dxxyy + dyyyy) * dxx)/2 - dxy*(dxxxy + dxyyy);
    return 2 * t * (dxx * dyy - dxy * dxy) + t * t * dhs;
}

bool HessnAff::Inverse2D(float m[2][2], float invm[2][2])
{
       float det = det2D(m);
       if (fabs(det) > 1e-8)
       {
            float divDet = 1 / det;
            invm[0][0] = divDet * m[1][1];
            invm[0][1] = -divDet * m[0][1];
            invm[1][0] = -divDet * m[1][0];
            invm[1][1] = divDet * m[0][0];
            return true;
        } else
            return false;
}


void HessnAff::affineAdapt(vector<KeyPoint *> &peaks, vector<vector<Image*> > &gaussPyramid, vector<vector<Image*> > &dxPyramid, vector<vector<Image*> > &dyPyramid)
{
    int times = 0, num = 0, nPxs = 0, cx = 0, cy = 0;
    int oct = 0, si = 0, xi = 0, yi = 0, ri = 0, x, y, wd = 0, ci = 0;
    unsigned int ioct;

    vector<KeyPoint*>::iterator it;
    vector<Cords>::iterator cit;
    vector<Cords> visited;
    vector<Board *> octFlags;
    Board *bdFlag = NULL;

    Image *dxImg = NULL, *dyImg = NULL, *gsImg = NULL;
    KeyPoint *crntKp = NULL;

    float fx = 0, fy = 0, dx = 0, dy = 0, px = 0, sx, sy = 0, dr = 0, r = 4.81;
    float itSigma = 0, ds = 0, tmpdt = 0, dt = 0, dxds, dyds, dta, theta;
    float U[2][2] = {0}, ivU[2][2] = {0}, nwU[2][2]={0}, tU[2][2];
    float M[2][2] = {0}, lm1 = 0, lm2 = 0, dlta = 0, sc = 0, tr = 0;
    float tmpU[2][2] = {0}, es[2] = {0}, block[5][5];
    bool _CNVG_ = false;

    for(ioct = 0; ioct < gaussPyramid.size(); ioct++)
    {
        gsImg  = gaussPyramid[ioct][0];
        bdFlag = new Board(gsImg->width, gsImg->height);
        octFlags.push_back(bdFlag);
    }
    ///cout<<"\t";
    for(it = peaks.begin(); it != peaks.end(); it++, num++)
    {
        crntKp  = *it;
        if(crntKp->KP == false)
        continue;

        oct     = crntKp->octIndex;
        si      = crntKp->scale;
        itSigma = crntKp->octSigma;
        U[0][0] = U[1][1] = 1;
        U[0][1] = U[1][0] = 0;
        _CNVG_  = false;
        r       = 4.5;
        times   = 0;
        bdFlag  = octFlags[oct];
        es[0]   = es[1] = 1;
        ds      = 0.08;
        dt      = 1.0;
        sx      = crntKp->sx; sy = crntKp->sy;
        dxImg   = dxPyramid[oct][si];
        dyImg   = dyPyramid[oct][si];
        gsImg   = gaussPyramid[oct][si];
        while(!_CNVG_ && times < 10)
        {
            M[0][0] = M[0][1] = M[1][0] = M[1][1] = 0;
            sc = (itSigma*3.0)/r;
            wd = (int)floor((sc*es[1]*es[0]*itSigma*3.0)/r);
            wd = wd*2+1;
            nPxs  = 0;
            visited.clear();
            for(theta = 0; theta < 6.32; theta += 0.45)
            {
                for(dr = 0; dr < 4.81; dr += 0.8)
                {
                   xi = dr*cos(theta);
                   yi = dr*sin(theta);
                   fx = sc*xi*U[0][0] + sc*yi*U[0][1];
                   fy = sc*xi*U[1][0] + sc*yi*U[1][1];
                   x  = (int)floor(sx + fx);
                   y  = (int)floor(sy + fy);
                   for(ri = -wd; ri <= wd; ri++)
                   {
                       cy = y + ri;
                       for(ci = -wd; ci <= wd; ci++)
                       {
                          cx = x + ci;
                          if(bdFlag->getPixel(cx, cy) == 0)
                          {
                             dx = dxImg->getPixel(cx, cy);
                             dy = dyImg->getPixel(cx, cy);
                             M[0][0] += dx*dx;
                             M[0][1] += dx*dy;
                             M[1][1] += dy*dy;
                             nPxs++;
                             bdFlag->setPixel(cx, cy, 1);
                             Cords pt; pt.x = cx; pt.y = cy;
                             visited.push_back(pt);
                          }
                       }
                    }
                }
            }

            for(cit = visited.begin(); cit != visited.end(); cit++)
            {
                  Cords & pt = *cit;
                  bdFlag->setPixel(pt.x, pt.y, 0);
            }
            visited.clear();
            for(yi = 0; yi < 5; yi++)
            {
               for(xi = 0; xi < 5; xi++)
               {
                    x  = xi - 2; y = yi - 2;
                    fx = sc*x*U[0][0] + sc*y*U[0][1];
                    fy = sc*x*U[1][0] + sc*y*U[1][1];
                    x  = crntKp->sx + fx;
                    y  = crntKp->sy + fy;
                    px = gsImg->getPixelBI(x, y);
                    block[yi][xi] = px;
               }
            }

            tmpdt = getdHds(block, itSigma);
            if(tmpdt*dt < 0)
            {
                ds = ds/2;
            }
            if(tmpdt > 0)
            {
               itSigma = itSigma + ds;
            }else{
               itSigma = itSigma - ds;
            }
            dt = tmpdt;
            this->getdRds(block, itSigma, dxds, dyds);
            sx = sx + ds * dxds;
            sy = sy + ds * dyds;
            ///cout<<"tm: "<<times<<"\tds: "<<ds<<"\tsig: "<<itSigma<<"\txy: "<<dxds<<"\t"<<dyds<<endl;

            M[1][0] = M[0][1];
            M[0][0] = M[0][0]/nPxs; M[0][1] = M[0][1]/nPxs;
            M[1][0] = M[1][0]/nPxs; M[1][1] = M[1][1]/nPxs;

            tmpU[0][0] = M[0][0]; tmpU[0][1] = M[0][1];
            tmpU[1][0] = M[1][0]; tmpU[1][1] = M[1][1];
            ///cout<<"M: "<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<nPxs<<endl;

            VMath::sqrtSymMat(tmpU[0][0], tmpU[0][1], tmpU[1][1], nwU, es);
            Inverse2D(U, ivU);
            tU[0][0] = nwU[0][0]*ivU[0][0] + nwU[1][0]*ivU[1][0];
            tU[0][1] = nwU[0][0]*ivU[0][1] + nwU[1][0]*ivU[1][1];
            tU[1][0] = nwU[1][0]*ivU[0][0] + nwU[1][1]*ivU[1][0];
            tU[1][1] = nwU[1][0]*ivU[0][1] + nwU[1][1]*ivU[1][1];
            tr   = tU[0][0] + tU[1][1];
            dlta = tU[0][0]*tU[1][1] - tU[1][0]*tU[0][1];
            lm1  = tr/2 + sqrt(tr*tr/4-dlta);
            lm2  = tr/2 - sqrt(tr*tr/4-dlta);
            ///cout<<"lm1: "<<lm1<<"\t"<<lm2<<endl;
            if(lm2/lm1 > 0.96)
            {
                _CNVG_ = true;
            }
            U[0][0] = nwU[0][0]; U[0][1] = nwU[0][1];
            U[1][0] = nwU[1][0]; U[1][1] = nwU[1][1];
            ///exit(0);
            times++;
        }///while(!_CNVG_)
        ///exit(0);
        if(times >= 10)
        {
           crntKp->KP = false;
        }
        ///cout<<num<<"\t"<<times<<"\t"<<lm2/lm1<<"\tdt = "<<dt<<endl;
        crntKp->a = U[0][0]; crntKp->b = U[0][1]; crntKp->c = U[1][1];
        ///cout<<crntKp->a<<"\t"<<crntKp->b<<"\t"<<crntKp->c<<endl;
        crntKp->octSigma = itSigma;
        crntKp->iscale   =  pow(2.0, crntKp->octIndex)*itSigma*HessnAff::mag;
        crntKp->dscale   =  pow(2.0, crntKp->octIndex)*itSigma;
        crntKp->x = (int) round(sx * pow(2.0, crntKp->octIndex));
        crntKp->y = (int) round(sy * pow(2.0, crntKp->octIndex));
        ///cout<<"\r\r\r\r\t"<<num<<"/"<<peaks.size();
    }
    Cleaner::releaseBoards(octFlags);
    ///cout<<"done\n";
    ///cout<<endl;
}

void HessnAff::affineAdapt2D(vector<KeyPoint *> &peaks, vector<vector<Image*> > &gaussPyramid, vector<vector<Image*> > &dxxPyramid,
                           vector<vector<Image*> > &dyyPyramid, vector<vector<Image*> > &dxyPyramid)
{
    int times = 0, num = 0, nPxs = 0, cx = 0, cy = 0;
    int oct = 0, si = 0, xi = 0, yi = 0, ri = 0, x, y, wd = 0, ci = 0;
    unsigned int ioct;

    vector<KeyPoint*>::iterator it;
    vector<Cords>::iterator cit;
    vector<Cords> visited;
    vector<Board *> octFlags;
    Board *bdFlag = NULL;

    Image *dxxImg = NULL, *dyyImg = NULL, *dxyImg = NULL, *gsImg = NULL;
    KeyPoint *crntKp = NULL;

    float fx = 0, fy = 0, dxx = 0, dyy = 0, dxy = 0, px = 0, sx, sy = 0, dr = 0, r = 4.81;
    float itSigma = 0, ds = 0, tmpdt = 0, dt = 0, dxds, dyds, dta, theta;
    float U[2][2] = {0}, ivU[2][2] = {0}, nwU[2][2]={0}, tU[2][2];
    float M[2][2] = {0}, lm1 = 0, lm2 = 0, dlta = 0, sc = 0, tr = 0;
    float tmpU[2][2] = {0}, es[2] = {0}, block[5][5];
    bool _CNVG_ = false;

    for(ioct = 0; ioct < gaussPyramid.size(); ioct++)
    {
        gsImg  = gaussPyramid[ioct][0];
        bdFlag = new Board(gsImg->width, gsImg->height);
        octFlags.push_back(bdFlag);
    }
    ///cout<<"\t";
    for(it = peaks.begin(); it != peaks.end(); it++, num++)
    {
        crntKp  = *it;
        if(crntKp->KP == false)
        continue;

        oct     = crntKp->octIndex;
        si      = crntKp->scale;
        itSigma = crntKp->octSigma;
        U[0][0] = U[1][1] = 1;
        U[0][1] = U[1][0] = 0;
        _CNVG_  = false;
        r       = 4.5;
        times   = 0;
        bdFlag  = octFlags[oct];
        es[0]   = es[1] = 1;
        ds      = 0.08;
        dt      = 1.0;
        sx      = crntKp->sx; sy = crntKp->sy;
        dxxImg  = dxxPyramid[oct][si];
        dyyImg  = dyyPyramid[oct][si];
        dxyImg  = dxyPyramid[oct][si];
        gsImg   = gaussPyramid[oct][si];

        while(!_CNVG_ && times < 10)
        {
            M[0][0] = M[0][1] = M[1][0] = M[1][1] = 0;
            sc = (itSigma*3.0)/r;
            wd = (int)floor((sc*es[1]*es[0]*itSigma*3.0)/r);
            wd = wd*2+1;
            nPxs  = 0;
            visited.clear();

            for(theta = 0; theta < 6.32; theta += 0.45)
            {
                for(dr = 0; dr < 4.81; dr += 0.8)
                {
                   xi = dr*cos(theta);
                   yi = dr*sin(theta);

                   fx = sc*xi*U[0][0] + sc*yi*U[0][1];
                   fy = sc*xi*U[1][0] + sc*yi*U[1][1];
                   x  = (int)floor(sx + fx);
                   y  = (int)floor(sy + fy);
                   for(ri = -wd; ri <= wd; ri++)
                   {
                       cy = y + ri;
                       for(ci = -wd; ci <= wd; ci++)
                       {
                          cx = x + ci;
                          if(bdFlag->getPixel(cx, cy) == 0)
                          {
                             dxx = dxxImg->getPixel(cx, cy);
                             dyy = dyyImg->getPixel(cx, cy);
                             dxy = dxyImg->getPixel(cx, cy);
                             M[0][0] += dxx*dxx;
                             M[0][1] += dxy*dxy;
                             M[1][1] += dyy*dyy;
                             nPxs++;
                             bdFlag->setPixel(cx, cy, 1);
                             Cords pt; pt.x = cx; pt.y = cy;
                             visited.push_back(pt);
                          }//else
                       }
                       //cout<<endl;
                    }
                }
            }

            for(cit = visited.begin(); cit != visited.end(); cit++)
            {
                  Cords & pt = *cit;
                  bdFlag->setPixel(pt.x, pt.y, 0);
            }
            visited.clear();
            for(yi = 0; yi < 5; yi++)
            {
               for(xi = 0; xi < 5; xi++)
               {
                    x  = xi - 2; y = yi - 2;
                    fx = sc*x*U[0][0] + sc*y*U[0][1];
                    fy = sc*x*U[1][0] + sc*y*U[1][1];
                    x  = crntKp->sx + fx;
                    y  = crntKp->sy + fy;
                    px = gsImg->getPixelBI(x, y);
                    block[yi][xi] = px;
               }
            }

            tmpdt = getdHds(block, itSigma);
            if(tmpdt*dt < 0)
            {
                ds = ds/2;
            }
            if(tmpdt > 0)
            {
               itSigma = itSigma + ds;
            }else{
               itSigma = itSigma - ds;
            }
            dt = tmpdt;
            this->getdRds(block, itSigma, dxds, dyds);
            sx = sx + ds * dxds;
            sy = sy + ds * dyds;

            M[1][0] = M[0][1];
            M[0][0] = M[0][0]/nPxs; M[0][1] = M[0][1]/nPxs;
            M[1][0] = M[1][0]/nPxs; M[1][1] = M[1][1]/nPxs;

            tmpU[0][0] = M[0][0]; tmpU[0][1] = M[0][1];
            tmpU[1][0] = M[1][0]; tmpU[1][1] = M[1][1];

            VMath::sqrtSymMat(tmpU[0][0], tmpU[0][1], tmpU[1][1], nwU, es);
            Inverse2D(U, ivU);
            tU[0][0] = nwU[0][0]*ivU[0][0] + nwU[1][0]*ivU[1][0];
            tU[0][1] = nwU[0][0]*ivU[0][1] + nwU[1][0]*ivU[1][1];
            tU[1][0] = nwU[1][0]*ivU[0][0] + nwU[1][1]*ivU[1][0];
            tU[1][1] = nwU[1][0]*ivU[0][1] + nwU[1][1]*ivU[1][1];
            tr   = tU[0][0] + tU[1][1];
            dlta = tU[0][0]*tU[1][1] - tU[1][0]*tU[0][1];
            lm1  = tr/2 + sqrt(tr*tr/4-dlta);
            lm2  = tr/2 - sqrt(tr*tr/4-dlta);
            if(lm2/lm1 > 0.96)
            {
                _CNVG_ = true;
            }
            U[0][0] = nwU[0][0]; U[0][1] = nwU[0][1];
            U[1][0] = nwU[1][0]; U[1][1] = nwU[1][1];
            times++;
        }///while(!_CNVG_)
        if(times >= 10)
        {
           crntKp->KP = false;
        }
        crntKp->a = U[0][0]; crntKp->b = U[0][1]; crntKp->c = U[1][1];
        crntKp->octSigma = itSigma;
        crntKp->iscale   =  pow(2.0, crntKp->octIndex)*itSigma*HessnAff::mag;
        crntKp->dscale   =  pow(2.0, crntKp->octIndex)*itSigma;
        crntKp->x = (int)round(sx * pow(2.0, crntKp->octIndex));
        crntKp->y = (int)round(sy * pow(2.0, crntKp->octIndex));
        crntKp->fx = sx * pow(2.0, crntKp->octIndex);
        crntKp->fy = sy * pow(2.0, crntKp->octIndex);
    }
    Cleaner::releaseBoards(octFlags);
}


vector<KeyPoint*> HessnAff::findOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GaussianOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint*> newkps;
    KeyPoint *newkp;

    float sigma, m, theta, degree, fx, fy;
    float U[2][2];
    int c, index;
    int indexa, indexb, indexc, j;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval,maxp;

    bool valid;
    float thetas[NOrient] = {0};
    float weight;

    for (i = 0; i < kps.size(); i++)
    {
        if(kps[i]->KP == false)
        continue;

        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *HessnAff::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size()/2;

        for(j = 0; j < NOrient; j++)
        {
            thetas[j] = 0;
        }
        U[0][0] = kps[i]->a; U[0][1] = kps[i]->b;
        U[1][0] = kps[i]->b; U[1][1] = kps[i]->c;

        for (int y = -c; y <= c; y++)
        {
            for (int x = -c; x <= c; x++)
            {
                if (sqrt((float) x*x + y*y) > 3.0 * sigma)
                    continue;

                fx = x*U[0][0] + y*U[0][1];
                fy = x*U[1][0] - y*U[1][1];

                valid = Filter::GetPixOrientation((int) (kps[i]->sx + fx + 0.5), (int) (kps[i]->sy + fy + 0.5),
                                                  GaussianOctaves[kps[i]->octave][kps[i]->scale], m, theta);

                if(valid)
                {
                    degree = theta / PI * 180.0 + 180.0;
                    index  = ((int) (degree / DEGREE));
                    weight = m*gmat[y + c][x + c];
                    index  = index%NOrient;
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
            VMath::SmoothHistogram(thetas, NOrient);

        maxval = VMath::maxVec(thetas, NOrient);

        for (int j = 0; j < NOrient; j++)
        {
            if (thetas[j] < maxval * NwKpThresh)
                continue;

            indexa = VMath::mod(j - 1, NOrient);
            indexb = j;
            indexc = VMath::mod(j + 1, NOrient);
            thetaa = thetas[indexa];
            thetab = thetas[indexb];
            thetac = thetas[indexc];

            if (!(thetab > thetaa && thetab > thetac))
                continue;

            maxp = VMath::getMaxParabola(-1, thetaa, 0, thetab, 1, thetac);

            if(thetas[j] == maxval)
            {
                kps[i]->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;
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
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;

                newkps.push_back(newkp);

            }
        }
    }
    return  newkps ;
}

bool HessnAff::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *oriImg = NULL, *tmpImg = NULL;

    vector<vector<Image*> >  GxxOctaves, GyyOctaves, GxyOctaves;
    vector<vector<Image*> >  GxOctaves, GyOctaves;
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
    tmpImg = Filter::ScaleImage(this->crntimg, 0, _SIGMA, INITSIGMA);
    this->crntimg = tmpImg;
    this->intImg = this->buildIntImage(this->crntimg);

    ///this->BuildOctaves(this->crntimg, _Dx_,  GxOctaves,  HessnAff::_SIGMA, HessnAff::SCALES);
    ///this->BuildOctaves(this->crntimg, _Dy_,  GyOctaves,  HessnAff::_SIGMA, HessnAff::SCALES);
    this->BuildOctaves(this->crntimg, _Dxx_, GxxOctaves, HessnAff::_SIGMA, HessnAff::SCALES);
    this->BuildOctaves(this->crntimg, _Dyy_, GyyOctaves, HessnAff::_SIGMA, HessnAff::SCALES);
    this->BuildOctaves(this->crntimg, _Dxy_, GxyOctaves, HessnAff::_SIGMA, HessnAff::SCALES);
    this->BuildOctaves(this->crntimg, _Conv2, GaussOctaves, HessnAff::_SIGMA, HessnAff::SCALES);
    HessOctaves = this->BuildHessOctaves(GxxOctaves, GyyOctaves, GxyOctaves);

    for(ioctave = 0; ioctave < HessOctaves.size(); ioctave++)
    {
        vector<KeyPoint *> peaks = this->FindPeaksScales(ioctave, HessOctaves[ioctave],
                                                         GxxOctaves[ioctave], GyyOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    Cleaner::releaseOctaves(HessOctaves);

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
        //AbstractDetector::topkSelect(kps,this->fix_kp_numb);
        break;
    }
    }

    affineAdapt2D(kps, GaussOctaves, GxxOctaves, GyyOctaves, GxyOctaves);
    ///affineAdapt(kps, GaussOctaves, GxOctaves, GyOctaves);

    extra_peaks = this->findOrientByGrad(kps, GaussOctaves);

    Cleaner::releaseOctaves(GxxOctaves);
    Cleaner::releaseOctaves(GyyOctaves);
    Cleaner::releaseOctaves(GxyOctaves);
    //Cleaner::releaseOctaves(GxOctaves);
    //Cleaner::releaseOctaves(GyOctaves);
    Cleaner::releaseOctaves(GaussOctaves);

    if(strcmp(dstfn, ""))
    {
        (this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
    }

    delete this->crntimg;
    this->crntimg = oriImg;

    if(strcmp(descfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn,this->resize_rate);
    }

    if(strcmp(dvfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }

    delete this->intImg;
    delete this->crntimg;
    this->intImg   = NULL;
    this->crntimg  = NULL;

    return true;
}

void HessnAff::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new HessnAff();
    mydetector->Init(config, "SIFT");
    /**/
    for(int i = 1; i <= 4; i++)
    {
        sprintf(vsrc, "/home/wlzhao/datasets/trec03/img%d.jpg",    i);
        sprintf(vdst, "/home/wlzhao/datasets/trec03/img%d.keys",   i);
        sprintf(vdesc, "/home/wlzhao/datasets/trec03/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }
    /**/
}

HessnAff::~HessnAff()
{
}
