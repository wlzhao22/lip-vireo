#include "dsurf.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "keypoint.h"
#include "descsift.h"
#include "intimage.h"
#include "kpdrawer.h"
#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

using namespace std;

const int DSURF::MaxOctaves  = 3;
const int DSURF::PatchMag    = 20;
const float DSURF::INITSIGMA = 0.5;
const float DSURF::SIGMA     = 1.2;
const int DSURF::SCALES      = 4;

const float DSURF::w0        = 0.81;
const float DSURF::THRESH    = 6.0;
const int DSURF::INTERP_KEYS = 1;
const float DSURF::LWTHRESH  = 0.035;

const int DSURF::PatchSize   = 41;
const float DSURF::mag       = 8; // 5
const int DSURF::kpnum0      = 750;
const float DSURF::PI0       = 3.14159265359f;

DSURF::DSURF()
{
    cout<<"Detector ................................ DSURF\n";
    this->featsBin = new float[64];
    this->descWin  = new float[DSURF::PatchSize*DSURF::PatchSize];
    this->DETECTOR = dsurf;
    this->AFF_OUT  = false;
    this->intImg   = NULL;
    this->sel_option = THRSH;

    Octs[0][0] = 9;  Octs[0][1] = 15;  Octs[0][2] = 21;  Octs[0][3] = 27;
    Octs[1][0] = 15; Octs[1][1] = 27;  Octs[1][2] = 39;  Octs[1][3] = 51;
    Octs[2][0] = 27; Octs[2][1] = 51;  Octs[2][2] = 75;  Octs[2][3] = 99;
}

bool DSURF::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    //if(paras.find(argv[0]) != paras.end())
    //    this->thresh = atof(paras["sigma"]);

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
        this->sel_option = TOPK;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    if(paras.find(argv[3]) != paras.end())
    {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option = DENS;
        cout<<"Topk .................................... "<<this->fix_kp_numb<<endl;
    }

    return true;
}

float DSURF::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
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

bool DSURF::InterpKey(int x, int y, int s, vector<Image *> & DImages, float * fx, float * fy, float * fs,float *dogVal)
{
    bool addkey = true;
    int moves_left = 5;
    int tx = x;
    int ty = y;
    int ts = s;

    float dx, dy, ds, val;
    bool updated;

    do
    {
        moves_left--;
        updated = false;

        val = InterpKeyStep(tx, ty, ts, DImages, &dx, &dy, &ds);

        if (val < LWTHRESH / (float) SCALESPEROCTAVE)
        {
            addkey = false;
            continue;
        }

        if (dx > 0.6 && tx < DImages[0]->width - 3)
        {
            tx++;
            updated = true;
        }
        else if (dx < -0.6 && tx > 3)
        {
            tx--;
            updated = true;
        }


        if (dy > 0.6 && ty < DImages[0]->height - 3)
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

void DSURF::FindPeaksOctaves(vector<vector<Image *> > &HessOctaves, vector<vector<Image *> > &LoGOctaves)
{
    vector<KeyPoint *> ps;
    for (unsigned int i = 0; i < HessOctaves.size(); i++)
    {
        ps = FindPeaksScales(i, HessOctaves[i], LoGOctaves[i]);
        kps.insert(kps.begin(), ps.begin(), ps.end());
        ps.clear();
    }

    return ;
}

vector<KeyPoint *> DSURF::FindPeaksScales(const int octave, vector<Image *> &HessImages, vector<Image *> &LoGImages)
{
    vector<KeyPoint *> peaks;
    Image * kpfound = new Image(HessImages[0]->width, HessImages[0]->height);
    float fx, fy, fs, funcVal, trace;
    int   x, y, border;
    unsigned int s;
    float contrast_thresh = (float) this->thresh;

    for (s = 1; s < HessImages.size()-1; s++)
    {
        border = Octs[octave][s]/2;
        for (y = border; y < (HessImages[0]->height - border); y++)
        {
            for (x = border; x < (HessImages[0]->width - border); x++)
            {
                funcVal = HessImages[s]->getPixel(x, y);
                if(funcVal <= 0)
                    continue;

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
                if(INTERP_KEYS)
                {
                    if (!InterpKey(x, y, s, HessImages, &fx, &fy, &fs, &funcVal))
                        continue;
                }
                fs  = fs < 0?0:fs;

                KeyPoint * peak = new KeyPoint();
                trace    = LoGImages[s]->getPixel(x, y);
                peak->x  = (int)round(fx);
                peak->y  = (int)round(fy);

                peak->dscale   = ((floor(Octs[octave][s]/6.0f) - 1) + fs)*DSURF::SIGMA;
                peak->gscale  = (floor(Octs[octave][s]/6.0f) - 1) + fs;
                peak->iscale  = peak->gscale*DSURF::mag;
                peak->funcVal = funcVal;
                peak->ori     = 0;
                peak->fscale  = fs;
                peak->sx      = fx;
                peak->sy      = fy;
                peak->octave  = octave;
                peak->KP      = true;
                peak->img_width   = this->crntimg->width;
                peak->div     = trace>0?1:-1;
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

bool DSURF::isSpatialPeak(Image *image, const int x, const int y)
{
    float center =image->getPixel(x, y);
    float p2 = 0;
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
        return false;
    }
}

bool DSURF::isScalePeak(Image * aimage, Image * bimage, Image * cimage, const int x, const int y)
{
    assert(aimage);
    assert(bimage);
    assert(cimage);

    Image *ims[3];
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
        return false;
    }
    return true;
}

void DSURF::BuildSURFOctaves(vector<vector<Image *> > &octaves, vector<vector<Image *> > &LoGOctaves)
{
    assert(intImg);
    int oct, sc, lobe, radius, x = 0, y = 0;
    double area, Dxx, Dyy, Dxy;
    unsigned int tpSize = 0;
    float  DetHess, trace;
    Image *crnt_img = NULL;
    Image *LoG_img  = NULL;

    for (oct = 0; oct < DSURF::MaxOctaves; oct++)
    {
        vector<Image *> HessianScales;
        vector<Image *> LoGScales;
        for(sc = 0 ; sc < DSURF::SCALESPEROCTAVE; sc++)
        {
            lobe   = pow(2, (oct+1))*(sc+1) + 1;
            area   = pow(3*lobe, 2);
            tpSize = Octs[oct][sc];
            radius = tpSize/2;
            crnt_img = new Image(this->intImg->width, this->intImg->height);
            LoG_img  = new Image(this->intImg->width, this->intImg->height);

            for(y = radius; y < (intImg->height - radius); y++)
            {
                for(x = radius; x < (intImg->width - radius); x++)
                {
                    Dxx     = this->getDxx(x, y, radius, tpSize);
                    Dyy     = this->getDyy(x, y, radius, tpSize);
                    Dxy     = this->getDxy(x, y, radius, tpSize);
                    Dxx     = Dxx/area;
                    Dyy     = Dyy/area;
                    Dxy     = Dxy/area;
                    DetHess = Dxx*Dyy - 0.81*Dxy*Dxy;
                    trace   = (Dxx + Dyy);
                    crnt_img->setPixel(x, y, DetHess);
                    LoG_img->setPixel(x, y, trace);
                }
                //exit(0);
            }
            HessianScales.push_back(crnt_img);
            LoGScales.push_back(LoG_img);
        }
        octaves.push_back(HessianScales);
        LoGOctaves.push_back(LoGScales);
    }
    return ;
}

void DSURF::FindOrientation(vector<KeyPoint*> &keyps)
{
    int  r = 0, c = 0, i = 0, j = 0, dx = 0, dy = 0, dr = 0;
    float resX[109], resY[109], Ang[109], sc = 1;
    const int id[13] = {6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6};
    float maxMag  = 0.f, orientation = 0.f, weight = 0.0f;
    float ang1 = 0.f, ang2 = 0.f, ang = 0, deltA = PI0/3.0;
    float sumX = 0.f, sumY = 0.f, mag = 0.0f;
    unsigned int k = 0, idx = 0;

    vector<KeyPoint*>::iterator vit;
    KeyPoint *crntKpt;

    for(vit = keyps.begin(); vit != keyps.end(); vit++)
    {
        crntKpt = *vit;
        memset(resX, 0, sizeof(float)*109);
        memset(resY, 0, sizeof(float)*109);
        memset(Ang,  0, sizeof(float)*109);

        idx = 0;
        sc = crntKpt->dscale;
        r  = (int)round(crntKpt->y);
        c  = (int)round(crntKpt->x);
        gMat = Filter::GaussianKernel2D(2.5*sc, 6);
        for(i = -6; i <= 6; i++)
        {
            for(j = -6; j <= 6; j++)
            {
                if(i*i + j*j < 36)
                {
                       weight = gMat[id[i+6]][id[j+6]];
                           dx = (int)round(j*sc);
                           dy = (int)round(i*sc);
                           dr = (int)round(2*sc);
                    resX[idx] = weight * this->getDx(c+dx, r+dy, dr, 2*dr);
                    resY[idx] = weight * this->getDy(c+dx, r+dy, dr, 2*dr);
                    ///Ang[idx]  = getAngle(j, i); ///it is not as good as the lower one
                     Ang[idx] = getAngle(resX[idx], resY[idx]);
                    idx++;
                }
            }
        }
        Cleaner::clear2DArray(gMat);

        sumX = sumY = mag = 0.0f;
        maxMag  = orientation = 0.0f;
        ang1 = ang2 = ang = 0;

        for(ang1 = 0; ang1 < PI2;  ang1 += 0.15f)
        {
            ang2 = ang1 + deltA;
            while(ang2 >= PI2)
            ang2 = ang2 - PI2;

            sumX = sumY = 0.f;
            for(k = 0; k < 109; ++k)
            {
                ang = Ang[k];
                if(ang1 < ang2 && ang1 < ang && ang < ang2)
                {
                    sumX += resX[k];
                    sumY += resY[k];
                }else if (ang2 < ang1 && ((ang > 0 && ang < ang2)
                                      || (ang > ang1 && ang < PI2)))
                {
                    sumX += resX[k];
                    sumY += resY[k];
                }
            }///for(k)
            mag = sumX*sumX + sumY*sumY;
            if (mag > maxMag)
            {
                maxMag = mag;
                orientation = getAngle(sumX, sumY);
            }
        }///for(Ang1)

        crntKpt->ori = orientation;
    }//for(vit)

    return ;
}

float DSURF::getAngle(float dx, float dy)
{
    if(dx == 0 && dy == 0)
    {
        return 0;
    }

    if(dx > 0 && dy >= 0)
        return atan(dy/dx);

    if(dx < 0 && dy >= 0)
        return PI0 - atan(-dy/dx);

    if(dx < 0 && dy < 0)
        return PI0 + atan(dy/dx);

    if(dx > 0 && dy < 0)
        return PI2 - atan(-dy/dx);

    return 0;
}

float DSURF::getDx(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral(x0, (y0 - rd0), rd0, ts)
                - this->intImg->boxIntegral((x0-rd0), (y0-rd0), rd0, ts);
    return val;
}

float DSURF::getDy(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral(x0-rd0, y0, ts, rd0)
                - this->intImg->boxIntegral((x0-rd0), (y0-rd0), ts, rd0);
    return val;
}

float DSURF::getDxx(const int x0, const int y0, const int rd0, const int ts)
{
    int rectx[4], recty[4];
    int d = ts/3;
    float posbox1, posbox2, negbox, val;
    rectx[0]  = rectx[1] = rectx[2] = rectx[3] = x0 - rd0;
    recty[0]  = recty[1] = recty[2] = recty[3] = y0 - rd0;
    rectx[0]  += 0;
    rectx[1]  += (d-1);
    rectx[2]  += 0;
    rectx[3]  += (d-1);
    recty[0]  += ((d+1)/2);
    recty[1]  += ((d+1)/2);
    recty[2]  += (ts-1-((d+1)/2));
    recty[3]  += (ts-1-((d+1)/2));

    posbox1    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox1   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    rectx[0]  += d;
    rectx[1]  += d;
    rectx[2]  += d;
    rectx[3]  += d;
    negbox     =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    negbox    +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    rectx[0]  += d;
    rectx[1]  += d;
    rectx[2]  += d;
    rectx[3]  += d;
    posbox2    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox2   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    val        = posbox2 + posbox1 - 2*negbox;

    return val;
}


float DSURF::getDyy(const int x0, const int y0, const int rd0, const int ts)
{
    int rectx[4], recty[4];
    int d = ts/3;
    float posbox1, posbox2, negbox, val;

    rectx[0]  = rectx[1] = rectx[2] = rectx[3] = x0 - rd0;
    recty[0]  = recty[1] = recty[2] = recty[3] = y0 - rd0;
    rectx[0]  += ((d+1)/2);
    rectx[1]  += ((d+1)/2);
    rectx[2]  += (ts-1-((d+1)/2));
    rectx[3]  += (ts-1-((d+1)/2));
    recty[0]  += 0;
    recty[1]  += (d-1);
    recty[2]  += 0;
    recty[3]  += (d-1);

    posbox1    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox1   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    recty[0]  += d;
    recty[1]  += d;
    recty[2]  += d;
    recty[3]  += d;
    negbox     =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    negbox    +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    recty[0]  += d;
    recty[1]  += d;
    recty[2]  += d;
    recty[3]  += d;
    posbox2    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox2   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    val        = posbox2 + posbox1 - 2*negbox;

    return val;
}

float DSURF::getDxy(const int x0, const int y0, const int rd0, const int ts)
{
    int rectx[4] = {0}, recty[4] = {0};
    int d = ts/3;
    float posbox1 = 0.0f, posbox2 = 0.0f, negbox1 = 0.0f, negbox2 = 0.0f, val = 0.0f;
    rectx[0]  = rectx[1] = rectx[2] = rectx[3] = x0 - rd0;
    recty[0]  = recty[1] = recty[2] = recty[3] = y0 - rd0;

    rectx[0]  += (d+1)/2 - 1;
    rectx[1]  += (rd0 - 1);
    rectx[1]  += (d+1)/2 - 1;
    rectx[3]  += (rd0 - 1);

    recty[0]  += (d+1)/2 - 1;
    recty[1]  += (d+1)/2 - 1;
    recty[2]  += (rd0 - 1);
    recty[3]  += (rd0 - 1);

    posbox1    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox1   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    rectx[0]  += (d+1);
    rectx[1]  += (d+1);
    rectx[2]  += (d+1);
    rectx[3]  += (d+1);
    negbox1    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    negbox1   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    rectx[0]  -= (d+1);
    rectx[1]  -= (d+1);
    rectx[2]  -= (d+1);
    rectx[3]  -= (d+1);
    recty[0]  += (d+1);
    recty[1]  += (d+1);
    recty[2]  += (d+1);
    recty[3]  += (d+1);
    negbox2    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    negbox2   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    rectx[0]  += (d+1);
    rectx[1]  += (d+1);
    rectx[2]  += (d+1);
    rectx[3]  += (d+1);
    posbox2    =  this->intImg->getPixel(rectx[0], recty[0]) - this->intImg->getPixel(rectx[1], recty[1]);
    posbox2   +=  this->intImg->getPixel(rectx[3], recty[3]) - this->intImg->getPixel(rectx[2], recty[2]);
    val        = posbox1 + posbox2 - negbox1 - negbox2;

    return val;
}

Image * DSURF::ScaleInitImage(Image * image)
{
    assert(image);
    Image *dst = NULL;

    if(RESIZE)
    {
        Image * im = Image::doubleSizeImage(image);
        dst = new Image(im->width, im->height);
        double sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA * 4);
        Filter::BlurImage(im, dst, sigma);
        delete im;
    }
    else
    {
        dst = new Image(image->width, image->height);
        double sigma = sqrt(SIGMA * SIGMA - INITSIGMA * INITSIGMA);
        Filter::BlurImage(image, dst, sigma);
    }

    return dst;
}

bool DSURF::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *oriImg = NULL, *tmpImg = NULL;

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
    tmpImg = Filter::ScaleImage(this->crntimg, 0, SIGMA, INITSIGMA);
    this->crntimg = tmpImg;
    vector<vector<Image*> > HessianOctaves;
    vector<vector<Image*> > LoGOctaves;

    this->intImg = this->buildIntImage(this->crntimg);
    this->BuildSURFOctaves(HessianOctaves, LoGOctaves);
    this->FindPeaksOctaves(HessianOctaves, LoGOctaves);

    this->FindOrientation(kps);

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

    if(this->myDescriptor != NULL)
    {
        delete this->crntimg;
        this->crntimg = oriImg;
    }

    if(strcmp(descfn,"")&&this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }

    if(strcmp(dvfn,"")&&this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }
    Cleaner::releaseOctaves(HessianOctaves);
    Cleaner::releaseOctaves(LoGOctaves);

    delete this->crntimg;
    delete this->intImg;
    this->intImg   = NULL;
    this->crntimg  = NULL;
    return true;
}

IntImage *DSURF::buildIntImage(Image *img0)
{
    unsigned int h = img0->height, w = img0->width;
    unsigned int i, j, loc;
    IntImage *intImg = new IntImage(w, h);
    double val = 0;
    int x, y;

    for(i = 0; i < h; i++)
    {
        y = i;
        for(j = 0; j < w; j++)
        {
            x = j;
            loc = i*w + j;     ///[x, y]
            val = intImg->getPixel(x, y-1) + intImg->getPixel(x-1, y)
                  + img0->pix[loc] - intImg->getPixel(x-1, y-1);
            intImg->setPixel(x, y, val);
        }
    }

    return intImg;
}

DSURF::~DSURF()
{
    vector<vector<float> >::iterator it;
    for(it = gMat.begin(); it != gMat.end(); it++)
    {
        vector<float> &crnt_vect = *it;
        crnt_vect.clear();
    }
    gMat.clear();
    delete [] featsBin;
    delete [] descWin;
}

void DSURF::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new DSURF();
    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 4; i++)
    {
        sprintf(vsrc,  "/home/wlzhao/datasets/trec03/img%d.jpg",   i);
        sprintf(vdst,  "/home/wlzhao/datasets/trec03/img%d.keys",  i);
        sprintf(vdesc, "/home/wlzhao/datasets/trec03/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }
}
