#include "hessian.h"

#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

/**constant initialization**/

const int Hessian::MaxOctaves    = 4;
const int Hessian::SCALES        = 6;
const float Hessian::_SIGMA      = 1.4;
const float Hessian::INITSIGMA   = 0.5;
const bool  Hessian::INTERP_KEYS = true;

const float Hessian::mag         = 3.5; //optimal 8
const int Hessian::BORDER        = 5;
const int Hessian::THRESH        = 600; //100

const int Hessian::DEGREE        = 10;
const int Hessian::NOrient       = 36;
const int Hessian::DEGPERBIN     = (360 / NOrient);
const float Hessian::NwKpThresh  = 0.8;

const float Hessian::cvtRatio    = 10.0f;

const int Hessian::maxIter       = 10;
const int Hessian::W0            = 10;
const float Hessian::err_c       = 0.80f;

/**
* I try to use a dynamic threshold, however, no significant improvement I can get,
it actually produces a few bigger number of interest points

@date Mar.2nd.2009

const float Hessian::thresh_ratio = 2.5f;
**/

/**
I modified the implementation completely. Now it performs as good as Hessian-Laplacian.
However, it is really Hessian detector

@date 11-May-2011
@author Wan-Lei Zhao
**/

/**
I modified the implementation a little bit, it performs better and faster

@date 11-May-2015
@author Wan-Lei Zhao
**/


Hessian::Hessian()
{
    cout<<"Detector ................................ Hessian\n";
    this->DETECTOR = hessian;
    this->intImg   = NULL;
    this->sel_option = THRSH;
}

bool Hessian::paramsCheck()
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

float Hessian::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
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

bool Hessian::InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal)
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

vector<vector<Image *> > Hessian::BuildLOGOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves)
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

vector<Image* > Hessian::BuildLoGScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales)
{
    Image *crntDxxImage = NULL, *crntDyyImage = NULL, *crntLoGImage = NULL;
    int x = 0, y = 0, width = 0, height = 0;
    float _sigma = 1.0f, lap = 0.0f;
    float dxx = 0.0f, dyy = 0.0f;
    vector<Image* > LoGScales;
    unsigned int nscale = 0;

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
                _sigma = Hessian::_SIGMA * pow(2.0,  nscale / (float) SCALES);
                lap = _sigma*_sigma*(dxx+dyy);
                ///lap = lap<0?(-1*lap):lap;
                crntLoGImage->setPixel(x,y,lap);
            }
        }
        LoGScales.push_back(crntLoGImage);
    }

    return LoGScales;
}

vector<vector<Image *> > Hessian::BuildHessOctaves(vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves)
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

vector<Image* > Hessian::BuildHessScales(vector<Image *>  &gxxScales,vector<Image *> &gyyScales,vector<Image *> &gxyScales)
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
                _sigma = Hessian::_SIGMA * pow(2.0,  nscale / (float) SCALES);
                det = _sigma*_sigma*(dxx*dyy - dxy*dxy);

                crntHessImg->setPixel(x, y, det);
            }
        }
        HessScales.push_back(crntHessImg);
    }

    return HessScales;
}

bool Hessian::isSpatialPeak(Image *image,const int x,const int y)
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

bool Hessian::isScalePeak(Image * aimage, Image * bimage,Image * cimage, const int x, const int y)
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

vector<KeyPoint *> Hessian::FindPeaksScales(const int octave, vector<Image*> &HessImages, vector<Image*> &GxxImages, vector<Image*> &GyyImages)
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
                peak->dscale  = Hessian::_SIGMA * pow(2.0, octave + fs / (float) SCALES);
                peak->iscale  = peak->dscale*Hessian::mag;
                peak->funcVal = funcVal;
                peak->ori     = 0;
                peak->scale   = s;
                peak->fscale  = fs;
                peak->gscale  = octave + fs/Hessian::SCALES;
                peak->sx      = x;
                peak->sy      = y;
                peak->octave  = octave;
                peak->KP      = true;
                peak->div     = trace > 0 ? 1:-1;
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

vector<KeyPoint*> Hessian::FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint*> newkps;
    KeyPoint *newkp = NULL;
    int indexa = 0, indexb = 0, indexc = 0, j = 0, c = 0, index = 0, x = 0, y = 0;
    float sigma = 0.0f, m = 0.0f, theta = 0.0f, degree = 0.0f;
    float maxval = 0.0f, maxp = 0.0f, weight = 0, Wi = 1.0f;
    float thetaa = 0.0f, thetab = 0.0f, thetac = 0.0f;
    float thetas[NOrient] = {0};
    unsigned int i = 0;
    bool valid = false;

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES)*Hessian::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
         c = gmat.size()/2;
        Wi = 3.0f*sigma;
        Wi = Wi > 1.0f ? Wi:1.0f;
        Wi = Wi*Wi;

        memset(thetas, 0, sizeof(float)*NOrient);
        for(y = -c; y <= c; y++)
        {
            for(x = -c; x <= c; x++)
            {
                if ((x*x + y*y) > Wi)
                    continue;

                valid = Filter::GetPixOrientation((int) (kps[i]->sx + x + 0.5), (int) (kps[i]->sy + y + 0.5),
                                                  GOctaves[kps[i]->octave][kps[i]->scale], m, theta);

                if(valid)
                {
                    degree = theta / PI * 180.0 + 180.0;
                    index  = ((int) (degree / DEGREE));
                    index  = index%NOrient;
                    weight = m*gmat[y + c][x + c];
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
        gmat.erase(gmat.begin(), gmat.end());

        for(j = 0; j < 6; j++)
        {
           VMath::SmoothHistogram(thetas, NOrient);
        }

        maxval = VMath::maxVec(thetas, NOrient);

        for (j = 0; j < NOrient; j++)
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

bool Hessian::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *oriImg = NULL, *tmpImg = NULL;

    vector<vector<Image*> >  GxxOctaves, GyyOctaves, GxyOctaves;
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

    this->BuildOctaves(this->crntimg, _Dxx_, GxxOctaves, Hessian::_SIGMA, Hessian::SCALES);
    this->BuildOctaves(this->crntimg, _Dyy_, GyyOctaves, Hessian::_SIGMA, Hessian::SCALES);
    this->BuildOctaves(this->crntimg, _Dxy_, GxyOctaves, Hessian::_SIGMA, Hessian::SCALES);
    this->BuildOctaves(this->crntimg, _Conv2, GaussOctaves, Hessian::_SIGMA, Hessian::SCALES);
    HessOctaves = this->BuildHessOctaves(GxxOctaves, GyyOctaves, GxyOctaves);

    for(ioctave = 0; ioctave < HessOctaves.size(); ioctave++)
    {
        vector<KeyPoint *> peaks = this->FindPeaksScales(ioctave, HessOctaves[ioctave],
                                                         GxxOctaves[ioctave], GyyOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    Cleaner::releaseOctaves(GxxOctaves);
    Cleaner::releaseOctaves(GyyOctaves);
    Cleaner::releaseOctaves(GxyOctaves);
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
        //AbstractDetector::topkSelect(kps,this->fix_kp_numb);
        break;
    }
    }

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

void Hessian::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new Hessian();
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

Hessian::~Hessian()
{
}
