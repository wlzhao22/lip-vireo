#include "log.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

/**constant initialization**/

const int LoG::MaxOctaves   = 5;
const int LoG::SCALES       = 6;
const float LoG::_SIGMA     = 1.4;
const bool  LoG::INTERP_KEYS = true;
const float LoG::mag        = 4.50; //4.5
const int LoG::BORDER       = 5;
const int LoG::THRESH       = 45; //55

const float LoG::INITSIGMA  = 0.5;

const int LoG::DEGREE       = 10;
const int LoG::NumOrient    = 36;
const int LoG::DEGPERBIN    = (360 / NumOrient);
const float LoG::NwKpThresh = 0.8;

const float LoG::cvtRatio   = 10.0f;
const int LoG::maxIter      = 10;
const int LoG::W0           = 10;
const float LoG::err_c      = 0.80f;

/**
*            lap = _sigma*_sigma*(dxx+dyy);
*            lap = _sigma*(dxx+dyy); //this one is better practically
*
**/

LoG::LoG()
{
    cout<<"Detector ................................ LoG\n";
    this->DETECTOR   = _log;
    this->sel_option = THRSH;
}

bool LoG::paramsCheck()
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
    }

    return true;
}

float LoG::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
{
    /// first derivative of D with respect to x, y, s
    float Dp[3] = {0};
    /// Hessian of D
    float Dpp[3][3] = {{0}};

    Dp[0] = (DI[s]->getPixel(x+1, y) - DI[s]->getPixel(x-1, y)) / 2.0; // Dx
    Dp[1] = (DI[s]->getPixel(x, y+1) - DI[s]->getPixel(x, y-1)) / 2.0; // Dy
    Dp[2] = (DI[s+1]->getPixel(x, y) - DI[s-1]->getPixel(x, y)) / 2.0; // Ds

    /**
      Hessian (3x3) matrix is defined as follows:
      Symmetric matrix

      Dxx Dxy Dxs
      Dyx Dyy Dys
      Dsx Dsy Dss
    **/

    Dpp[0][0] = (DI[s]->getPixel(x+1, y) + DI[s]->getPixel(x-1, y)
                 - 2.0 * DI[s]->getPixel(x, y));

    Dpp[1][1] = (DI[s]->getPixel(x, y+1) + DI[s]->getPixel(x, y-1)
                 - 2.0 * DI[s]->getPixel(x, y));

    /// Dzz
    Dpp[2][2] = (DI[s+1]->getPixel(x, y) + DI[s-1]->getPixel(x, y)
                 - 2.0 * DI[s]->getPixel(x, y));


    /// Dxy = Dyx
    Dpp[0][1] = Dpp[1][0] = (DI[s]->getPixel(x+1, y+1) - DI[s]->getPixel(x-1, y+1)
                             - DI[s]->getPixel(x+1, y-1) + DI[s]->getPixel(x-1, y-1)) / 4.0;

    /// Dxs = Dsx
    Dpp[0][2] = Dpp[2][0] = (DI[s+1]->getPixel(x+1, y) - DI[s+1]->getPixel(x-1, y)
                             - DI[s-1]->getPixel(x+1, y) + DI[s-1]->getPixel(x-1, y)) / 4.0;

    /// Dys = Dsy
    Dpp[1][2] = Dpp[2][1] = (DI[s+1]->getPixel(x, y+1) - DI[s+1]->getPixel(x, y-1)
                             - DI[s-1]->getPixel(x, y+1) + DI[s-1]->getPixel(x, y-1)) / 4.0;


    float invDpp[3][3];
    VMath::mInv33(Dpp, invDpp);

    /// Solve for delta positions
    *dx = 0;
    for (int i = 0; i < 3; i++)
        *dx -= invDpp[0][i] * Dp[i];

    *dy = 0;
    for (int i = 0; i < 3; i++)
        *dy -= invDpp[1][i] * Dp[i];

    *ds = 0;
    for (int i = 0; i < 3; i++)
        *ds -= invDpp[2][i] * Dp[i];

    float val = DI[s]->getPixel(x, y);
    val += 0.5 * (Dp[0] * *ds + Dp[1] * *dy + Dp[2] * *ds);

    return fabs(val);
}

bool LoG::InterpKey(int x, int y, int s, vector<Image *> & LoGImages, float * fx, float * fy, float * fs,float *dogVal)
{
    bool addkey = true;
    int moves_left = 5;
    int tx = x;
    int ty = y;
    int ts = s;

    float dx, dy, ds,val;
    bool updated;
    float contrast_thresh = this->thresh*0.8;

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

vector<vector<Image *> > LoG::BuildLOGOctaves(vector<vector<Image *> > & GxOctaves,vector<vector<Image *> > & GyOctaves)
{
    vector<vector<Image *> > LoGOctaves;
    vector<Image *> gxScales;
    vector<Image *> gyScales;
    vector<Image *> LoGScales;

    for (unsigned int i = 0; i < GxOctaves.size(); i++)
    {
        gxScales = GxOctaves[i];
        gyScales = GyOctaves[i];
        LoGScales = BuildLoGScales(gxScales,gyScales);
        LoGOctaves.push_back(LoGScales);
    }

    return LoGOctaves;
}

vector<Image* > LoG::BuildLoGScales(vector<Image *>  &gxScales,vector<Image *> &gyScales)
{
    vector<Image* > LoGScales;
    Image *crntDxxImage = NULL, *crntDyyImage = NULL, *crntLoGImage = NULL;
    int x = 0, y = 0, width = 0, height = 0;
    float _sigma = 1.0f, lap = 0;
    unsigned int nscale = 0;
    float dxx = 0, dyy = 0;

    for(nscale = 0; nscale < gxScales.size(); nscale++)
    {
        crntDxxImage = gxScales[nscale];
        crntDyyImage = gyScales[nscale];
        width  = crntDxxImage->width;
        height = crntDxxImage->height;
        crntLoGImage = new Image(crntDxxImage->width,crntDxxImage->height);
        _sigma = _SIGMA * pow(2.0,  nscale / (float) SCALES);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dxx = crntDxxImage->getPixel(x,y);
                dyy = crntDyyImage->getPixel(x,y);
                ///lap = _sigma*_sigma*(dxx+dyy);
                lap = _sigma*(dxx+dyy);
                crntLoGImage->setPixel(x,y,lap);
            }
        }
        LoGScales.push_back(crntLoGImage);
    }
    return LoGScales;
}

bool LoG::adaptAffine(vector<vector<Image *> > &LoGOctaves, vector<vector<Image *> > & GxxOctaves,vector<vector<Image *> > & GyyOctaves,vector<vector<Image *> >& GxyOctaves)
{
    vector<KeyPoint*>::iterator it;
    KeyPoint *crntKp = NULL;
    Image *crnt_Gyy = NULL, *crnt_Gxx = NULL, *crnt_Gxy = NULL;
    float u[2][2], ei[2], mi[2][2];
    float a, b, c, xn, yn;
    int  is, x, y, ioctave = 0;
    unsigned int i = 0;

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

bool LoG::isPeak(Image * aimage, Image * bimage,Image * cimage, short &valley, int x, int y)
{
    assert(aimage);
    assert(bimage);
    assert(cimage);


    Image *ims[3] = {NULL};

    ims[0] = aimage;
    ims[1] = bimage;
    ims[2] = cimage;

    float center = bimage->getPixel(x, y);

    if (center > 0.0)
    {

        for (int si = 0; si < 3; si++)
        {
            for (int yi = y - 1; yi <= y + 1; yi++)
            {
                for (int xi = x - 1; xi <= x + 1; xi++)
                {
                    float p2 = ims[si]->getPixel(xi, yi);

                    if (center < p2)
                        return false;

                }
            }
        }
        valley = 1;
        return true;
    }else
    {
        for (int si = 0; si < 3; si++)
        {
            for (int yi = y - 1; yi <= y + 1; yi++)
            {
                for (int xi = x - 1; xi <= x + 1; xi++)
                {
                    float p2 = ims[si]->getPixel(xi, yi);

                    if (center > p2)
                        return false;
                }
            }
        }
        valley = -1;
        return true;
    }
}

vector<KeyPoint *> LoG::FindPeaksScales(int octave, vector<Image *> & LoGImages)
{
    vector<KeyPoint *> peaks;
    Image *kpfound = new Image(LoGImages[0]->width, LoGImages[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0;
    short valley;

    float contrast_thresh = this->thresh;///(float) NUM_SCALES;

    /**
    if(INTERP_KEYS)
        contrast_thresh *= 0.8;
    **/

    for (unsigned int s = 1; s < LoGImages.size() - 1; s++)
    {
        for (int y = BORDER; y < LoGImages[0]->height - BORDER; y++)
        {

            for (int x = BORDER; x < LoGImages[0]->width - BORDER; x++)
            {
                funcVal = LoGImages[s]->getPixel(x, y);
                funcVal = fabs(funcVal);

                if (funcVal <= contrast_thresh)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                if (!isPeak(LoGImages[s-1], LoGImages[s], LoGImages[s+1], valley, x, y))
                    continue;

                /**
                if (isEdgePeak(x, y, LoGImages[s]))
                	continue;
                **/

                fx = x;  fy = y;  fs = s;

                /**/
                if (INTERP_KEYS)
                {
                    if (!InterpKey(x, y, s, LoGImages, &fx, &fy, &fs,&funcVal))
                        continue;
                }
                /**/

                fs  = fs < 0?0:fs;

                KeyPoint * peak = new KeyPoint();

                peak->x = (int)round(fx * pow(2.0, octave));
                peak->y = (int)round(fy * pow(2.0, octave));

                peak->dscale  = LoG::_SIGMA * pow(2.0, octave + fs / (float) SCALES);
                peak->iscale  = peak->dscale*LoG::mag;
                peak->gscale  = octave + fs/LoG::SCALES;
                peak->funcVal = funcVal;
                peak->ori     = 0;
                peak->div     = (int)valley;

                peak->scale   = s;
                peak->fscale  = fs;
                peak->sx      = fx;
                peak->sy      = fy;
                peak->octave  = octave;
                peak->KP      = true;
                peak->img_width = this->crntimg->width;

                leveli_kps.push_back(peak);
                kpfound->setPixel((int)(fx + 0.5), (int)(fy + 0.5), 1);
            }
        }
        peaks.insert(peaks.begin(),leveli_kps.begin(),leveli_kps.end());
        leveli_kps.clear();
    }
    //this->crntimg = tmpimg;
    delete kpfound;
    return peaks;
}

vector<KeyPoint *> LoG::FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves)
{
    vector<KeyPoint * > newkps;
    vector<vector<float> > gmat;
    KeyPoint *newkp;

    float sigma, m, theta, degree;
    int   c, index, j, x, y;
    int   indexa, indexb, indexc;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval, maxp;
    float weight, angle, curl = 0;

    bool valid;
    float thetas[NumOrient] = {0};

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size()/2;

        for(j = 0; j < NumOrient; j++)
        {
            thetas[j] = 0;
        }
        for (y = -c; y <= c; y++)
        {
            for (x = -c; x <= c; x++)
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
                    weight = m*gmat[y + c][x + c];
                    thetas[index] +=  weight;

                    if(y ==0 && x ==0)
                        continue;

                    angle = theta - atan2(y, x);

                    if(angle >=0)
                    {
                        if(angle <=PI)
                        {
                            curl -= fabs(sin(angle))*weight;
                        }
                        else
                        {
                            curl += fabs(sin(angle))*weight;
                        }
                    }
                    else
                    {
                        if(angle >= -PI )
                        {
                            curl += fabs(sin(angle))*weight;
                        }
                        else
                        {
                            curl -= fabs(sin(angle))*weight;
                        }
                    }
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

        if(curl > 0.0f)
        {
            kps[i]->flip = 1;
        }
        else
        {
            kps[i]->flip = -1;
        }

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
                if(this->mydesc == ERIFT || this->mydesc == NERIFT ||
                        this->mydesc == NSPIN || this->mydesc == SPIN)
                {
                    newkp->KP = false;
                }
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NumOrient - PI;
                newkps.push_back(newkp);
            }
        }
    }
    return newkps;
}

bool LoG::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn)
{
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *oriImg = NULL,*tmpImg = NULL;

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
    this->intImg  =  AbstractDetector::buildIntImage(this->crntimg);

    vector<vector<Image*> >  GxxOctaves, GyyOctaves, GxyOctaves;
    vector<vector<Image *> > LoGOctaves;
    vector<vector<Image *> > GaussOctaves;
    vector<Image *> LoGScales;
    vector<KeyPoint *> peaks;
    unsigned int ioctave;

    this->BuildOctaves(this->crntimg, _Dxx_, GxxOctaves, LoG::_SIGMA, LoG::SCALES);
    this->BuildOctaves(this->crntimg, _Dyy_, GyyOctaves, LoG::_SIGMA, LoG::SCALES);
    this->BuildOctaves(this->crntimg, _Dxy_, GxyOctaves, LoG::_SIGMA, LoG::SCALES);

    LoGOctaves = this->BuildLOGOctaves(GxxOctaves,GyyOctaves);

    for(ioctave = 0; ioctave < LoGOctaves.size(); ioctave++)
    {
        peaks = this->FindPeaksScales(ioctave, LoGOctaves[ioctave]);
        this->kps.insert(this->kps.begin(),peaks.begin(),peaks.end());
        peaks.clear();
    }

    this->BuildOctaves(this->crntimg, _Conv2, GaussOctaves, LoG::_SIGMA, LoG::SCALES);

    if(this->AFF_OUT)
    {
        this->adaptAffine(LoGOctaves, GxxOctaves, GyyOctaves, GxyOctaves);

    }

    vector<KeyPoint *> nwkps;
    nwkps = this->FindOrientByGrad(kps, GaussOctaves);

    this->kps.insert(kps.begin(), nwkps.begin(), nwkps.end());
    nwkps.clear();

    Cleaner::releaseOctaves(GaussOctaves);

    delete this->crntimg;
    this->crntimg = oriImg;

    switch(this->sel_option)
    {
    case 0:
    {
        AbstractDetector::topkSelect(kps, this->fix_kp_numb);
        break;
    }
    case 1:
    {
        AbstractDetector::topkEqDnSelect(kps,this->crntimg->width, this->crntimg->height, this->fix_kp_numb);
        break;
    }
    default:
    {
        //AbstractDetector::topkSelect(kps,this->fix_kp_numb);
        break;
    }
    }

    Cleaner::releaseOctaves(GxxOctaves);
    Cleaner::releaseOctaves(GyyOctaves);
    Cleaner::releaseOctaves(GxyOctaves);
    Cleaner::releaseOctaves(LoGOctaves);

    if(strcmp(dstfn, ""))
    {
        (this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
    }

    if(strcmp(descfn, "")&&this->myDescriptor != NULL)
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

void LoG::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new LoG();
    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 6; i++)
    {
        sprintf(vsrc,  "/home/wlzhao/datasets/vgg/graf/graf%d.jpg",   i);
        sprintf(vdst,  "/home/wlzhao/datasets/vgg/graf/graf%d.keys",  i);
        sprintf(vdesc, "/home/wlzhao/datasets/vgg/graf/graf%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }

}
