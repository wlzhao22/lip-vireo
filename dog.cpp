#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "keypoint.h"
#include "descsift.h"
#include "cleaner.h"
#include "vmath.h"
#include "dog.h"

using namespace std;

const int DoG::MaxOctaves      = 5;
const float DoG::INITSIGMA     = 0.5;
const float DoG::_SIGMA        = 1.6;
const int DoG::SCALESPEROCTAVE = 6;

const int DoG::NumOrient      = 36;
const int DoG::DEGPERBIN      = (360 / NumOrient);
const float DoG::NwKpThresh   = 0.8;
const float DoG::LOW_CONTRAST_THRESH= 0.0350;
const float DoG::EDGE_PEAK_THRESH   = 10.0;
const float DoG::EDGE_PEAK_L        = 0.1;

const int DoG::INTERP_KEYS          = 1;

const int DoG::BORDER               = 5;

const int DoG::featLen              = 128;
const int DoG::ORIENTATION          = 8;
const int DoG::GRID                 = 4;
const int DoG::DSize                = GRID*GRID*ORIENTATION;
const int DoG::DEGREE               = 10;
const int DoG::PatchSize            = 41;
const float DoG::mag                = 8; //5
const int DoG::PatchMag             = 20;

const int DoG::kpnum0 = 750;

DoG::DoG()
{
    cout<<"Detector ................................ DoG\n";
    this->featsBin   = new float[128];
    this->descWin    = new float[PatchSize*PatchSize];
    this->DETECTOR   = dog;
    this->sel_option = THRSH;
    this->AFF_OUT    = false;
}

bool DoG::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    if(this->paras.find(argv[1]) != this->paras.end())
    {
        this->thresh = atof(paras["thresh"]);
    }else{
        this->thresh = LOW_CONTRAST_THRESH;
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

float DoG::InterpKeyStep(int x, int y, int s, vector<Image *> & DI, float * dx, float * dy, float * ds)
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


    float invDpp[3][3], val;
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

    val = DI[s]->getPixel(x, y);
    val += 0.5 * (Dp[0] * *ds + Dp[1] * *dy + Dp[2] * *ds);

    return fabs(val);
}

bool DoG::InterpKey(int x, int y, int s, vector<Image *> & DImages, float * fx, float * fy, float * fs,float *dogVal)
{
    bool addkey = true;
    int moves_left = 5;
    int tx = x;
    int ty = y;
    int ts = s;

    float dx, dy, ds,val;
    bool updated;

    do
    {
        moves_left--;
        updated = false;

        val = InterpKeyStep(tx, ty, ts, DImages, &dx, &dy, &ds);

        if (val < (float) LOW_CONTRAST_THRESH / (float) SCALESPEROCTAVE)
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

bool DoG::isPeak(Image * aimage, Image * bimage,Image * cimage, unsigned char &valley, int x, int y)
{

    assert(aimage);
    assert(bimage);
    assert(cimage);

    Image * ims[3];

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
        valley = 0;
    }

    else
    {

        /**/
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
        valley = 1;
        /**/
        //return false;
    }

    return true;
}

bool DoG::isEdgePeak(int x, int y, Image * image)
{
    assert(image);

    float r = EDGE_PEAK_THRESH;
    float b = pow(r + 1, 2) / r;
    float c = pow((EDGE_PEAK_L+1),2)/r;

    float D = image->getPixel(x, y);
    float Dxx = image->getPixel(x + 1, y) + image->getPixel(x - 1, y) - 2*D;
    float Dyy = image->getPixel(x, y + 1) + image->getPixel(x, y - 1) - 2*D;
    float Dxy = ((image->getPixel(x + 1, y + 1) - image->getPixel(x - 1, y + 1))
                 - (image->getPixel(x + 1, y - 1) - image->getPixel(x - 1, y - 1))) / 4.0;

    float TrH = Dxx + Dyy;
    float DetH = Dxx*Dyy - Dxy * Dxy;

    float a = TrH * TrH / DetH;

    if (a < b && a >= c)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void DoG::FindPeaksOctaves(vector<vector<Image *> > & DOctaves,vector<vector<Image*> > &GImage)
{
    vector<KeyPoint *> ps;
    for (unsigned int i = 0; i < DOctaves.size(); i++)
    {
        ps = FindPeaksScales(i, DOctaves[i],GImage[i]);
        kps.insert(kps.begin(), ps.begin(), ps.end());
        ps.clear();
    }

    return ;
}

vector<KeyPoint *> DoG::FindPeaksScales(int octave, vector<Image *> & DImages,vector<Image*> &GImage)
{
    vector<KeyPoint *> peaks;
    Image * kpfound = new Image(DImages[0]->width, DImages[0]->height);
    float fx = 0.0f, fy = 0.0f, fs = 0.0f, dogVal = 0.0f;
    unsigned char valley = 0;
    unsigned int s = 0;
    int x = 0, y = 0;

    float contrast_thresh = (float) this->thresh / (float) SCALESPEROCTAVE;

    if (INTERP_KEYS)
    {
        contrast_thresh *= 0.8;
    }

    for (s = 1; s < DImages.size() - 1; s++)
    {
        for (y = BORDER; y < DImages[0]->height - BORDER; y++)
        {
            for (x = BORDER; x < DImages[0]->width - BORDER; x++)
            {

                if (fabs(DImages[s]->getPixel(x, y)) <= contrast_thresh)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                if (!isPeak(DImages[s-1], DImages[s], DImages[s+1], valley, x, y))
                    continue;

                if (isEdgePeak(x, y, DImages[s]))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                dogVal = DImages[s]->getPixel(x, y);

                if (INTERP_KEYS)
                {
                    if (!InterpKey(x, y, s, DImages, &fx, &fy, &fs,&dogVal))
                        continue;
                }

                fs  = fs < 0?0:fs;

                KeyPoint * peak = new KeyPoint();
                peak->x      = (int)round(fx * pow(2.0, octave));
                peak->y      = (int)round(fy * pow(2.0, octave));
                peak->octave = octave;
                peak->scale  = s;
                peak->dscale = DoG::_SIGMA * pow(2.0, octave + fs / (float) SCALESPEROCTAVE);
                peak->iscale = peak->dscale*DoG::mag;
                peak->gscale = octave + fs / (float) SCALESPEROCTAVE;
                peak->funcVal = dogVal;
                peak->div     = dogVal > 0?1:-1;
                peak->img_width = this->crntimg->width;

                /**/
                if (this->RESIZE)
                {
                    peak->iscale /= 2.0;
                    //peak->gscale /=2.0;
                }
                /**/

                peak->scale  = s;
                peak->fscale = fs;

                peak->sx = fx;
                peak->sy = fy;
                peak->KP = true;
                peak->ori = 0;
                peak->octave = octave;

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

void DoG::BuildOctaves(Image * image,vector<vector<Image *> > &octaves)
{
    assert(image);
    int dim = min(image->height, image->width);
    int numoctaves = int (log((double) dim) / log(2.0)) - 2;
    numoctaves =numoctaves>MaxOctaves?MaxOctaves:numoctaves;
    Image * timage = image->clone();
    vector<Image *> scales;

    for (int i = 0; i < numoctaves; i++)
    {
        scales = BuildGaussianScales(timage);
        octaves.push_back(scales);
        Image * simage  = Image::halfSizeImage(scales[SCALESPEROCTAVE]);

        delete timage;
        timage = simage;

    }
    delete timage;
}

vector<Image *> DoG::BuildGaussianScales(Image * image)
{
    assert(image);
    vector<Image *> GScales;

    double k = pow(2, 1.0/(float)SCALESPEROCTAVE);
    GScales.push_back(image->clone());
    float sigma1,sigma2,sigma;
    int UPBOUND;
    Image *dstImg = NULL;
    UPBOUND = SCALESPEROCTAVE + 2;

    for (int i =  1; i < UPBOUND; i++)
    {
        dstImg = new Image(image->width, image->height);
        sigma1 = pow(k, i - 1) * DoG::_SIGMA;
        sigma2 = pow(k, i) * DoG::_SIGMA;
        sigma  = sqrt(sigma2*sigma2 - sigma1*sigma1);
        Filter::BlurImage(GScales[GScales.size() - 1], dstImg, sigma);
        GScales.push_back(dstImg);
    }
    return GScales;
}

vector<vector<Image *> > DoG::BuildDOGOctaves(vector<vector<Image *> > & GOctaves)
{
    vector<vector<Image *> > DOctaves;
    vector<Image *> scales;

    for (unsigned int i = 0; i < GOctaves.size(); i++)
    {
        scales = GOctaves[i];
        DOctaves.push_back(BuildDOGScales(GOctaves[i]));
    }

    return DOctaves;
}

vector<Image *> DoG::BuildDOGScales(vector<Image *> &LImages)
{
    vector<Image *> DImages;

    for (unsigned int i = 1; i < LImages.size(); i++)
    {
        Image * dog = new Image(LImages[i]->width, LImages[i]->height);
        LImages[i]->sub(LImages[i-1], dog);
        DImages.push_back(dog);
    }

    return DImages;
}

Image * DoG::ScaleInitImage(Image * image)
{
    assert(image);
    Image * dst;

    if (RESIZE)
    {
        Image * im = Image::doubleSizeImage(image);
        dst = new Image(im->width, im->height);
        double sigma = sqrt(DoG::_SIGMA * DoG::_SIGMA - INITSIGMA * INITSIGMA * 4);
        Filter::BlurImage(im, dst, sigma);
        delete im;
    }
    else
    {
        dst = new Image(image->width, image->height);
        double sigma = sqrt(DoG::_SIGMA * DoG::_SIGMA - INITSIGMA * INITSIGMA);
        Filter::BlurImage(image, dst, sigma);
    }

    return dst;
}


vector<KeyPoint *> DoG::FindOrientByGrad(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves)
{
    vector<KeyPoint * > newkps;
    vector<vector<float> > gmat;
    KeyPoint *newkp;
    float sigma, m, theta, degree;
    int c, index, j;
    int indexa, indexb, indexc;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval, maxp;

    bool valid;
    float thetas[NumOrient] = {0};
    float thetas0[NumOrient] = {0};
    float weight, ratio;

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALESPEROCTAVE) * DoG::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size()/2;
        for(j = 0; j < NumOrient; j++)
        {
            thetas[j] = 0;
        }
        ratio = 0;
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
                    weight = m * gmat[y + c][x + c];
                    thetas[index] += weight;
                }
            }
        }
        memcpy(thetas0, thetas, NumOrient*sizeof(float));
        vector<vector<float> >::iterator it;
        vector<float> crntvect;
        for(it = gmat.begin(); it != gmat.end(); it++)
        {
            crntvect = *it;
            crntvect.clear();
        }
        gmat.erase(gmat.begin(),gmat.end());

        for (j = 0; j < 6; j++)
            VMath::SmoothHistogram(thetas, NumOrient);

        maxval = VMath::maxVec(thetas, NumOrient);

        for (j = 0; j < NumOrient; j++)
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
                if(this->mydesc == FIND || this->mydesc == NFIND/**|| this->mydesc == NFIFT|| this->mydesc == FIFT**/)
                {
                    kps[i]->flip = AbstractDetector::findFlip(thetas0, j, NumOrient, ratio);
                }
            }
            else
            {
                /**/
                newkp = new KeyPoint();
                memcpy(newkp, kps[i], sizeof(KeyPoint));
                if(this->mydesc == ERIFT || this->mydesc == NERIFT ||
                        this->mydesc == NSPIN || this->mydesc == SPIN)
                {
                    newkp->KP = false;
                }
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NumOrient - PI;
                if(this->mydesc == FIND || this->mydesc == NFIND/**|| this->mydesc == NFIFT|| this->mydesc == FIFT**/)
                {
                    newkp->flip = AbstractDetector::findFlip(thetas0, j, NumOrient, ratio);
                }
                newkps.push_back(newkp);
                /**/
            }
        }
    }
    return newkps;
}

bool DoG::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
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
    tmpImg = Filter::ScaleImage(this->crntimg, 0, DoG::_SIGMA, INITSIGMA);
    this->crntimg = tmpImg;

    this->intImg = AbstractDetector::buildIntImage(this->crntimg);

    vector<vector<Image *> > GOctaves;
    BuildOctaves(this->crntimg,GOctaves);
    vector<vector<Image *> > DoGOctaves = BuildDOGOctaves(GOctaves);

    FindPeaksOctaves(DoGOctaves,GOctaves);

    vector<KeyPoint *> newkps;
    newkps = FindOrientByGrad(kps, GOctaves);

    kps.insert(kps.begin(), newkps.begin(), newkps.end());
    newkps.clear();
    stable_sort(kps.begin(),kps.end(),KeyPoint::keypCompF);

    /** only keep top 'kpnum0' keypoints according to DoG function**/
    switch(this->sel_option)
    {
    case 0:
    {
        AbstractDetector::topkSelect(kps, this->fix_kp_numb);
        break;
    }
    case 1:
    {
        AbstractDetector::topkEqDnSelect(kps,this->crntimg->width,this->crntimg->height, this->fix_kp_numb);
        break;
    }
    default:
    {
        break;
    }
    }

    Cleaner::releaseOctaves(DoGOctaves);

    if(strcmp(dstfn,""))
    {
        (this->*saveKpts)(this->kps,this->kps.size(),dstfn,this->resize_rate, this->AFF_OUT);
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

    Cleaner::releaseOctaves(GOctaves);

    delete this->crntimg;
    delete this->intImg;
    this->crntimg = NULL;
    this->intImg  = NULL;
    return true;
}

void DoG::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new DoG();

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

DoG::~DoG()
{
    delete [] featsBin;
    delete [] descWin;
}
