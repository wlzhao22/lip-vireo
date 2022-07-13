#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>

#include "keypoint.h"
#include "cleaner.h"
#include "filter.h"
#include "harris.h"
#include "vmath.h"

using namespace std;

const float Harris::sigma0  = 1.2f;
const float Harris::_SIGMA  = 1.4f;
const float Harris::k0      = 0.06f;
const int Harris::BORDER    = 8;
const int Harris::SCALES    = 4;
const int Harris::THRESH    = 720;
const float Harris::mag     = 4.0;
const float Harris::dfactor = 1.125f;

const int Harris::DEGREE    = 10;
const int Harris::NumOrient = 36;
const int Harris::DEGPERBIN = (360 / NumOrient);
const float Harris::NwKpThresh = 0.8;

Harris::Harris()
{
    cout<<"Detector ................................ Harris\n";
    this->DETECTOR   = harris;
    this->sel_option = THRSH;
    this->AFF_OUT    = false;
}

bool Harris::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    this->sigma   = Harris::sigma0;

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

bool Harris::isSpatialPeak(Image *image, const int x, const int y)
{
    float center =image->getPixel(x, y);
    int xi, yi;
    float p2;
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

vector<KeyPoint*> Harris::FindPeaksScales(const int octave, vector<Image*> HarImages)
{
    Image *kpfound = new Image(HarImages[0]->width, HarImages[0]->height);
    float contrast_thresh = this->thresh;
    float fx = 0.0f, fy = 0.0f, fs = 0.0f, funcVal = 0.0f;
    vector<KeyPoint *> peaks;
    unsigned int s = 0;
    int x = 0, y = 0;

    for (s = 1; s < HarImages.size()-1; s++)
    {
        for (y = BORDER; y < (HarImages[0]->height - BORDER); y++)
        {
            for (x = BORDER; x < (HarImages[0]->width - BORDER); x++)
            {
                funcVal = HarImages[s]->getPixel(x, y);
                if(funcVal <= 0)
                continue;

                funcVal = sqrt(funcVal);

                if(funcVal <= contrast_thresh)
                    continue;

                if(kpfound->getPixel(x, y) == 1)
                    continue;

                if(!isSpatialPeak(HarImages[s], x, y))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                KeyPoint * peak = new KeyPoint();

                peak->x = (int)round(fx * pow(2.0, octave));
                peak->y = (int)round(fy * pow(2.0, octave));

                peak->dscale    = Harris::_SIGMA * pow(2.0, octave + fs / (float) Harris::SCALES);
                peak->iscale    = peak->dscale*Harris::mag;
                peak->funcVal   = funcVal;
                peak->img_width = this->crntimg->width;
                peak->ori       = 0;
                peak->scale     = s;
                peak->fscale    = fs;
                peak->gscale    = octave + fs/Harris::SCALES;
                peak->sx        = fx;
                peak->sy        = fy;
                peak->octave    = octave;
                peak->KP        = true;
                peak->div       = 1;

                leveli_kps.push_back(peak);
                kpfound->setPixel((int)(fx + 0.5), (int)(fy + 0.5), 1);
            }
        }
        peaks.insert(peaks.begin(), leveli_kps.begin(), leveli_kps.end());
        leveli_kps.clear();
    }
    delete kpfound;
    return peaks;
}

vector<vector<Image *> > Harris::BuildHarrisOctaves(vector<vector<Image *> > & GxOctaves, vector<vector<Image *> > & GyOctaves)
{
    vector<vector<Image *> > HarrisOctaves;
    vector<Image *> gxScales, gyScales;
    vector<Image *> HarrisScales;

    for (unsigned int i = 0; i < GxOctaves.size(); i++)
    {
        gxScales = GxOctaves[i];
        gyScales = GyOctaves[i];
        HarrisScales = BuildHarrisScales(gxScales, gyScales);
        HarrisOctaves.push_back(HarrisScales);
    }

    return HarrisOctaves;
}

vector<Image* > Harris::BuildHarrisScales(vector<Image *>  &gxScales, vector<Image *> &gyScales)
{
    Image *crntDxImg = NULL, *crntDyImg = NULL, *crntDxDyImg = NULL, *crntHarrisImg = NULL;
    float dx2 = 0, dy2 = 0, Det = 0, trace = 0, harris = 0;
    float _sigma = 0, dsigma = 0, dsigma2 = 0;
    int x = 0, y = 0, width = 0, height = 0;
    float dx = 0, dy = 0, dxdy = 0;
    vector<Image* > HarrisScales;
    unsigned int  nscale = 0;
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
        crntDxDyImg = new Image(crntDxImg->width, crntDxImg->height);
        _sigma  = Harris::_SIGMA * pow(2.0, (nscale+1.0f) / (float) Harris::SCALES);
        dsigma  = _sigma*Harris::dfactor;
        dsigma2 = dsigma*dsigma;

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx   = crntDxImg->getPixel(x, y);
                dy   = crntDyImg->getPixel(x, y);
                dxdy = dx*dy*dsigma2;
                crntDxDyImg->setPixel(x, y, dxdy);
            }
        }

        crntDxImg->exp(2);
        crntDyImg->exp(2);
        crntDxImg->multiply(dsigma2);
        crntDyImg->multiply(dsigma2);
        Gkern = Filter::GaussianKernel1D(_sigma);
        Filter::Convolve1DWidth(Gkern,  crntDxImg, tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,    crntDxImg);
        Filter::Convolve1DWidth(Gkern,  crntDyImg, tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,    crntDyImg);
        Filter::Convolve1DWidth(Gkern,  crntDxDyImg,tmpimg);
        Filter::Convolve1DHeight(Gkern, tmpimg,      crntDxDyImg);
        Gkern.clear();
        crntHarrisImg = new Image(crntDxImg->width, crntDxImg->height);

        for(x = 0; x < width; x++)
        {
            for(y = 0; y < height; y++)
            {
                dx2    = crntDxImg->getPixel(x, y);
                dy2    = crntDyImg->getPixel(x, y);
                dxdy   = crntDxDyImg->getPixel(x, y);
                Det    = dx2*dy2 - dxdy*dxdy;
                trace  = dx2 + dy2;
                trace  = trace*trace;
                harris = Det - Harris::k0*trace;
                crntHarrisImg->setPixel(x, y, harris);
            }
        }

        HarrisScales.push_back(crntHarrisImg);
        delete crntDxDyImg;
    }

    delete tmpimg;

    return HarrisScales;
}


vector<KeyPoint*> Harris::FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint*> newkps;
    KeyPoint *newkp = NULL;

    float sigma, m, theta, degree;
    int c, index;
    int indexa, indexb, indexc, j;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval,maxp;

    bool valid = true;
    float thetas[NumOrient] = {0};
    float weight = 0.0f;
    float dsigma = Harris::_SIGMA*Harris::dfactor;

    for (i = 0; i < kps.size(); i++)
    {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *dsigma;
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

bool Harris::KeypointBuild(const char *fn, const char *dstfn,
                           const char *descfn, const char *dvfn)
{
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    this->crntimg = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    Image *tmpImg = NULL, *oriImg = NULL;

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

    int kpnum = 0;
    unsigned int ioctave  = 0;
    int width    = this->crntimg->width;
    int height   = this->crntimg->height;
    int octave   = width>height?height:width;
    float dsigma = Harris::dfactor*Harris::_SIGMA;

    vector< vector<Image*> > GxOctaves;
    vector< vector<Image*> > GyOctaves;
    vector< vector<Image*> > GaussOctaves;
    vector<KeyPoint *>::iterator it;
    vector<KeyPoint *> peaks;
    KeyPoint *crntPt = NULL;
    octave  = (int)floor(VMath::lgx(octave, 2) - 3);
    oriImg  = this->crntimg->clone();

    this->BuildOctaves(this->crntimg, _Dx_,    GxOctaves,    dsigma, Harris::SCALES);
    this->BuildOctaves(this->crntimg, _Dy_,    GyOctaves,    dsigma, Harris::SCALES);
    this->BuildOctaves(this->crntimg, _Conv2, GaussOctaves, dsigma, Harris::SCALES);

    vector< vector<Image*> > HarrisOctaves = this->BuildHarrisOctaves(GxOctaves, GyOctaves);
    Cleaner::releaseOctaves(GxOctaves);
    Cleaner::releaseOctaves(GyOctaves);


    for(ioctave = 0; ioctave < HarrisOctaves.size(); ioctave++)
    {
        peaks = this->FindPeaksScales(ioctave, HarrisOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    Cleaner::releaseOctaves(HarrisOctaves);

    switch(this->sel_option)
    {
    case 0:
    {
        AbstractDetector::topkSelect(kps, this->fix_kp_numb);
        break;
    }
    case 1:
    {
        AbstractDetector::topkEqDnSelect(kps, this->crntimg->width,
                                            this->crntimg->height, this->fix_kp_numb);
        break;
    }
    default:
    {
        break;
    }
    }

    vector<KeyPoint*> extra_peaks = this->FindOrientByGrad(kps, GaussOctaves);
    ///this->FindOrientation(kps);
    Cleaner::releaseOctaves(GaussOctaves);


    if(extra_peaks.size() > 0)
    {
        this->kps.insert(kps.begin(), extra_peaks.begin(), extra_peaks.end());
        extra_peaks.clear();
    }

    delete this->crntimg;
    this->crntimg = oriImg;

    if(strcmp(dstfn, ""))
    {
        (this->*saveKpts)(this->kps, kpnum, dstfn, this->resize_rate, this->AFF_OUT);
    }

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntPt = *it;
        crntPt->iscale = crntPt->iscale*Harris::mag;
    }

    if(strcmp(descfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kpnum, kps, descfn, this->resize_rate);
    }
    if(strcmp(dvfn,"") && this->myDescriptor!=NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kpnum, kps, dvfn, this->resize_rate);
    }

    delete oriImg;
    return true;
}

void Harris::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new Harris();
    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 4; i++)
    {
        sprintf(vsrc,  "/home/wlzhao/datasets/trec03/img%d.jpg",   i);
        sprintf(vdst,  "/home/wlzhao/datasets/trec03/img%d.keys",  i);
        sprintf(vdesc, "/home/wlzhao/datasets/trec03/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, vdst, vdesc, "");
    }
}
