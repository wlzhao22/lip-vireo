#include "dense.h"

#include "cleaner.h"
#include "filter.h"
#include "vmath.h"


#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

/**Constants Initialization**/

const int Dense::MaxOctaves   = 1;
const int Dense::SCALES       = 1;
const float Dense::_SIGMA     = 1.1;
const float Dense::INITSIGMA  = 0.5;
const float Dense::mag        = 12.5f; //12;
const int Dense::Delta        = 5; //10

const int Dense::THRESH       = 11;

const int Dense::DEGREE       = 10;
const int Dense::NumOrient    = 36;
const int Dense::DEGPERBIN    = (360 / NumOrient);
const float Dense::NwKpThresh = 0.8;
const float Dense::k0         = 0.06f;

/**
Dense sampling has been widely accepted as an important
feature type for image classfication tasks. People believe
that for image categorization tasks, simple dense sampling
outperforms any kind of detectors. I have one implementation
which tries to sample points in multiple octave and multiple
scales.

@author Wanlei Zhao
@date 22/Mar/2011
**/


Dense::Dense()
{
    cout<<"Detector ................................ Dense\n";
    this->DETECTOR = dense;
    this->AFF_OUT  = false;
    this->BORDER   = (int)ceil(Dense::_SIGMA*Dense::mag);
}

bool Dense::paramsCheck()
{
    const char *argv[2] = {"step", "sigma"};

    if(paras.find(argv[0]) != paras.end())
    {
        this->step = atof(paras["step"]);
        cout<<"Step size ............................... "<<this->step<<endl;
    }else{
        this->step = Dense::Delta;
        cout<<"Step size ............................... "<<Dense::Delta<<" (by default)\n";
    }

    if(paras.find(argv[1]) != paras.end())
    {
        this->sigma0 = atof(paras["sigma"]);
        cout<<"sigma ................................... "<<this->sigma0<<endl;
    }else{
        this->sigma0 = Dense::_SIGMA;
        cout<<"sigma ................................... "<<Dense::_SIGMA<<" (by default)\n";
    }

    return true;
}

vector<KeyPoint *> Dense::generatePeaksScales(const int octave, vector<Image *> & GScales)
{
    vector<KeyPoint *> peaks;

    Image * kpfound = new Image(GScales[0]->width, GScales[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0;
    int x = 0, y = 0, ux = 0, uy = 0;
    unsigned int delta = this->step;
    unsigned int s = 0;

	for (s = 0; s < GScales.size(); s++)
	{
		for (y = BORDER; y < (GScales[0]->height - BORDER); y += delta)
		{
			for (x = BORDER; x < (GScales[0]->width - BORDER); x += delta)
			{
			    funcVal = GScales[s]->getPixel(x, y);

				fx = x; fy = y; fs = s;

				KeyPoint *peak = new KeyPoint();

			    peak->x = (int)round(fx * pow(2.0, octave));
			    peak->y = (int)round(fy * pow(2.0, octave));

				peak->dscale  = this->sigma0 * pow(2.0, octave + fs / (float) SCALES);
				peak->iscale  = peak->dscale*Dense::mag;
				peak->funcVal = funcVal;

				peak->scale   = s;
				peak->fscale  = fs;
                peak->gscale  = VMath::lgx(peak->iscale, 2.0f);
				peak->sx      = fx;
				peak->sy      = fy;
				peak->octave  = octave;
				peak->KP      = true;
				peak->ori     = 0;
				peak->a       = 1;
				peak->b       = 0;
				peak->c       = 1;
				ux = floor(fx+0.5);
				uy = floor(fy+0.5);

				leveli_kps.push_back(peak);
                kpfound->setPixel(ux, uy, 1);
			}
		}
		peaks.insert(peaks.begin(),leveli_kps.begin(),leveli_kps.end());
		leveli_kps.clear();
		///cout<<"scales: "<<s<<"\t"<<GScales[0]->height<<endl;
	}
	delete kpfound;

	return peaks;
}

void Dense::BuildOctaves(Image * image, GBlur blur_opt, vector<vector<Image *> > &Ggoctaves, const float _SIGMA0)
{
    assert(image);
    int dim        = min(image->height, image->width);
    int numoctaves = 1;
    while((dim/2 - 41) > 0)
    {
        numoctaves++;
        dim = dim/2;
    }
    numoctaves     = numoctaves>MaxOctaves?MaxOctaves:numoctaves;
    Image *timage  = image->clone();
    Image *simage  = NULL;
    vector<Image *> imgScales;
    int i = 0;

    for (i = 0; i < numoctaves; i++)
    {
        imgScales = BuildGaussianScales(timage, blur_opt, i, _SIGMA0);
        Ggoctaves.push_back(imgScales);
        Filter::BlurImage(timage, _SIGMA0);
        simage = Image::halfSizeImage(timage);
        delete timage;
        timage = simage;
    }
    delete timage;
    timage = NULL;
}

vector<Image*> Dense::BuildGaussianScales(Image * image, GBlur blur_opt, const int noctave,const float _SIGMA0)
{
    assert(image);
    vector<Image*> Scales;

    double k = pow(2, 1.0/(float)Dense::SCALES);
    float sigma1, sigma2, sigma;
    Image *dggImg = NULL;

    vector<float> Gkern;
    vector<float> dkern;

    for (int i =  1; i <= Dense::SCALES; i++)
    {
        dggImg = new Image(image->width, image->height);
        sigma1 = pow(k, i - 1) *_SIGMA0;
        sigma2 = pow(k, i) *_SIGMA0;

        if(blur_opt == _Conv2)
        {
            sigma = sqrt(sigma2*sigma2 - sigma1*sigma1);
            if(Scales.size() == 0)
            {
                Filter::BlurImage(image, dggImg, sigma);
            }
            else
            {
                Filter::BlurImage(Scales[Scales.size()-1], dggImg, sigma);
            }
        }else if(blur_opt == _Dx_)
        {
            dkern = Filter::GaussianDxKernel1D(sigma2);
            Filter::Convolve1DWidth(dkern, image, dggImg);
            dkern.clear();
        }else if(blur_opt == _Dy_)
        {
            dkern = Filter::GaussianDxKernel1D(sigma2);
            Filter::Convolve1DHeight(dkern, image, dggImg);
            dkern.clear();
        }
        Scales.push_back(dggImg);
    }//end-for(int i)

    return Scales;
}

vector<KeyPoint *> Dense::FindOrientation(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves)
{
	vector<KeyPoint * > newkps;
	vector<vector<float> > gmat;
	KeyPoint *newkp = NULL;

	float sigma, m, theta, degree;
	int c, index, indexa, indexb, indexc, j;
	float thetaa, thetab, thetac;
	unsigned int i;
	float maxval, maxp;
	bool valid;
	float thetas[NumOrient] = {0};
	float weight;

	for (i = 0; i < kps.size(); i++)
	{
		sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *this->sigma0;
		Filter::GaussianKernel2D(sigma, gmat);
		c     = gmat.size()/2;

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

     for(j = 0; j < 6; j++)
     VMath::SmoothHistogram(thetas, NumOrient);

     maxval = VMath::maxVec(thetas, NumOrient);
     for(j = 0; j < NumOrient; j++)
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
            }else
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
	return newkps;
}

bool Dense::KeypointBuild(const char *fn,const char *dstfn,const char *descfn, const char* dvfn)
{
    this->crntimg  = new Image(fn);

    if(!this->crntimg->isActive())
      return false;

    AbstractDetector::releaseKpList(this->kps);
    Image *oriImg = NULL, *tmpImg = NULL;

    vector<vector<Image *> > GaussOctaves;
    vector<KeyPoint *> peaks;
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

    this->BuildOctaves(this->crntimg, _Conv2, GaussOctaves, this->sigma0);
    ///cout<<GaussOctaves.size()<<endl;

    for(ioctave = 0; ioctave < GaussOctaves.size(); ioctave++)
    {
        peaks = this->generatePeaksScales(ioctave, GaussOctaves[ioctave]);
        this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    /**/
    peaks = this->FindOrientation(kps, GaussOctaves);
    this->kps.insert(this->kps.begin(), peaks.begin(), peaks.end());
    peaks.clear();
    /**/

    Cleaner::releaseOctaves(GaussOctaves);
    if(peaks.size() > 0)
    {
        this->kps.insert(kps.begin(), peaks.begin(), peaks.end());
        peaks.clear();
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
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }

    if(strcmp(dvfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, dvfn, this->resize_rate);
    }

    delete this->crntimg;
    return true;
}

void Dense::test()
{
    const char *config = "/home/wlzhao/bin/etc/lip-vireo.conf";
    char vsrc[128], vdst[128], vdesc[128];
    AbstractDetector *mydetector = new Dense();
    mydetector->Init(config, "SIFT");

    for(int i = 1; i <= 4; i++)
    {
        sprintf(vsrc,  "/home/wlzhao/datasets/trec03/img%d.jpg",   i);
        sprintf(vdst,  "/home/wlzhao/datasets/trec03/img%d.keys",  i);
        sprintf(vdesc, "/home/wlzhao/datasets/trec03/img%d.pkeys", i);
        mydetector->KeypointBuild(vsrc, "", vdesc, ""); //vdesc
    }
}
