#include "descspin.h"
#include "viewboard.h"
#include "vmath.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

const int DescSPIN::Ints      = 10;
const int DescSPIN::Dist      = 10;
const float DescSPIN::alpha   = 1.8f;
const float DescSPIN::belta   = 3.6f;
const int DescSPIN::PatchMag  = 20;
const int DescSPIN::CLR_DEPTH = 255;
const int DescSPIN::PatchSize = 41;

DescSPIN::DescSPIN(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];
    this->featsBin = new float[DescSPIN::Ints*DescSPIN::Dist];
    this->pcaFeat = new float[DescSPIN::Ints*DescSPIN::Dist];
    this->crntImg  = NULL;

    switch(desc)
    {
    case SPIN:
    {
        cout<<"Descriptor .............................. un-normalized SPIN\n";
        this->featLen = DescSPIN::Ints*DescSPIN::Dist;
        break;
    }
    case NSPIN:
    {
        cout<<"Descriptor .............................. normalized SPIN\n";
        this->featLen = DescSPIN::Ints*DescSPIN::Dist;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = DescSPIN::Ints*DescSPIN::Dist;
    }
    }
    this->descoption = desc;
}

int DescSPIN::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    ofstream outStrm(descfn);

    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }

    int kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }
        kpnumb++;
    }

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
        if(properties[i])
            nproperty += KP_PROP_SZ[i];
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

        this->getDescPatch(crntKpt, this->descWin, PatchSize);
        this->getSPINDescriptor(this->descWin);
        if(descoption == NSPIN)
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

    return kpnumb;
}

int DescSPIN::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
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

int DescSPIN::getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sine = 0, cosine = 1, gray = 0;
    float rPtx = 0, rPty = 0, ptx = 0, pty = 0;
    if(keyp->_type == FLM)
    {
        sine   = 0;
        cosine = 1;
    }else{
        sine   = sin(keyp->sori);
        cosine = cos(keyp->sori);
    }

    int radius = (int)floor(Size/2.0);
    float sc = (2*keyp->iscale)/Size;

    int i = 0, j = 0, irow = 0;
    float scx = 0, scy = 0;
    if(keyp->_type == FLM)
    {
        scx = scy = sc;
    }else{
        scx = keyp->e1/radius;
        scy = keyp->e2/radius;
    }

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            ptx = scx*j;
            pty = scy*i;
            rPtx = (ptx*cosine + pty*sine);
            rPty = (pty*cosine - ptx*sine);
            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }

    return 1;
}


int DescSPIN::getWinSize(KeyPoint *key,float &sizeratio)
{
    float scale = 1.5*SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);
    //float scale = key->fscale;
    int patchsize = (int) (PatchMag * scale);
    // make odd
    patchsize /= 2;
    patchsize = patchsize * 2 + 1;
    sizeratio = (float) patchsize / (float) PatchSize;
    return patchsize;
}


int DescSPIN::getSPINDescriptor(Image *win, const float sizeratio)
{
    int i, j, k, m;
    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = 0.0f;
    }

    float alpha2, belta2;
    int radius = PatchSize/2;
    float intense, _intense, dist, _dist, step_dist, step_ints;
    unsigned int  _idx;

    step_dist = PatchMag/DescSPIN::Dist;
    step_ints = CLR_DEPTH/DescSPIN::Ints;
    alpha2 = -2.0f*DescSPIN::alpha*DescSPIN::alpha;
    belta2 = -2.0f*DescSPIN::belta*DescSPIN::belta;
    for(i = 0; i < PatchSize; i++)
    {
        for(j = 0; j< PatchSize; j++)
        {
            intense = win->getPixel(i,j);
            dist = (i - radius)*(i - radius)  + (j - radius)*(j - radius);
            dist = sqrt(dist);
            for(m = 0; m < DescSPIN::Dist; m++)
            {
                _dist = dist - m*step_dist;
                _dist = (_dist*_dist)/alpha2;
                for(k = 0; k < DescSPIN::Ints; k++)
                {
                    _intense = intense - k*step_ints;
                    _intense = (_intense*_intense)/belta2;
                    _idx = m*DescSPIN::Dist + k;
                    this->featsBin[_idx] =  exp(_intense+_dist);
                }
            }
        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = round(this->featsBin[i]);
    }
    return 0;
}

int DescSPIN::getSPINDescriptor(const float *myWin)
{
    int i, j, k, m;
    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = 0.0f;
    }

    float alpha2, belta2;
    int pos1, pos2, radius = PatchSize/2;
    float  intense, _intense, dist, _dist, step_dist, step_ints;
    unsigned int  _idx;

    step_dist = PatchSize/DescSPIN::Dist;
    step_ints = DescSPIN::CLR_DEPTH/DescSPIN::Ints;
    alpha2 = -2.0f*DescSPIN::alpha*DescSPIN::alpha;
    belta2 = -2.0f*DescSPIN::belta*DescSPIN::belta;

    for(i = 0; i < PatchSize; i++)
    {
        pos1 = i*PatchSize;
        for(j = 0; j< PatchSize; j++)
        {
            pos2 = pos1 + j;
            intense = myWin[pos2];
            dist = (i - radius)*(i - radius)  + (j - radius)*(j - radius);
            dist = sqrt(dist);
            for(k = 0; k < DescSPIN::Ints; k++)
            {
                _intense = intense - k*step_ints;
                _intense = (_intense*_intense)/belta2;
                for(m = 0; m < DescSPIN::Dist; m++)
                {
                    _dist = dist - m*step_dist;
                    _dist = (_dist*_dist)/alpha2;
                    _idx = k*DescSPIN::Ints + m;
                    this->featsBin[_idx] +=  exp(_intense+_dist);
                }
            }
        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = round(this->featsBin[i]);
    }

    return 0;
}


DescSPIN::~DescSPIN()
{

    delete [] featsBin;
    delete [] pcaFeat;
}
