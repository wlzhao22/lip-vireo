#include "descrift.h"

#include "index_template.h"
#include "viewboard.h"
#include "vmath.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

const int DescRIFT::Ori = 8;
const int DescRIFT::BLCK = 4;
const int DescRIFT::Dist = 8; //6;
const int DescRIFT::DistF = 4; //6;
const int DescRIFT::PatchMag = 20;
const int DescRIFT::factor = 2;
const int DescRIFT::PatchSize = 41;
const float DescRIFT::Blck_Size = PI2/DescRIFT::BLCK;
const float DescRIFT::Theta_per_bin = PI2/DescRIFT::Ori;
const float DescRIFT::Dist_per_bin = (PatchMag+1)/DescRIFT::Dist;
const float DescRIFT::alpha = 1.0f;

const float Dist_Map[9]  = {0.0f, 7.78f, 11.0f, 13.40f, 15.559f, 17.395f, 19.05f, 21.06f, 27.0f};
const float RDist_Map[5] = {0.0f, 5.0f, 10.0f, 15.0f, 20.0f};

/**********************Log for Color-SIFT**************************************
I tried color sift ( by Jan-Mark Geusebroek, in "Performance
evaluation of local colour invariants", CVIU'09). However failed. The performance
 is slightly worse than the current version. The reason behind this is that it
 becomes grey-invariant which reduces the distinctiveness of local patch.

date: 26/Jan/2010
***************************************************************************/

DescRIFT::DescRIFT(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];

    switch(desc)
    {
    case ERIFT:
    {
        cout<<"Descriptor .............................. un-normalized ERIFT\n";
        this->featLen = DescRIFT::factor *DescRIFT::Ori*DescRIFT::Dist;
        computeDesc = &DescRIFT::getERIFTDescriptor;
        break;
    }
    case NERIFT:
    {
        cout<<"Descriptor .............................. normalized ERIFT\n";
        this->featLen = DescRIFT::factor *DescRIFT::Ori*DescRIFT::Dist;
        computeDesc = & DescRIFT::getERIFTDescriptor;
        break;
    }
    case RIFT:
    {
        cout<<"Descriptor .............................. un-normalized RIFT\n";
        this->featLen = DescRIFT::Ori*DescRIFT::DistF;
        computeDesc = & DescRIFT::getRIFTDescriptor;
        break;
    }case NRIFT:
    {
        cout<<"Descriptor .............................. normalized RIFT\n";
        this->featLen = DescRIFT::Ori*DescRIFT::DistF;
        computeDesc = & DescRIFT::getRIFTDescriptor;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = DescRIFT::factor *DescRIFT::Ori*DescRIFT::Dist;
    }
    }
    this->descoption = desc;

    this->featsBin = new float[this->featLen];
    this->pcaFeat = new float[this->featLen];
    this->crntImg = NULL;

    /**initialize nn map**/
}


int DescRIFT::buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    int kpnumb = 0;

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
         if(properties[i])
         nproperty += KP_PROP_SZ[i];
    }

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }
        kpnumb++;
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

        this->getNormDescPatch(crntKpt, this->descWin, PatchSize);
        (*this.*computeDesc)(this->descWin);

        if(this->descoption == NERIFT || this->descoption == NRIFT)
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

int DescRIFT::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
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

int DescRIFT::getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size)
{
    assert(this->crntImg);
    float sine = 0, cosine = 1;
    float gray;
    float rPtx, rPty;

    int radius = (int)floor(Size/2.0);
    float sc   = (2*keyp->iscale)/Size;
    int i, j, irow = 0;
    ///vector<float> kern;
    ///Filter::GaussianKernel1D(sc,AbstractDescriptor::SmthRatio, kern);

    cosine = cosine*sc;
    sine   = sine*sc;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            rPtx = j*cosine + i*sine;
            rPty = i*cosine - j*sine;
            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }

    ///kern.clear();
    return 1;
}


int DescRIFT::getWinSize(KeyPoint *key,float &sizeratio)
{
    float scale = 1.5*SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);
    int patchsize = (int) (PatchMag * scale);
    // make odd
    patchsize /= 2;
    patchsize = patchsize * 2 + 1;
    sizeratio = (float) patchsize / (float) PatchSize;
    return patchsize;
}

int DescRIFT::getRIFTDescriptor(const float *myWin)
{
    int i, j;
    float a, b, mag;
    float dx, dy, angle, angle1, angle2, a_r, dist;
    int pos, row, next_y, prev_y, radius = PatchMag;
    unsigned int idx, idx_d, idx_a;

    memset(this->featsBin, 0, sizeof(float)*this->featLen);

    for(i = 1; i < PatchSize-1; i++)
    {
        row = i*PatchSize;
        for(j = 1; j < PatchSize-1; j++)
        {
            pos    = row + j;
            next_y = row + PatchSize;
            prev_y = row - PatchSize;
            a      = myWin[ pos + 1]   - myWin[pos - 1]; ///dx
            b      = myWin[next_y + j] - myWin[prev_y + j];///dy

            mag    = sqrt(a*a +  b*b);
            angle1 = atan2(b, a);

            dx     = j - radius;
            dy     = i - radius;
            idx_d  = idx_table[pos];
            dist   = sqrt(dx*dx+dy*dy);
            if(dist <= 5)
            {
                idx_d = 0;
            }else if(dist <= 10)
            {
                idx_d = 1;
            }else if(dist < 15)
            {
                idx_d = 2;
            }else if(dist <= 20)
            {
                idx_d = 3;
            }else{
                continue;
            }

            angle2 = atan2(dy, dx);
            angle = angle1 - angle2;

            while(angle < 0.0)
                angle += PI2;

            while(angle >= PI2)
                angle -= PI2;

            angle = angle/DescRIFT::Theta_per_bin;
            idx_a = (int)floor(angle);
            a_r   = angle - idx_a;
            idx_a = idx_a < 0?(DescRIFT::Ori-1):idx_a;

            idx   = idx_d*DescRIFT::Ori + idx_a;

            this->featsBin[idx] += mag*(1 - a_r);
            idx = idx_d*DescRIFT::Ori + (idx_a + 1)%DescRIFT::Ori;
            this->featsBin[idx] += mag*a_r;

        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = round(this->featsBin[i]);
    }
    return 1;
}


/**this works, don't touch it*/
int DescRIFT::getERIFTDescriptor(const float *myWin)
{
    int i, j, k;
    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = 0.0f;
    }

    float a, b, mag, mag_sin, mag_cos;
    int pos, row, next_y, prev_y, radius = PatchMag;
    float dx, dy;
    float angle, angle1, angle2, angle3;
    float a_r, r1;
    unsigned int  _idx, idx_d, idx_a;
    int next_dim = DescRIFT::Ori*DescRIFT::Dist;

    for(i = 1; i < PatchSize-1; i++)
    {
        row = i*PatchSize;
        for(j = 1; j < PatchSize-1; j++)
        {
            pos    = row + j;
            next_y = row + PatchSize;
            prev_y = row - PatchSize;
            a      = myWin[ pos + 1] - myWin[pos - 1]; ///dx
            b      = myWin[next_y + j] - myWin[prev_y + j];///dy

            mag    = sqrt(a*a +  b*b)/4;
            angle1 = atan2(b, a) + PI2;
            angle3 = atan2(b, -1*a) + PI2;
            mag_sin = fabs(sin(angle1)*mag);
            mag_cos = fabs(cos(angle1)*mag);
            dx     = j - radius;
            dy     = i - radius;
            idx_d  = idx_table[pos];

            angle2 = atan2(dy, dx) - PI/2+PI2;
            angle  = angle1 - angle2;

            while(angle < 0.0)
                angle += PI2;

            while(angle >= PI2)
                angle -= PI2;

            angle = angle/DescRIFT::Theta_per_bin;
            idx_a = (int)floor(angle);
            a_r = angle - idx_a;

            for(k = 0; k <= 1; k++)
            {
                r1 = mag_sin*((k ==0)?(1-a_r):a_r);
                _idx = idx_d*DescRIFT::Ori + (k+idx_a);
                this->featsBin[_idx] += r1;
            }

            for(k = 0; k <= 1; k++)
            {
                r1 = mag_cos*((k ==0)?(1-a_r):a_r);
                _idx = next_dim + idx_d*DescRIFT::Ori + (k+idx_a);
                this->featsBin[_idx] += r1;
            }
            angle = (angle3 - angle2);

            while(angle < 0.0)
                angle += PI2;

            while(angle >= PI2)
                angle -= PI2;

            angle = angle/DescRIFT::Theta_per_bin;
            idx_a = (int)floor(angle);
            a_r = angle - idx_a;

            for(k = 0; k <= 1; k++)
            {
                r1 = mag_sin*((k ==0)?(1-a_r):a_r);
                _idx = idx_d*DescRIFT::Ori + (k+idx_a);
                this->featsBin[_idx] += r1;
            }

            for(k = 0; k <= 1; k++)
            {
                r1 = mag_cos*((k ==0)?(1-a_r):a_r);
                _idx = next_dim + idx_d*DescRIFT::Ori + (k+idx_a);
                this->featsBin[_idx] += r1;
            }
        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = round(this->featsBin[i]);
    }

    return 0;
}


DescRIFT::~DescRIFT()
{

    delete [] featsBin;
    delete [] pcaFeat;
}
