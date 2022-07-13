#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

#include "viewboard.h"
#include "descfift.h"
#include "cleaner.h"
#include "vmath.h"

using namespace std;

const int DescFIFT::NumOrient   = 8;
const int DescFIFT::GRID        = 4;
const int DescFIFT::DSize       = GRID*GRID*NumOrient;
const int DescFIFT::PSIFTLen    = 36;
const int DescFIFT::PatchSize   = 41;
const int DescFIFT::PatchMag    = 20;

DescFIFT::DescFIFT(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];
    this->featsBin = new float[DSize];
    this->pcaFeat = new float[PSIFTLen];

    this->crntImg    = NULL;
    ///this->properties[_flip_] = true;

    ///gMat  = DescFIFT::GaussianWeight2D(PatchSize);
    cgmat = DescFIFT::GaussianWeight2D(PatchSize);

    switch(desc)
    {
    case FIFT:
    {
        cout<<"Descriptor .............................. un-normalized F-SIFT\n";
        this->featLen = DSize;
        break;
    }
    case NFIFT:
    {
        cout<<"Descriptor .............................. normalized F-SIFT\n";
        this->featLen = DSize;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = DSize;
    }
    }
    this->descoption = desc;
}

vector<vector<float> > DescFIFT::GaussianWeight2D(const int winSize)
{
    int   c     = winSize/2, i = 0;
    int   dim   = c*2 + 1;
    float sigma = c/2;
    vector<vector<float> > mat;
    for (i = 0; i < dim; i++)
    {
        vector<float> row(dim);
        mat.push_back(row);
    }

    float s2 = sigma * sigma;
    float v;
    int   x, y;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }
    VMath::normMat(mat);
    return mat;
}

vector<vector<float> > DescFIFT::GaussianWeight2D(const float sigma, const int radius)
{
    assert (sigma > 0);
    int c = radius;
    int dim = c*2 + 1;
    vector<vector<float> > mat;
    for (int i = 0; i < dim; i++)
    {
        vector<float> row(dim);
        mat.push_back(row);
    }

    float v = 0, s2 = sigma * sigma;
    int x = 0, y = 0;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }

    Filter::normalizeMat(mat);
    return mat;
}

int DescFIFT::buildDescriptor(const int kpnum0, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;
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
    float *buffer = new float[kpnumb*this->featLen];
    memset(buffer, 0, kpnumb*this->featLen*sizeof(float));

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }

    float *output_feat = NULL;

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }

        this->calcPatchCurl(crntKpt, this->descWin, PatchSize);

        if(crntKpt->flip == 0)
        {
            crntKpt->KP = false;
            continue;
        }

        this->calcDescPatch(crntKpt, this->descWin, PatchSize);

        if(this->getFIFTDescriptor(this->descWin, crntKpt->iscale))
        {
            output_feat = this->featsBin;
            if(this->descoption == NFIFT)
            {
                VMath::l2norm(output_feat, this->featLen);
            }
            memcpy(buffer+kpnumb*this->featLen, this->featsBin, sizeof(float)*this->featLen);
            kpnumb++;
        }
        else
        {
            crntKpt->KP = false;
        }
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
        output_feat = buffer+kpnumb*this->featLen;
        if(this->_out_format == _vgg_fmrt)
        {
            this->saveDescVGG(crntKpt, output_feat, this->featLen, resize_rate, outStrm);
        }
        else if(this->_out_format == _vireo_fmrt)
        {
            this->saveDescVireo(crntKpt, output_feat, kpnumb, this->featLen, resize_rate, outStrm);
        }
        kpnumb++;
    }

    outStrm.close();
    delete [] buffer;
    buffer = NULL;
    return kpnumb;
}

int DescFIFT::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
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
        this->calcPatchCurl(crntKpt, this->descWin, PatchSize);

        if(crntKpt->flip == 0)
        {
            crntKpt->KP =false;
            continue;
        }

        this->calcDescPatch(crntKpt, this->descWin, PatchSize);
        myview->addPatch(this->descWin, PatchSize, PatchSize);
    }
    myview->saveView(dvfn);
    delete myview;
    return 0;
}

int DescFIFT::getWinSize(KeyPoint *key,float &sizeratio)
{
    float scale = 1.5*SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);
    int patchsize = (int) (PatchMag * scale);
    patchsize /= 2;
    patchsize = patchsize * 2 + 1;

    sizeratio = (float) patchsize / (float) PatchSize;

    return patchsize;
}

int DescFIFT::calcPatchCurl(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sc = (2*keyp->iscale)/Size;
    int i, j;
    float dx, dy, gray;
    float rPtx,rPty;
    float sine, cosine;

    sine   = sin(keyp->ori);
    cosine = cos(keyp->ori);

    int radius = (int)floor(Size/2.0);
    int irow = 0;
    sine   = sine*sc;
    cosine = cosine*sc;

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

    int pos, next_y, prev_y;
    float a, b, theta, curl = 0.0f, m, angle, weight;
    float Radius = PatchSize/2;

    for(i = 1; i < PatchSize -1; i++)
    {
        irow = i*PatchSize;
        for(j = 1; j < PatchSize-1; j++)
        {
            pos = irow + j;
            a = myWin[ pos + 1] - myWin[pos - 1]; //dx
            next_y = irow + PatchSize;
            prev_y = irow - PatchSize;
            b = myWin[next_y + j] - myWin[prev_y + j];//dy
            m = sqrt(a*a +  b*b);
            theta = atan2(b,a);
            dx = j - Radius;
            dy = i - Radius;

            if(dx == 0.0 && dy == 0.0)
            {
                continue;
            }

            /***theoritically speaking, weighting below makes more sense*/
            weight = m*cgmat[i][j];
            angle  = theta - atan2(dy, dx);

            while(angle >= PI2)
                angle = angle - PI2;

            while(angle < 0)
                angle = angle+PI2;

            /**/
            if(angle < PI)
            {
                curl -= fabs(sin(angle))*weight;
            }
            else
            {
                curl += fabs(sin(angle))*weight;
            }
            /**/
        }
    }

    if(curl > 0)
    {
        keyp->flip = 1;
    }
    else
    {
        keyp->flip = -1;
    }

    return keyp->flip;
}

int DescFIFT::calcDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sine = 0.0f, cosine = 0.0f, gray = 0.0f, rPtx = 0.0f, rPty = 0.0f;
    float sc = (2*keyp->iscale)/Size;

    if(keyp->flip == -1)
    {
        if(keyp->ori > 0)
        {
            sine      = sin(PI - keyp->ori);
            cosine    = cos(PI - keyp->ori);
            keyp->ori = PI - keyp->ori;
        }
        else
        {
            sine      = -1*sin(PI + keyp->ori);
            cosine    = cos(PI + keyp->ori);
            keyp->ori = -1*(PI + keyp->ori);
        }
    }
    else
    {
        sine   = sin(keyp->ori);
        cosine = cos(keyp->ori);
    }

    int   radius = (int)floor(Size/2.0);
    int   i = 0, j = 0, irow = 0;
    short flip = keyp->flip ;
    sine   = sine*sc;
    cosine = cosine*sc;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            rPtx = j*cosine + i*sine;
            rPty = i*cosine - j*sine;

            rPtx = flip*rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j] = gray;
        }
    }

    return 1;
}

int DescFIFT::getFIFTDescriptor(const float *myWin, const float sc)
{
    int i = 0, j = 0;
    memset(this->featsBin, 0, sizeof(float)*DSize);

    float angle = 0, mag = 0;
    float histWidth = PatchSize/(GRID + 0.0f);
    float Theta_per_bin = NumOrient/PI2;
    float rbin, cbin, obin;
    float rdr, rdc, rda, a, b;
    int rbin0, cbin0, obin0, ir, ic, ia;
    int rloc, cloc, oloc, row, bloc;
    int crnt, prev, next;
    float dx, dy, da;

    for(i = 1; i < PatchSize-1; i++)
    {
        rbin = i/histWidth;
        rbin = rbin - 0.5;
        rbin0 = (int)floor(rbin);
        dy = rbin - rbin0;

        for(j = 1; j < PatchSize-1; j++)
        {
            cbin = j/histWidth;
            cbin = cbin - 0.5;
            cbin0 = (int)floor(cbin);
            dx = cbin - cbin0;

            crnt = i*PatchSize;
            a = myWin[crnt + j + 1] - myWin[crnt + j - 1]; //dx
            next = crnt + PatchSize;
            prev = crnt - PatchSize;
            b = myWin[next + j] - myWin[prev +j];//dy
            mag = sqrt(a*a +  b*b);
            angle = atan2(b, a);
            angle = angle < 0? (angle + PI2):angle;

            obin  = angle*Theta_per_bin;
            obin0 = (int) floor(obin);
            da = obin - obin0;

            mag = mag*cgmat[i][j];

            for(ir = 0; ir <= 1; ir++)
            {
                rloc = rbin0 + ir;

                if(rloc >= GRID || rloc < 0)
                    continue;

                rdr = mag*((ir == 0)?(1 - dy):dy);
                row = rloc*GRID;

                for(ic = 0; ic <= 1; ic++)
                {
                    cloc = cbin0 + ic;

                    if(cloc >= GRID || cloc < 0)
                        continue;

                    rdc  = rdr*((ic == 0)?(1 - dx):dx);
                    bloc = (row + cloc)*NumOrient;

                    for(ia = 0; ia <= 1; ia++)
                    {
                        oloc = (obin0 + ia)%NumOrient;
                        rda = rdc*((ia == 0)?(1 - da):da);
                        this->featsBin[bloc+oloc] += rda;
                    }
                }
            }//end-for, end of interpolation
        }//end inner for
    }//end outer for

    ///bool GOOD = VMath::sqrtSIFTNorm(this->featsBin, DSize);
    bool GOOD = VMath::SIFTNorm(this->featsBin, DSize);
    Cleaner::clear2DArray(gMat);

    for(i = 0; i < DSize; i++)
    {
        this->featsBin[i] = round(this->featsBin[i]);
    }

    if(GOOD)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

DescFIFT::~DescFIFT()
{
    vector<vector<float> >::iterator iter;
    for(iter = gMat.begin(); iter != gMat.end(); iter++)
    {
        iter->clear();
    }
    gMat.clear();
    for(iter = cgmat.begin(); iter != cgmat.end(); iter++)
    {
        iter->clear();
    }
    cgmat.clear();
    delete [] featsBin;
    delete [] pcaFeat;
}
