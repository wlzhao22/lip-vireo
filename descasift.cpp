#include "descasift.h"
#include "viewboard.h"
#include "vstring.h"
#include "vmath.h"

#include <fstream>
#include <cstdlib>
#include <cmath>


const int DescASIFT::NumOrient = 8;
const int DescASIFT::GRID        = 4;
const int DescASIFT::DSize       = GRID*GRID*NumOrient;
const int DescASIFT::DEGREE      = (360/NumOrient);
const int DescASIFT::PSIFTLen    = 36;
const int DescASIFT::PatchSize   = 41;
const int DescASIFT::PatchMag    = 20;

DescASIFT::DescASIFT(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];
    this->featsBin = new float[DSize];
    this->pcaFeat  = new float[PSIFTLen];

    this->crntImg    = NULL;

    gMat = DescASIFT::GaussianWeight2D(PatchSize);

    switch(desc)
    {
    case ASIFT:
    {
        cout<<"Descriptor .............................. un-normalized ASIFT\n";
        this->featLen = DSize;
        break;
    }
    case NASIFT:
    {
        cout<<"Descriptor .............................. normalized ASIFT\n";
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


vector<vector<float> > DescASIFT::GaussianWeight2D(const int winSize)
{
    int c       = winSize/2;
    float sigma = winSize/2;
    vector<vector<float> > mat;
    for (int i = 0; i < winSize; i++)
    {
        vector<float> row(winSize);
        mat.push_back(row);
    }

    float v, s2 = sigma * sigma;
    int   x, y;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }
    return mat;
}

int DescASIFT::calcDescPatch(KeyPoint *keyp, float *myWin, const int Size, const bool _NORM_)
{
    assert(this->crntImg);
    float gray = 0, rPtx = 0, rPty = 0, cosine = 0, sine = 0;//rotated x and y;
    float sine0 = 0, cosine0 = 0, ux = 0, uy = 0, ptx = 0, pty = 0;

    int   i = 0, j = 0, irow = 0;
    int   radius = (int)floor(Size/2.0);
    float scx = (keyp->iscale*keyp->e1)/(radius);
    float scy = (keyp->iscale*keyp->e2)/(radius);

    sine0   = sin(keyp->ori);
    cosine0 = cos(keyp->ori);
    ///sine0   = 0; cosine0 = 1;
    sine    = sin(keyp->sori);
    cosine  = cos(keyp->sori);

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            ux   = j*cosine0 + i*sine0;
            uy   = i*cosine0 - j*sine0;
            ptx  = scx*ux; //column
            pty  = scy*uy; //row
            rPtx = ptx*cosine + pty*sine;
            rPty = pty*cosine - ptx*sine;

            rPtx = rPtx + keyp->x;
            rPty = rPty + keyp->y;
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }
    return 1;
}

int DescASIFT::getSIFTDescriptor(const float *myWin)
{
    int i = 0, j = 0;
    for(i = 0; i < DSize; i++)
    {
        this->featsBin[i] = 0.0f;
    }

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

            mag = mag*gMat[i][j];

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

    bool GOOD = VMath::SIFTNorm(this->featsBin, DSize);

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

int DescASIFT::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
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
    float *buffer = new float[kpnumb*this->featLen];
    memset(buffer, 0, kpnumb*this->featLen*sizeof(float));

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

    float *output_feat = NULL;

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }

        ///this->getDescAffPatch(crntKpt, this->descWin, PatchSize);
        this->calcDescPatch(crntKpt, this->descWin, PatchSize, true);
        if(this->getSIFTDescriptor(this->descWin))
        {
            output_feat = this->featsBin;
            if(this->descoption == NASIFT)
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


int DescASIFT::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
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

    ViewBoard *fgrd_view = new ViewBoard(kpnumb, PatchSize);
    kpnumb = 0;

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }
        this->calcDescPatch(crntKpt, this->descWin, PatchSize, true);
        fgrd_view->addPatch(this->descWin, PatchSize, PatchSize);
        kpnumb++;
    }

    fgrd_view->saveView(dvfn);
    delete fgrd_view;

    return kpnumb;
}

DescASIFT::~DescASIFT()
{
    vector<vector<float> >::iterator iter;
    for(iter = gMat.begin(); iter != gMat.end(); iter++)
    {
        iter->clear();
    }
    gMat.clear();
    delete [] featsBin;
    delete [] pcaFeat;
}
