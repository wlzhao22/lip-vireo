#include "descfind.h"

#include "viewboard.h"
#include "filter.h"
#include "vmath.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

const int DescFIND::ORIENTATION = 8;
const int DescFIND::GRID        = 4;
const int DescFIND::DSize       = GRID*GRID*ORIENTATION;
const int DescFIND::DEGREE      = (360/ORIENTATION);
const int DescFIND::PSIFTLen    = 36;
const int DescFIND::PatchSize   = 41;
const int DescFIND::PatchMag    = 20;
const int DescFIND::SIFT_MAG    = 280;
const float DescFIND::inv_sqrt2 = 0.707106781f;

DescFIND::DescFIND(DESC desc)
{
    this->descWin = new float[PatchSize*PatchSize];
    this->featsBin = new float[81];
    this->pcaFeat = new float[81];

    gmat = DescFIND::GaussianWeight2D(PatchSize);

    this->crntImg = NULL;
    switch(desc)
    {
    case NFIND:
    {
        cout<<"Descriptor .............................. normalized FIND\n";
        this->featLen = 81;
        break;

    }
    case FIND:
    {
        cout<<"Descriptor .............................. unnormalized FIND\n";
        this->featLen = 81;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = 81;
    }
    }
    this->descoption = desc;
}

vector<vector<float> > DescFIND::GaussianWeight2D(const int winSize)
{
    int c       = winSize/2;
    float sigma = winSize/6;

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
            v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }
    return mat;
}

int DescFIND::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
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

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
         if(properties[i])
         nproperty += KP_PROP_SZ[i];
    }
    float *output_feat = NULL;

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

        if(this->getFINDDescriptor(this->descWin, crntKpt->flip) == 1)
        {
            if(this->descoption == NFIND)
            {
                output_feat = this->featsBin;
                VMath::l2norm(output_feat, this->featLen);
            }
            else if(this->descoption == FIND)
            {
                output_feat = this->featsBin;
            }
            if(output_feat != NULL)
            {
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
        }
        else
        {
            crntKpt->KP = false;
        }
        /**/
    }

    outStrm.close();
    return kpnumb;
}


int DescFIND::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *dvfn, const float resize_rate)
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
        this->getDescPatch(crntKpt, this->descWin, PatchSize);
        myview->addPatch(this->descWin, PatchSize, PatchSize);
    }
    myview->saveView(dvfn);
    delete myview;
    return 0;
}

int DescFIND::getFINDDescriptor(const float *myWin, const int flip)
{
    int i, j, idxi, idxj, si, sj;

    for(i = 0; i < this->featLen; i++)
    {
        this->featsBin[i] = 0.0f;
    }

    float angle, mag, weight;
    float radius = PatchSize/2;
    radius = radius*radius;
    float Theta_per_bin = ORIENTATION/PI2;
    float obin, d_o;
    int   obin0;
    int   oloc[2];
    int   bloc;
    int   crnt, prev, next;
    float a, b, wghs[2];

    for(i = -19; i < 19; i++)
    {
        idxi  = i + PatchMag;
        crnt  = idxi*PatchSize;
        si    = abs(i);

        for(j = -19; j < 19; j++)
        {
            idxj = j + PatchMag;
            sj   = abs(j);

            a     = myWin[crnt + idxj + 1] - myWin[crnt + idxj - 1]; //dx
            next  = crnt + PatchSize;
            prev  = crnt - PatchSize;
            b     = myWin[next + idxj] - myWin[prev + idxj];//dy

            mag   = sqrt(a*a +  b*b);
            angle = atan2(b, a);

            while(angle < 0.0)
                angle += PI2;

            while(angle >= PI2)
                angle -= PI2;

            obin    = angle*Theta_per_bin;
            obin0   = (int)floor(obin);
            d_o     = obin - obin0;

            weight  = gmat[idxi][idxj];
            mag     = mag*weight*DescFIND::SIFT_MAG;
            wghs[0] = mag*(1-d_o);
            wghs[1] = mag*d_o;
            oloc[0] = (obin0)%ORIENTATION;
            oloc[1] = (obin0+1)%ORIENTATION;

            if(i <= 0 && j <= 0)
            {
                //bin2
                bloc = 1*(ORIENTATION+1);
                this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
            }
            else if(i > 0 && j <= 0)
            {
                //bin8
                bloc = 7*(ORIENTATION+1);
                this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
            }
            else if(i <= 0 && j > 0)
            {
                //bin1
                bloc = 0*(ORIENTATION+1);
                this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
            }
            else
            {
                //bin9
                bloc = 8*(ORIENTATION+1);
                this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
            }

            if(sj <= 5 && sj <= 5)
            {
                if(i <= 0 && j <= 0)
                {
                    //bin3
                    bloc = 2*(ORIENTATION+1);
                    this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                    this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
                }
                else if(i > 0 && j <= 0)
                {
                    //bin7
                    bloc = 6*(ORIENTATION+1);
                    this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                    this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
                }
                else if(i <= 0 && j > 0)
                {
                    //bin4
                    bloc = 3*(ORIENTATION+1);
                    this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                    this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
                }
                else if(i > 0 && j > 0)
                {
                    //bin6
                    bloc = 5*(ORIENTATION+1);
                    this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                    this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
                }
            }

            if(si <= 2 && sj <= 2)
            {
                //bin4
                bloc = 4*(ORIENTATION+1);
                this->featsBin[bloc+oloc[0]] = this->featsBin[bloc+oloc[0]] + wghs[0];
                this->featsBin[bloc+oloc[1]] = this->featsBin[bloc+oloc[1]] + wghs[1];
            }
            /**/
        }//end inner for
    }//end outer for

    /**/
    for(i = 0; i < this->featLen; i += 9)
    {
        this->featsBin[i+ORIENTATION] = this->featsBin[i]*DescFIND::inv_sqrt2;
        this->featsBin[i]             = this->featsBin[i+ORIENTATION];
    }

    if(flip == -1)
    {
        for(i = 0; i < this->featLen; i++)
        {
            this->pcaFeat[this->featLen-i-1] = this->featsBin[i];
        }
        memcpy(this->featsBin, this->pcaFeat, sizeof(float)*this->featLen);
    }
    /**/
    return 1;
}

DescFIND::~DescFIND()
{
    vector<vector<float> >::iterator iter;
    for(iter = gmat.begin(); iter != gmat.end(); iter++)
    {
        iter->clear();
    }
    gmat.clear();


    delete [] featsBin;
    delete [] pcaFeat;
}

void DescFIND::test()
{
}

