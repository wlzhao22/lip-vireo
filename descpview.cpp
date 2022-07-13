#include "descpview.h"

#include "viewboard.h"
#include "filter.h"

#include <fstream>
#include <cstdlib>
#include <cmath>

const int DescPView::PatchSize   = 41;
const int DescPView::PatchMag    = 20;
const int DescPView::FIFT_MAG    = 280;

DescPView::DescPView(DESC desc)
{
    this->descWin  = new float[PatchSize*PatchSize];

    cgmat = DescPView::GaussianWeight2D(PatchSize);

    this->crntImg    = NULL;
    ///this->properties[_flip_] = true;

    switch(desc)
    {
    case PVIEW:
    {
        cout<<"Descriptor .............................. un-normalized Patch view\n";
        this->featLen = 0;
        break;
    }
    case NPVIEW:
    {
        cout<<"Descriptor .............................. normalized Patch view\n";
        this->featLen = 0;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = 0;
    }
    }
    this->descoption = desc;
}

vector<vector<float> > DescPView::GaussianWeight2D(const int winSize)
{
    int   c     = winSize/2, i;
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
            v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }

    return mat;
}


int DescPView::calcDescPatch(KeyPoint *keyp, float *myWin, const int Size, const int id)
{
    assert(this->crntImg);
    int xt = 0, yt = 0;
    float sc = (2*keyp->iscale)/Size;
    float sine, cosine, gray, rPtx,rPty, ptx, pty;//rotated x and y;
    float sine0, cosine0, ux, uy;

    sine   = sin(keyp->sori);
    cosine = cos(keyp->sori);
    if(this->descoption == PVIEW)
    {
        sine0    = 0;
        cosine0  = 1;
    }else if(this->descoption == NPVIEW)
    {
        sine0   = sin(keyp->ori);
        cosine0 = cos(keyp->ori);
    }else{
        cout<<"Wrong descriptor option!\n";
        exit(0);
    }

    int   radius = (int)floor(Size/2.0);
    int   i, j, irow = 0;
    vector<float> kern;
    Filter::GaussianKernel1D(sc, AbstractDescriptor::SmthRatio, kern);

    int c = radius;
    float scx = keyp->iscale/c;
    float scy = keyp->iscale/c;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            ux   = j*cosine0 + i*sine0;
            uy   = i*cosine0 - j*sine0;
            ptx  = scx*ux; //column
            pty  = scy*uy; //row
            rPtx = (ptx*cosine + pty*sine);
            rPty = (pty*cosine - ptx*sine);

            rPtx = rPtx + keyp->x;
            rPty = rPty + keyp->y;
            xt   = (int)round(rPtx);
            yt   = (int)round(rPty);
            gray = crntImg->get2DConVal(xt, yt, kern);
            myWin[irow*Size + radius + j]  = gray;
        }
    }

    kern.clear();
    return 1;
}


int DescPView::calcPatchCurl(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sc = (2*keyp->iscale)/Size;
    int xt, yt, i, j;
    float dx, dy, gray;
    float rPtx,rPty;
    float sine, cosine;

    sine   = sin(keyp->ori);
    cosine = cos(keyp->ori);

    int radius = (int)floor(Size/2.0);
    int irow = 0;
    vector<float> kern;
    Filter::GaussianKernel1D(sc, AbstractDescriptor::SmthRatio, kern);

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
            xt   = (int)round(rPtx);
            yt   = (int)round(rPty);
            gray = crntImg->get2DConVal(xt, yt, kern);
            myWin[irow*Size + radius + j]  = gray;
        }
    }

    kern.clear();

    int pos, next_y, prev_y;
    float a, b, theta, curl =0.0f, m, angle, weight;
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
            weight = cgmat[i][j];
            weight = m*weight*DescPView::FIFT_MAG;

            angle = theta - atan2(dy, dx);

            while(angle >= PI2)
            angle = angle - PI2;

            while(angle < 0)
            angle = angle+PI2;

            if(angle < PI)
            {
                curl -= fabs(sin(angle))*weight;
            }else
            {
                curl += fabs(sin(angle))*weight;
            }
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

int DescPView::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;
    //this->properties[_div_]   = false;
    //this->properties[_flip_]  = false;

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }
        kpnumb++;
    }

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }
    outStrm.close();

    ViewBoard *myview = new ViewBoard(kpnumb, PatchSize);

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(crntKpt->KP == false)
        {
            continue;
        }

        this->calcDescPatch(crntKpt, this->descWin, PatchSize, kpnumb);
        if(this->descoption == NPVIEW)
        {
            this->calcPatchCurl(crntKpt, descWin, PatchSize);
        }
        myview->addPatch(this->descWin, PatchSize, PatchSize);
    }
    myview->saveView(descfn);
    delete myview;
    return kpnumb;
}

DescPView::~DescPView()
{
    vector<vector<float> >::iterator iter;
    for(iter = cgmat.begin(); iter != cgmat.end(); iter++)
    {
        iter->clear();
    }
    cgmat.clear();
}
