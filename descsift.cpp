#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>

#include "viewboard.h"
#include "descsift.h"
#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

using namespace std;

const int DescSIFT::NumOrient = 8;
const int DescSIFT::GRID      = 4;
const int DescSIFT::DSize     = GRID*GRID*NumOrient;
const int DescSIFT::PSIFTLen  = 36;
const int DescSIFT::PatchSize = 41;
const int DescSIFT::PatchMag  = 20;

/**
According to my experiment on Oxford5k with VLAD and Covariant-VLAD,
performing square-root on raw SIFT (rrSIFT) before normalization performs
much better than performing square-root on normalized SIFT (rnsift). The re-
sults are follows.

OX5K	    	Feat.	mAP
HessrrSIFT  	vlad*   0.333
HessrrSIFT  	cvlad   0.392
HessrrFIFT  	vlad*   0.321
HessrrFIFT  	cvlad   0.395
HesaffrnSIFT	vlad*   0.171
HesaffrnSIFT	cvlad   0.244
HessrnSIFT  	vlad*   0.207
HessrnSIFT  	cvlad	0.269
DSURFrrSIFT     vlad*   0.365
DSURFrrSIFT     cvlad   0.465

Hesaff=Hessian-Affine, which is from VGG and INRIA
Hess=Hessian, which from this LIP-VIREO package
FIFT=Flip invariant SIFT

@date:      Apr.-4-2014
@author:    Wan-Lei Zhao
**/

DescSIFT::DescSIFT(DESC desc)
{
    this->descWin = new float[PatchSize*PatchSize];
    this->featsBin = new float[DSize];
    this->pcaFeat = new float[PSIFTLen];

    gMat = DescSIFT::GaussianWeight2D(PatchSize);

    this->crntImg = NULL;
    switch(desc)
    {
    case NSIFT:
    {
        cout<<"Descriptor .............................. normalized SIFT\n";
        this->featLen = DSize;
        break;

    }
    case SIFT:
    {
        cout<<"Descriptor .............................. unnormalized SIFT\n";
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

vector<vector<float> > DescSIFT::GaussianWeight2D(const float sigma, const int radius)
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

vector<vector<float> > DescSIFT::GaussianWeight2D(const int winSize)
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
    ///VMath::normMat(mat);
    return mat;
}

int DescSIFT::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;
    ofstream outStrm(descfn);
    assert(this->crntImg);
    float *tmpFeat = new float[this->featLen];

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
    float *buffer = new float[kpnumb*this->featLen];
    memset(buffer, 0, kpnumb*this->featLen*sizeof(float));

    int i = 0, nproperty = 2;
    for(i = 0; i < IDetector::Numb_PROP; i++)
    {
        if(properties[i])
        {
            nproperty += KP_PROP_SZ[i];
        }
    }

    float *output_feat = NULL;

    assert(this->featsBin);

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        assert(crntKpt);
        if(!crntKpt->KP)
        {
            continue;
        }

        this->getDescPatch(crntKpt, this->descWin, PatchSize);

        if(this->getSIFTDescriptor(this->descWin, crntKpt->iscale) == 1)
        {
            /**it works, but not necessary*/
            ///VMath::sift2csift(this->featsBin, this->featLen, tmpFeat);
            ///memcpy(this->featsBin, tmpFeat, this->featLen*sizeof(float));
            if(this->descoption == NSIFT)
            {
                output_feat = this->featsBin;
                VMath::l2norm(output_feat, this->featLen);
            }
            else if(this->descoption == SIFT)
            {
                output_feat = this->featsBin;
            }
            if(output_feat != NULL)
            {
                memcpy(buffer+kpnumb*this->featLen, this->featsBin, sizeof(float)*this->featLen);
                kpnumb++;
            }
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
        ///outStrm<<kpnumb<<" "<<4<<" "<<32<<endl;
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
    delete [] tmpFeat;
    tmpFeat = NULL;

    return kpnumb;
}

int DescSIFT::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
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
    myview->saveView(descfn);
    delete myview;
    return 0;
}

int DescSIFT::getNormDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    int xt, yt, i, j;
    int radius = (int)floor(Size/2.0);
    float sine   = sin(keyp->ori);
    float cosine = cos(keyp->ori);
    float sc = keyp->iscale/radius;
    float lx, ly, dx, dy, gray;
    sine   = sc*sine;
    cosine = sc*cosine;
    float ptx, pty, rptx, rpty, tempval;
    vector<float> kern = Filter::GaussianKernel1D(sc);

    int irow = 0;
    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            ptx = keyp->e1*j; ///column
            pty = keyp->e2*i; ///row
            rptx = ptx*cosine - pty*sine;
            rpty = ptx*sine + pty*cosine;
            rptx = rptx + keyp->x; ///column
            rpty = rpty + keyp->y;  ///row
            lx =  floor(rptx);  ///column
            ly =  floor(rpty);  ///row
            xt =  (int)round(rptx);
            yt =  (int)round(rpty);
            dx =  rptx - lx;
            dy =  rpty - ly;

            tempval = crntImg->get2DConVal(xt, yt, kern);
            gray    = tempval*(1-dx)*(1-dy);
            tempval = crntImg->get2DConVal(xt+1, yt, kern);
            gray   += tempval*(1-dy)*dx;
            tempval = crntImg->get2DConVal(xt, yt+1, kern);
            gray   += tempval*(1-dx)*dy;
            tempval = crntImg->get2DConVal(xt+1, yt+1, kern);
            gray    = tempval*dy*dx;
            myWin[irow*Size + radius + j]  = gray;
        }
    }
    kern.clear();
    return 1;
}

bool DescSIFT::getPixGradient(const int x, const int y, const float *pixels,
                              const int size, float &mag, float &theta)
{
    assert(pixels);

    if (x < 1 || y < 1 || x > size - 2 || y > size - 2)
        return false;

    int posi1 = y*size;
    float a = pixels[posi1+x + 1] - pixels[posi1+x - 1]; //dx
    int posi2 = posi1+size;
    posi1 =posi1 - size;
    float b = pixels[posi2+x] - pixels[posi1 +x];//dy
    mag = sqrt(a*a +  b*b);
    theta = atan2(b, a);
    return true;
}

bool DescSIFT::getPixGradient(const int x, const int y, Image *win,
                              const float sizeratio, float &mag, float &theta)
{
    assert(win);
    float x1, x2, y1, y2, gx, gy;
    int size = win->width;
    if (x < 1 || y < 1 || x > size - 2 || y > size - 2)
        return false;

    x1 = win->getPixelBI((float) (x+1) * sizeratio, (float) y * sizeratio);
    x2 = win->getPixelBI((float) (x-1) * sizeratio, (float) y * sizeratio);
    y1 = win->getPixelBI((float) x * sizeratio, (float) (y + 1) * sizeratio);
    y2 = win->getPixelBI((float) x * sizeratio, (float) (y - 1) * sizeratio);
    gx = x1 - x2;
    gy = y1 - y2;

    mag = gx*gx + gy*gy;
    mag = sqrt(mag);
    theta = atan2(gy, gx);
    return true;
}


int DescSIFT::getSIFTDescriptor(const float *myWin, const float sc)
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
                        this->featsBin[bloc+oloc] = this->featsBin[bloc+oloc] + rda;
                    }
                }
            }//end-for, end of interpolation
        }//end inner for
    }//end outer for

    ///bool GOOD = VMath::sqrtSIFTNorm(this->featsBin, DSize);
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
        return -1;
    }
}

void DescSIFT::getPatch(const char*srcimg, const char *srcfn, const char *dstpatch)
{
    Image* myimg = new Image(srcimg);
    Image* patchImg = new Image(81, 81);
    ifstream *inStrm = new ifstream(srcfn);
    int num;
    int x, y, xi, yi, x_flr, y_flr;
    float v, a, b, c, sc, score;
    float tx, ty, dx, dy, gray;
    int counter = 0;

    (*inStrm)>>v;
    (*inStrm)>>num;
    cout<<v<<"\t"<<num<<endl;
    while(!inStrm->eof())
    {
        if(counter < v)
        {
            (*inStrm)>>x;
            (*inStrm)>>y;
            (*inStrm)>>a;
            (*inStrm)>>b;
            (*inStrm)>>c;
            (*inStrm)>>sc;
            (*inStrm)>>score;
            if(x == 613 && y == 89)
            {
                cout<<"selected \n";
                break;
            }
            counter++;
        }
        else
        {
            break;
        }
    }
    inStrm->close();

    sc = sc/40.5f;
    float m[2][2];

    cout<<a<<"\t"<<b<<endl;
    cout<<b<<"\t"<<c<<endl;
    cout<<x<<"\t"<<y<<endl;
    cout<<sc<<endl;
    //a = 1; b = 0; c = 1; sc = 1;
    m[0][0] = sc*a;
    m[0][1] = sc*b;
    m[1][0] = sc*b;
    m[1][1] = sc*c;

    for(yi = -40; yi < 41; yi++)
    {
        for(xi = -40; xi < 41; xi++)
        {
            tx = m[0][0]*xi + m[0][1]*yi;
            ty = m[1][0]*xi + m[1][1]*yi;
            tx = tx + x;
            ty = ty + y;
            x_flr = floor(tx);
            y_flr = floor(ty);
            dx = tx - x_flr;
            dy = ty - y_flr;
            gray = myimg->getPixel(x_flr, y_flr)*(1-dx)*(1-dy);
            gray += myimg->getPixel(x_flr+1,y_flr)*dx*(1-dy);
            gray += myimg->getPixel(x_flr, y_flr+1)*(1-dx)*dy;
            gray += myimg->getPixel(x_flr+1,y_flr+1)*dx*dy;
            patchImg->setPixel(xi+40, yi+40, gray);
        }
        //cout<<endl;
    }

    patchImg->setPixel(40, 40, 254);

    patchImg->save(dstpatch);

}

DescSIFT::~DescSIFT()
{
    vector<vector<float> >::iterator iter;
    for(iter = gMat.begin(); iter != gMat.end(); iter++)
    {
        iter->clear();
    }
    gMat.clear();

    delete [] featsBin;
    featsBin = NULL;
    delete [] pcaFeat;
    pcaFeat = NULL;
}

void DescSIFT::test()
{
    const char* srcfn    = "e:/wlzhao/matlab/mser/graf_img1.keys";
    const char* srcimg   = "e:/wlzhao/matlab/mser/graf_img1.jpg";
    const char* dstpatch = "e:/wlzhao/matlab/mser/patch_img1_8.pgm";
    DescSIFT::getPatch(srcimg, srcfn, dstpatch);
}
