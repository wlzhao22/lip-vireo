#include "abstractdescriptor.h"
#include "viewboard.h"
#include "filter.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

const float AbstractDescriptor::SIGMA         = 1.6f;
const float AbstractDescriptor::SmthRatio     = 0.35f; //3.0f
const int AbstractDescriptor::SCALESPEROCTAVE = 3;

bool AbstractDescriptor::setupImage(AbstractImage *srcImg)
{
    assert(srcImg);
    this->crntImg = srcImg;

    return true;
}


bool AbstractDescriptor::setupImage(CImage *srcImg)
{
    assert(srcImg);
    this->crntImg = srcImg;
    return true;
}

bool AbstractDescriptor::setupOctaves(vector<vector<Image*> > &ImgOctaves)
{
    ImageOctaves = ImgOctaves;
    return true;
}

int AbstractDescriptor::findFlip(const float thetas[], const int dbin, const int nbin, float &ratio)
{
    double lefta  = 0.0f;
    double righta = 0.0f;
    int idx = 0;
    int num = ceil((nbin+0.0f)/4.0);
    for(int k = dbin+1, i= 0; i < num; k++, i++)
    {
        idx     = k%nbin;
        righta += thetas[idx];
    }
    for(int k = dbin-1, i= 0; i < num; k--, i++)
    {
        idx    = (k+nbin)%nbin;
        assert(idx < nbin && idx >=0);
        lefta += thetas[idx];
    }

    ratio = righta > lefta? (lefta/righta):(righta/lefta);
    if(righta > lefta)
    {
        return  1;
    }
    else
    {
        return -1;
    }
}

int AbstractDescriptor::getCurl(KeyPoint *keyp,float *myWin,const int Size)
{
    assert(this->crntImg);
    float sc = (2*keyp->iscale)/Size;
    int xt, yt, i, j;

    float sine, cosine;
    if(this->descoption == NSPIN || this->descoption == SPIN)
    {
        sine = 0;
        cosine = 1;
    }
    else
    {
        sine = sin(keyp->ori);
        cosine = cos(keyp->ori);
    }
    float dx, dy, gray, rPtx,rPty;

    int radius = (int)floor(Size/2.0);
    int irow = 0;
    vector<float> kern;
    Filter::GaussianKernel1D(sc, AbstractDescriptor::SmthRatio, kern);
    float u[2][2];

    u[0][0] = (keyp->a*cosine - keyp->b*sine)*sc;
    u[0][1] = (keyp->b*cosine - keyp->c*sine)*sc;
    u[1][0] = (keyp->a*sine + keyp->b*cosine)*sc;
    u[1][1] = (keyp->b*sine + keyp->c*cosine)*sc;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            rPtx = j*u[0][0] + i*u[0][1];
            rPty = j*u[1][0] + i*u[1][1];
            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row

            /***/
            xt =  (int)round(rPtx);
            yt =  (int)round(rPty);
            gray = crntImg->get2DConVal(xt, yt, kern);
            myWin[irow*Size + radius + j]  = gray;
            /***/

            /**
            gray = crntImg->getPixelBI(rPtx,rPty);
            **/
        }
    }

    kern.clear();

    int pos, next_y, prev_y;
    float a, b, theta, curl = 0.0f, m, angle, weight;
    float Radius = Size/2;

    for(i = 1; i < Size -1; i++)
    {
        irow = i*Size;
        for(j = 1; j < Size-1; j++)
        {
            pos = irow + j;
            a = myWin[ pos + 1] - myWin[pos - 1]; //dx
            next_y = irow + Size;
            prev_y = irow - Size;
            b = myWin[next_y + j] - myWin[prev_y + j];//dy
            m = sqrt(a*a +  b*b);
            theta = atan2(b,a);
            dx = j - Radius;
            dy = i - Radius;
            if(myWin[pos] != 0)
            {
                weight = m/myWin[pos];
            }
            else
            {
                weight = m;
            }

            angle = theta - atan2(dy, dx);

            if(angle >=0)
            {
                if(angle <=PI)
                {
                    curl -= fabs(sin(angle))*weight;
                }
                else
                {
                    curl += fabs(sin(angle))*weight;
                }

            }
            else
            {
                if(angle >= -PI )
                {
                    curl += fabs(sin(angle))*weight;
                }
                else
                {
                    curl -= fabs(sin(angle))*weight;
                }
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


int AbstractDescriptor::getDescPatch1(KeyPoint *keyp, float *myWin, const int Size)
{
    int i = 0, j = 0, irow = 0;
    float sine = 0.0f, cosine = 0.0f, cx, cy;
    if(this->descoption == NSPIN || this->descoption == SPIN)
    {
        sine   = 0;
        cosine = 1;
    }
    else
    {
        sine = sin(keyp->ori);
        cosine = cos(keyp->ori);
    }
    float gray = 0.0f, rPtx = 0.0f, rPty = 0.0f;
    int radius = (int)floor(Size/2.0f);
    float sc = keyp->iscale/radius;
    sine    = sine*sc;
    cosine  = cosine*sc;
    float U[2][2] = {0};
    U[0][0] = keyp->a; U[0][1] = keyp->b;
    U[1][0] = keyp->b; U[1][1] = keyp->c;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {

            cx = j*U[0][0] + i*U[0][1];
            cy = j*U[1][0] + i*U[1][1];

            rPtx = cx*cosine - cy*sine;
            rPty = cy*cosine + cx*sine;

            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }
    //kern.clear();
    return 1;
}

int AbstractDescriptor::getDescPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    float gray = 0.0f, rPtx = 0.0f, rPty = 0.0f;
    int radius = (int)floor(Size/2.0f);
    float sc = keyp->iscale/radius, cx, cy;
    int i = 0, j = 0, irow = 0;
    float sine = 0.0f, cosine = 0.0f;
    float U[2][2] = {0};
    if(this->descoption == NSPIN || this->descoption == SPIN)
    {
        sine   = 0;
        cosine = 1;
    }
    else
    {
        sine = sin(keyp->ori);
        cosine = cos(keyp->ori);
    }

    sine    = sine*sc;
    cosine  = cosine*sc;

    U[0][0] = keyp->a; U[0][1] = keyp->b;
    U[1][0] = keyp->b; U[1][1] = keyp->c;

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            cx = j*U[0][0] + i*U[0][1];
            cy = j*U[1][0] + i*U[1][1];

            //rPtx = j*cosine + i*sine;
            //rPty = i*cosine - j*sine;
            rPtx = cx*cosine + cy*sine;
            rPty = cy*cosine - cx*sine;

            rPtx = rPtx + keyp->x; //column
            rPty = rPty + keyp->y;  //row
            gray = crntImg->getPixelBI(rPtx, rPty);
            myWin[irow*Size + radius + j]  = gray;
        }
    }
    //kern.clear();
    return 1;
}


int AbstractDescriptor::getDescAffPatch(KeyPoint *keyp, float *myWin, const int Size)
{
    assert(this->crntImg);
    float sine, cosine, gray, rPtx, rPty, ptx, pty;//rotated x and y;
    float sine0, cosine0, ux, uy;

    sine   = sin(keyp->sori);
    cosine = cos(keyp->sori);

    if(this->descoption == NSPIN || this->descoption == SPIN)
    {
        sine0  = 0;
        cosine0 = 1;
    }
    else
    {
        sine0   = sin(keyp->ori);
        cosine0 = cos(keyp->ori);
    }

    int   i = 0, j = 0, irow = 0;
    int   radius = (int)floor(Size/2.0);
    float scx = (keyp->iscale*keyp->e1)/radius;
    float scy = (keyp->iscale*keyp->e2)/radius;

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

int AbstractDescriptor::getDescPatch3C(KeyPoint *keyp,float *myWin,const int Size)
{
    assert(this->crntImg);
    float sc = (2*keyp->iscale)/Size;
    int xt,yt;
    //keyp->ori = 0;
    float sine = sin(keyp->ori);
    float cosine = cos(keyp->ori);
    float rPtx,rPty;
    float pixel[3];

    int radius = (int)floor(Size/2.0);
    int i,j;
    int irow = 0;
    vector<float> kern;
    Filter::GaussianKernel1D(sc, AbstractDescriptor::SmthRatio, kern);

    for(i = -radius; i <= radius; i++)
    {
        irow = i + radius;
        for(j = -radius; j <= radius; j++)
        {
            rPtx = j*cosine-i*sine;
            rPty = j*sine+i*cosine;
            rPtx = rPtx*sc + keyp->x; //column
            rPty = rPty*sc + keyp->y;  //row

            xt =  (int)round(rPtx);
            yt =  (int)round(rPty);

            crntImg->get2DConVal(xt, yt, pixel, kern);//*(1-dx)*(1-dy);
            myWin[irow+ radius + j]    = pixel[0];
            myWin[irow+ radius + j+1]  = pixel[1];
            myWin[irow+ radius + j +2] = pixel[2];
        }
    }
    kern.clear();

    return 0;
}

void AbstractDescriptor::saveDescVGG(const KeyPoint* crnt_kpt, const float *feature, const int dim,
                                     const float resize, ofstream &outStrm)
{
    if(!outStrm.is_open())
    {
        cout<<"File is not open!\n";
        return ;
    }
    assert(feature);
    float a, b, c;
    float mat[2][2];
    float e1 = 1.0f/crnt_kpt->dscale;
    e1       = e1*e1;
    outStrm<<crnt_kpt->x*resize<<" "<<crnt_kpt->y*resize<<" ";

    if(this->DETECTOR == corner)
    {
        mat[0][0] = cos(crnt_kpt->sori);
        mat[0][1] = sin(crnt_kpt->sori);
        mat[1][0] = -1*sin(crnt_kpt->sori);
        mat[1][1] = cos(crnt_kpt->sori);
        a         = mat[0][0]*mat[0][0]*crnt_kpt->e1 + mat[0][1]*mat[0][1]*crnt_kpt->e2;
        b         = mat[0][0]*mat[0][1]*crnt_kpt->e1 - mat[0][1]*mat[1][0]*crnt_kpt->e2;
        c         = mat[1][0]*mat[1][0]*crnt_kpt->e1 + mat[0][0]*mat[0][0]*crnt_kpt->e2;
        outStrm<<a<<" "<<b<<" "<<c<<endl;
    }
    else
    {
        outStrm<<e1<<" "<<0<<" "<<e1;
    }
    for(int i = 0; i < dim; i++)
    {
        outStrm<<" "<<feature[i];
    }
    outStrm<<endl;

    return ;
}

void AbstractDescriptor::saveDescTrain(const KeyPoint* crnt_kpt, const float *feature, const int dim,
                                       const float resize, ofstream &outStrm)
{
    if(!outStrm.is_open())
    {
        cout<<"File is not open!\n";
        return ;
    }
    assert(feature);
    int i;
    for(i = 0; i < dim-1; i++)
    {
        outStrm<<feature[i]<<" ";
    }
    outStrm<<feature[i]<<endl;

    return ;
}

void AbstractDescriptor::saveDescVireo(const KeyPoint* crnt_kpt, const float *feature, const int idx0,
                                       const int dim, const float resize, ofstream &outStrm)
{
    if(!outStrm.is_open())
    {
        cout<<"File is not open!\n";
        return ;
    }
    assert(feature);

    int i = 0, tmpx = (int)round(crnt_kpt->x*resize);
    int tmpy = (int)round(crnt_kpt->y*resize);
    outStrm<<tmpx<<" "<<tmpy;
    outStrm<<" "<<crnt_kpt->gscale;
    outStrm<<" "<<crnt_kpt->ori<<endl;
    /**
    if(properties[_flip_])
    {
        outStrm<<" "<<crnt_kpt->flip;
    }else{
        outStrm<<" 0";
    }
    **/

    /**/
    for(i = 1; i <= dim; i++)
    {
        if(i%12)
        {
            outStrm<<setprecision(4)<<feature[i-1]<<" ";
        }
        else
        {
            outStrm<<setprecision(4)<<feature[i-1]<<endl;
        }
    }
    if(dim%12)
    {
        outStrm<<endl;
    }
    /**/
    return ;
}
