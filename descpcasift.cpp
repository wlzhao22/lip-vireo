#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "descpcasift.h"
#include "viewboard.h"
#include "vstring.h"
#include "filter.h"
#include "vmath.h"

using namespace std;

const unsigned int DescPCASIFT::PCALen   = 36;
const unsigned int DescPCASIFT::PatchMag = 20;

DescPCASIFT::DescPCASIFT(DESC desc, const char *pcamapfn)
{
    this->featsBin  = NULL;
    this->pcaFeat  = new float[PCALen];
    this->gradients = new float[GPLEN];
    this->descWin   = new float[PatchSize*PatchSize];
    this->crntImg   = NULL;
    strcpy(pcafn,   pcamapfn);

    switch(desc)
    {
    case NPCASIFT:
    {
        cout<<"Descriptor .............................. normalized PCA-SIFT\n";
        this->featLen = PCALen;
        computeDesc = &DescPCASIFT::getPCASIFTDesc;
        this->featsBin = this->pcaFeat;
        this->paramsCheck();
        break;
    }
    case PCASIFT:
    {
        cout<<"Descriptor .............................. unnormalized PCA-SIFT\n";
        this->featLen = PCALen;
        computeDesc = &DescPCASIFT::getPCASIFTDesc;
        this->featsBin = this->pcaFeat;
        this->paramsCheck();
        break;
    }
    case PCAPATCH:
    {
        cout<<"Descriptor .............................. unnormalized Patches for PCA-SIFT\n";
        this->featLen  = GPLEN;
        computeDesc    = &DescPCASIFT::makePatches;
        this->featsBin = this->gradients;
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
        this->featLen = PCALen;
    }
    }
    this->descoption = desc;
}

void DescPCASIFT::loadPCAMatrix(const char *pcafn)
{
    assert(pcafn);
    FILE * fp = fopen(pcafn, "rb");
    float val;

    if(fp == NULL)
    {
        cout<<"File "<<pcafn<<" can not open\n";
        exit(1);
    }

    for (int i = 0; i < GPLEN; i++)
    {
        if (fscanf(fp, "%f", &val) != 1)
        {
            printf("Invalid pca vector file (avg).\n");
            exit(1);
        }
        avgs[i] = val;
    }

    /**** read in vector, transposed*/
    for (int i = 0; i < GPLEN; i++)
    {
        for (unsigned int j = 0; j < PCALen; j++)
        {
            if (fscanf(fp, "%f", &val) != 1)
            {
                printf("Invalid pca vector file (eig).\n");
                exit(1);
            }
            if (j < EPCALEN)
                eigs[j][i] = val;
        }
    }
    fclose(fp);
}

bool DescPCASIFT::paramsCheck()
{
    cout<<"Loading PCA mapping matrix ... ";
    this->loadPCAMatrix(pcafn);
    cout<<"done.\n";
    return true;
}

int DescPCASIFT::buildDescriptor(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;
    unsigned int kpnumb = 0;
    ofstream outStrm;
    char dstfn[1024];
    if(this->_out_format == _train_fmrt)
    {
        strcpy(dstfn, "");
        VString::parseDIR(dstfn, descfn);
        strcat(dstfn, "gpatches.txt");
        outStrm.open(dstfn, ios::app);
    }
    else
    outStrm.open(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(1);
    }

    int nproperty = 2 + IDetector::Numb_PROP;

    kpnumb = 0;
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

        this->getDescPatch(crntKpt, this->descWin, PatchSize);

        (*this.*computeDesc)(this->descWin);

        if(this->descoption == NPCASIFT)
        {
            featsBin = this->pcaFeat;
            VMath::l2norm(featsBin, this->featLen);
        }

        if(this->_out_format == _vgg_fmrt)
        {
            this->saveDescVGG(crntKpt, this->featsBin, this->featLen, resize_rate, outStrm);
        }
        else if(this->_out_format == _vireo_fmrt)
        {
            this->saveDescVireo(crntKpt, this->featsBin, kpnumb, this->featLen, resize_rate, outStrm);
        }else if(this->_out_format == _train_fmrt)
        {
            this->saveDescTrain(crntKpt, this->featsBin, this->featLen, resize_rate, outStrm);
        }
        kpnumb++;
    }

    outStrm.close();
    return kpnumb;
}

int DescPCASIFT::buildPatchView(const int kpnum, vector<KeyPoint*> &kps, const char *descfn, const float resize_rate)
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
    myview->saveView(descfn);
    delete myview;
    return kpnumb;
}

int DescPCASIFT::getNormDescPatch(KeyPoint *keyp,float *myWin,const int Size)
{
    assert(this->crntImg);
    int radius = (int)floor(Size/2.0f);
    float sc = keyp->iscale/radius;
    float gray;
    int i, j, irow, xt, yt;
    float rPtx,rPty;
    float sine, cosine;
    sine = sin(keyp->ori);
    cosine = cos(keyp->ori);
    sine    = sine*sc;
    cosine  = cosine*sc;
    vector<float> kern;
    Filter::GaussianKernel1D(sc,AbstractDescriptor::SmthRatio, kern);

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
            //gray = crntImg->getPixelBI(rPtx, rPty);
            //myWin[irow*Size + radius + j]  = gray;
        }
    }
    kern.clear();
    return 1;
}

int DescPCASIFT::getWinSize(KeyPoint *key,float &sizeratio)
{
    float scale = 1.5*SIGMA * pow(2.0, key->fscale / (float) SCALESPEROCTAVE);
    int patchsize = (int) (PatchMag * scale);
    patchsize /= 2;
    patchsize = patchsize * 2 + 1;

    sizeratio = (float) patchsize / (float) PatchSize;

    return patchsize;
}

int DescPCASIFT::makePatches(const float *myWin)
{
    unsigned int x, y, count = 0;
    unsigned int row, brow, arow;
    float x1,  x2, y1, y2, gx, gy;

    memset(this->gradients, 0, this->GPLEN*sizeof(float));

    //cout<<"i am here\n";

    for (y = 1; y < PatchSize - 1; y++)
    {
        row  = PatchSize*y;
        brow = PatchSize*(y-1);
        arow = PatchSize*(y+1);

        for (x = 1; x < PatchSize - 1; x++)
        {
            x1 = myWin[row+x+1];
            x2 = myWin[row+x-1];

            y1 = myWin[arow+x];
            y2 = myWin[brow+x];

            gx = x1 - x2;
            gy = y1 - y2;
            gradients[count] = gx;
            gradients[count + 1] = gy;
            count += 2;
        }
    }
    VMath::l1norm(gradients, GPLEN, 100);
    return 1;
}

int DescPCASIFT::getPCASIFTDesc(const float *myWin)
{
    int x, y, count = 0;
    unsigned int row, brow, arow;
    float x1, x2, y1, y2, gx, gy;

    memset(this->gradients, 0, this->GPLEN*sizeof(float));

    for (y = 1; y < PatchSize - 1; y++)
    {
        row  = PatchSize*y;
        brow = PatchSize*(y-1);
        arow = PatchSize*(y+1);

        for (x = 1; x < PatchSize - 1; x++)
        {
            x1 = myWin[row+x+1];
            x2 = myWin[row+x-1];

            y1 = myWin[arow+x];
            y2 = myWin[brow+x];

            gx = x1 - x2;
            gy = y1 - y2;
            gradients[count] = gx;
            gradients[count + 1] = gy;
            count += 2;
        }
    }

    VMath::l1norm(gradients, GPLEN, 100);

    for(x = 0; x < GPLEN; x++)
    {
        gradients[x] = gradients[x] - avgs[x];
    }

    for(y = 0; y < EPCALEN; y++)
    {
        this->pcaFeat[y] = 0;

        for(x = 0; x < GPLEN; x++)
        {
            this->pcaFeat[y] += eigs[y][x] * gradients[x];
        }
    }

    VMath::l2norm(this->pcaFeat, EPCALEN);

    return 1;
}

void DescPCASIFT::test()
{
}

