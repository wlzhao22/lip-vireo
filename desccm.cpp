#include <iostream>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "desccm.h"
#include "vmath.h"

const int DescCM::CM_NUM = 5;
const int DescCM::PatchSize = 41;

DescCM::DescCM(DESC desc)
{
    cout<<"With 'Color Moment' Decriptor.\n";
    this->descWin = new float[PatchSize*PatchSize];
    this->featLen = 0;
    CMfeats = NULL;
    block = NULL;
    descPatch = NULL;

    switch(desc)
    {
    case CM:
    {
        cout<<"Descriptor ...................... un-normalized color-moments\n";
        break;
    }
    case NCM:
    {
        cout<<"Descriptor ...................... normalized color-moments\n";
        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
    }
    }
    this->descoption = desc;

    this->bl_size  = 0;
    this->desc_patch_size = 0;
}

bool DescCM::reShape()
{
    cout<<"Input Color image is strongly recommended!\n";

    int channel = this->crntImg->channel;
    int nwfeatLen = CM_NUM*CM_NUM*3*channel;

    if(this->featLen !=  nwfeatLen || CMfeats == NULL);
    {
        this->featLen = nwfeatLen;
        CMfeats = new float[this->featLen];
    }

    int bl_height,bl_width;

    bl_width = DescCM::PatchSize/CM_NUM;
    bl_height = DescCM::PatchSize/CM_NUM;

    int nwbl_size = bl_width*bl_height*channel;

    if(nwbl_size != this->bl_size || block == NULL)
    {
        bl_size = nwbl_size;
        block = new float[bl_size];
    }

    int nwdesc_patch_size = PatchSize*PatchSize*channel;
    if(desc_patch_size != nwdesc_patch_size || this->descPatch == NULL)
    {
        desc_patch_size = nwdesc_patch_size;
        this->descPatch = new float[desc_patch_size];
    }

    return true;
}


bool DescCM::reShape(const int channel,const int width,const int height)
{
    cout<<"Input Color image is strongly recommended!\n";

    int nwfeatLen = CM_NUM*CM_NUM*3*channel;

    if(this->featLen !=  nwfeatLen || CMfeats == NULL);
    {
        this->featLen = nwfeatLen;
        CMfeats = new float[this->featLen];
    }

    int bl_height,bl_width;

    bl_width = width/CM_NUM;
    bl_height = height/CM_NUM;

    int nwbl_size = bl_width*bl_height*channel;

    if(nwbl_size != this->bl_size || block == NULL)
    {
        bl_size = nwbl_size;
        block = new float[bl_size];
    }

    int nwdesc_patch_size = PatchSize*PatchSize*channel;
    if(desc_patch_size != nwdesc_patch_size || this->descPatch == NULL)
    {
        desc_patch_size = nwdesc_patch_size;
        this->descPatch = new float[desc_patch_size];
    }

    return true;
}

int DescCM::buildDescriptor(const int kpnum,vector<KeyPoint*> &kps,const char *descfn,const float resize_rate)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntKpt;

    int kpnumb = 0;
    this->reShape();

    int nproperty = 2 + IDetector::Numb_PROP;

    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }
        kpnumb++;
    }

    ofstream outStrm(descfn);
    if(!outStrm.is_open())
    {
        cout<<"Target file '"<<descfn<<"' cannot open for write!\n";
        exit(0);
    }

    if(this->_out_format == _vgg_fmrt)
        outStrm<<this->featLen<<endl<<kpnumb<<endl;
    else if (this->_out_format == _vireo_fmrt)
        outStrm<<kpnumb<<" "<<nproperty<<" "<<this->featLen<<endl;

    kpnumb = 0;
    for(it = kps.begin(); it != kps.end(); it++)
    {
        crntKpt = *it;
        if(!crntKpt->KP)
        {
            continue;
        }

        this->getDescPatch3C(crntKpt,this->descPatch,PatchSize);
        this->getLocalColorMoment(this->descPatch);
        if(this->descoption == NCM)
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

/**
* Calculate first 3 color moment in a 3 by 3 grid for 41x41 normalized patch
*/

void DescCM::getLocalColorMoment(const float *patch)
{
    float mean[3];
    float std_dev[3];
    float skew[3];
    float Lab[3];
    float RGB[3];

    int bl_height, bl_width;
    bl_width = DescCM::PatchSize/CM_NUM;
    bl_height = DescCM::PatchSize/CM_NUM;

    int i,j,row,column;
    int counter = 0,loc = 0,size = 0;
    int sx = 0,sy = 0, ex = 0, ey = 0;
    int pix_loc,row_loc;

    memset(CMfeats,0,sizeof(float)*this->featLen);

    for(row = 0; row <CM_NUM; row++)
    {
        for(column = 0; column <CM_NUM; column++)
        {
            sx = column*bl_width;
            sy =  row*bl_height;
            ex =  (column+1)*bl_width;
            ex = ex<DescCM::PatchSize?ex:DescCM::PatchSize;
            ey = (row+1)*bl_height;
            ey = ey<DescCM::PatchSize?ey:DescCM::PatchSize;
            counter = 0;
            for(j=sy; j<ey; j++)
            {
                row_loc = j*DescCM::PatchSize;
                for(i=sx; i<ex; i++)
                {
                    pix_loc = row_loc + i;
                    RGB[0] = this->descPatch[pix_loc];
                    RGB[1] = this->descPatch[pix_loc+1];
                    RGB[2] = this->descPatch[pix_loc+2];
                    DescCM::RGB2LUV(RGB,Lab);
                    block[counter] = Lab[0];
                    block[counter+1] = Lab[1];
                    block[counter+2] = Lab[2];
                    counter = counter + 3;
                }
            }

            //calcu first 3 moments
            //calcu mean for Lab
            mean[0] = 0;//for L
            mean[1] = 0; //for alpha
            mean[2] = 0; //for belta

            for(i=0; i<counter/3; i++)
            {
                mean[0] = mean[0] + block[i];
                mean[1] = mean[1] + block[i+1];
                mean[2] = mean[2] + block[i+2];
            }

            size = counter/3;
            mean[0] = mean[0]/size;
            mean[1] = mean[1]/size;
            mean[2] = mean[2]/size;

            CMfeats[loc+0] = (float)mean[0]; //for L
            CMfeats[loc+1] = (float)mean[1]; //for alpha
            CMfeats[loc+2] = (float)mean[2]; //for belta


            std_dev[0] = 0;//for L
            std_dev[1] = 0; //for alpha
            std_dev[2] = 0; //for belta
            skew[0] = 0;//for L
            skew[1] = 0; //for alpha
            skew[2] = 0; //for belta
            //calcu std_devi and skew
            for(i=0; i<counter/3; i++)
            {
                std_dev[0] = std_dev[0] + (block[i]-mean[0])*(block[i]-mean[0]);
                std_dev[1] = std_dev[1] + (block[i+1]-mean[1])*(block[i+1]-mean[1]);
                std_dev[2] = std_dev[2] + (block[i+2]-mean[2])*(block[i+2]-mean[2]);
                skew[0] = skew[0] + pow((block[i]-mean[0]),3.0);
                skew[1] = skew[1] + pow((block[1+i]-mean[1]),3.0);
                skew[2] = skew[2] + pow((block[2+i]-mean[2]),3.0);
            }

            skew[0] = skew[0] < 0?0:skew[0];
            skew[1] = skew[1] < 0?0:skew[1];
            skew[2] = skew[2] < 0?0:skew[2];
            CMfeats[loc + 3] = (float)sqrt(std_dev[0]/size);
            CMfeats[loc + 4] = (float)sqrt(std_dev[1]/size);
            CMfeats[loc + 5] = (float)sqrt(std_dev[2]/size);
            CMfeats[loc + 6] = (float)pow(skew[0]/size,1.0f/3.0f);
            CMfeats[loc + 7] = (float)pow(skew[1]/size,1.0f/3.0f);
            CMfeats[loc + 8] = (float)pow(skew[2]/size,1.0f/3.0f);
            loc = loc + 9;
        }
    }

}


void DescCM::getLocalColorMoment(CImage *myimg,const int height,const int width)
{
    float mean[3];
    float std_dev[3];
    float skew[3];
    float Lab[3], RGB[3];

    int bl_height = 0, bl_width = 0;
    bl_width  = width/CM_NUM;
    bl_height = height/CM_NUM;

    ///cout<<width<<"\t"<<height<<endl;

    int i, j, row, column;
    int counter = 0, loc = 0, size = 0;
    int sx = 0, sy = 0, ex = 0, ey = 0;
    memset(CMfeats, 0, sizeof(float)*this->featLen);

    for(row = 0; row <CM_NUM; row++)
    {
        for(column = 0; column <CM_NUM; column++)
        {
            sx = column*bl_width;
            sy =  row*bl_height;
            ex =  (column+1)*bl_width;
            ex = ex<width?ex:width;
            ey = (row+1)*bl_height;
            ey = ey<height?ey:height;
            counter = 0;
            for(j = sy; j < ey; j++)
            {
                for(i = sx; i < ex; i++)
                {
                    myimg->getPixel(i, j, RGB);
                    DescCM::RGB2LUV(RGB, Lab);
                    block[counter]   = Lab[0];
                    block[counter+1] = Lab[1];
                    block[counter+2] = Lab[2];
                    counter = counter + 3;
                }
            }


            //calcu first 3 moments
            //calcu mean for Lab
            mean[0] = 0;//for L
            mean[1] = 0; //for alpha
            mean[2] = 0; //for belta

            for(i = 0; i < counter/3; i++)
            {
                mean[0] = mean[0] + block[i];
                mean[1] = mean[1] + block[i+1];
                mean[2] = mean[2] + block[i+2];
            }

            size = counter/3;
            mean[0] = mean[0]/size;
            mean[1] = mean[1]/size;
            mean[2] = mean[2]/size;

            CMfeats[loc+0] = (float)mean[0]; //for L
            CMfeats[loc+1] = (float)mean[1]; //for alpha
            CMfeats[loc+2] = (float)mean[2]; //for belta

            //cout<<"bug 1.02\n";

            std_dev[0] = 0;//for L
            std_dev[1] = 0; //for alpha
            std_dev[2] = 0; //for belta
            skew[0] = 0;//for L
            skew[1] = 0; //for alpha
            skew[2] = 0; //for belta
            //calcu std_devi and skew
            for(i=0; i<counter/3; i++)
            {
                std_dev[0] = std_dev[0] + (block[i]-mean[0])*(block[i]-mean[0]);
                std_dev[1] = std_dev[1] + (block[i+1]-mean[1])*(block[i+1]-mean[1]);
                std_dev[2] = std_dev[2] + (block[i+2]-mean[2])*(block[i+2]-mean[2]);
                skew[0] = skew[0] + pow((block[i]-mean[0]),3.0);
                skew[1] = skew[1] + pow((block[1+i]-mean[1]),3.0);
                skew[2] = skew[2] + pow((block[2+i]-mean[2]),3.0);
            }

            //cout<<"bug 1.04\n";

            skew[0] = skew[0] < 0?0:skew[0];
            skew[1] = skew[1] < 0?0:skew[1];
            skew[2] = skew[2] < 0?0:skew[2];
            CMfeats[loc + 3] = (float)sqrt(std_dev[0]/size);
            CMfeats[loc + 4] = (float)sqrt(std_dev[1]/size);
            CMfeats[loc + 5] = (float)sqrt(std_dev[2]/size);
            CMfeats[loc + 6] = (float)pow(skew[0]/size,1.0f/3.0f);
            CMfeats[loc + 7] = (float)pow(skew[1]/size,1.0f/3.0f);
            CMfeats[loc + 8] = (float)pow(skew[2]/size,1.0f/3.0f);
            loc = loc + 9;
        }
    }
}

void DescCM::RGB2LUV(const float RGB[],float LUV[])
{
    float L, M, S, l, m, s;
    float r = 0, g = 0, b = 0;

    r = RGB[0];  g = RGB[1];  b = RGB[2];

    l = 0.3811f*r + 0.5783f*g + 0.0402f*b;
    m = 0.1967f*r + 0.7244f*g + 0.0782f*b;
    s = 0.0241f*r + 0.1288f*g + 0.8444f*b;

    if(l <= 0)
    {
        l = 2.0f;
    }

    if(m <= 0)
    {
        m = 2.0f;
    }

    if(s <= 0)
    {
        s = 2.0f;
    }

    L = VMath::lgx(l,10);
    M = VMath::lgx(m,10);
    S = VMath::lgx(s,10);

    LUV[0] = 0.57735f*L  + 0.57735f*M   + 0.57735f*S;
    LUV[1] = 0.408248f*L + 0.408248f*M   + -0.816496f*S;
    LUV[2] = 0.707107f*L + -0.707107f*M + 0.0f*S;
    return ;
}

DescCM::~DescCM()
{
    delete [] CMfeats;
    delete [] block;
    delete [] descPatch;
}

void DescCM::test()
{
    const char *img1 = "e:/siftlab/src/trec03_img1.bmp";
    CImage *myimg = new CImage(img1);

    DescCM *mycm = new DescCM(CM);
    mycm->reShape(3,myimg->height,myimg->width);
    clock_t start = clock();
    mycm->getLocalColorMoment(myimg,myimg->height,myimg->width);
    clock_t end = clock();
    cout<<end- start<<endl;

}

