#include "abstractdetector.h"
#include "descpcasift.h"
#include "descasift.h"
#include "descsift.h"
#include "descljet.h"
#include "descspin.h"
#include "descrift.h"
#include "descfift.h"
///#include "descgloh.h"
#include "descsurf.h"
#include "descfind.h"
#include "desccm.h"

#include "scriptparser.h"
#include "kpdrawer.h"
#include "vstring.h"
#include "cleaner.h"
#include "filter.h"
#include "config.h"
#include "iotool.h"
#include "kernel.h"
#include "prjmat.h"
#include "vmath.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cmath>

const int AbstractDetector::SCALESPEROCTAVE = 3;
const int AbstractDetector::MaxOctaves = 5;
const int AbstractDetector::BORDER     = 4;
const int AbstractDetector::KNN        = 9;
const float AbstractDetector::minDist  = 3.0*1.415;
const int AbstractDetector::EPS        = 0;
const float AbstractDetector::Ec       = 0.90f;
const float AbstractDetector::El       = 6.0f;
const float AbstractDetector::alpha    = 0.06f;
const float AbstractDetector::deltaSc0 = 1.2f;
const float AbstractDetector::deltaSc1 = 0.1f;
const int AbstractDetector::kp_num0    = 750;
const int AbstractDetector::imgSzBound = 1500;
const int AbstractDetector::nOrient0   = 32;
const int AbstractDetector::nScale0    = 16;
const float AbstractDetector::PI0      = 3.1415926;
const float AbstractDetector::PI4      = 12.5663706;
const int AbstractDetector::KN0        = 8;
const int AbstractDetector::Ki0        = 8;


AbstractDetector::AbstractDetector()
{
    INIT = false;
    this->myDescriptor = NULL;
    this->crntimg   = NULL;

    this->sel_option = TOPK;

    for(int i = 0; i < 5; i++)
    {
        this->kp_property[i] = false;
    }

    this->AFF_OUT     = false;
    this->CIRCLE_OUT  = false;
    this->out_FMRT    = _vireo_fmrt;
    this->fix_kp_numb = AbstractDetector::kp_num0;

    strcpy(detectorsID[0],  "_hessian");
    strcpy(detectorsID[1],  "_dog");
    strcpy(detectorsID[2],  "_harris");
    strcpy(detectorsID[3],  "_non");
    strcpy(detectorsID[4],  "_log");
    strcpy(detectorsID[5],  "_harlap");
    strcpy(detectorsID[6],  "_hesslap");
    strcpy(detectorsID[7],  "_mser");
    strcpy(detectorsID[8],  "_dense");
    strcpy(detectorsID[9], "_dsurf");
    strcpy(detectorsID[10], "_msurf");
    strcpy(detectorsID[11], "_corner");
    strcpy(detectorsID[12], "_hesaff");
    strcpy(detectorsID[13], "_hesaff");

}

bool AbstractDetector::Init(const char *scriptfn, const char *descOpt)
{
    char tmp[20];
    if(!strcmp(scriptfn, ""))
    {
        cout<<"Configure file is required (-c conf)!\n";
        exit(0);
    }

    ScriptParser::getDtTask(this->paras, scriptfn);
    this->sel_option = THRSH;

    this->paramsCheck();

    this->AFF_OUT = false;
    if(this->paras.find("affine") != this->paras.end())
    {
        strcpy(tmp, paras["affine"]);
        VString::toLower(tmp);
        if(!strcmp(tmp, "yes"))
        {
            if(this->DETECTOR == hessian || this->DETECTOR == harlap || this->DETECTOR == hesslap ||
                    this->DETECTOR == _log)
            {
                cout<<"Affine adaption ......................... yes\n";
                this->AFF_OUT= true;
            }
        }
    }

    this->RESIZE = false;
    this->resize_rate = 1.0f;

    if(!strcmp(descOpt, "SIFT"))
    {
        this->myDescriptor = new DescSIFT(SIFT);
        this->mydesc = SIFT;
    }
    else if(!strcmp(descOpt, "NSIFT"))
    {
        this->myDescriptor = new DescSIFT(NSIFT);
        this->mydesc = NSIFT;
    }
    else if(!strcmp(descOpt, "NLJET"))
    {
        this->myDescriptor = new DescLJet(NLJET);
        this->mydesc = NLJET;
    }
    else if(!strcmp(descOpt, "LJET"))
    {
        this->myDescriptor = new DescLJet(LJET);
        this->mydesc = LJET;
    }
    else if(!strcmp(descOpt, "PCAPATCH"))
    {
        this->myDescriptor = new DescPCASIFT(PCAPATCH, "");
        this->mydesc = PCAPATCH;
    }
    else if(!strcmp(descOpt, "SPIN"))
    {
        this->myDescriptor = new DescSPIN(SPIN);
        this->mydesc = SPIN;
    }
    else if(!strcmp(descOpt, "NSPIN"))
    {
        this->myDescriptor = new DescSPIN(NSPIN);
        this->mydesc = NSPIN;
    }
    else if(!strcmp(descOpt, "ERIFT"))
    {
        this->myDescriptor = new DescRIFT(ERIFT);
        this->mydesc = ERIFT;
    }
    else if(!strcmp(descOpt, "NERIFT"))
    {
        this->myDescriptor = new DescRIFT(NERIFT);
        this->mydesc = NERIFT;
    }
    else if(!strcmp(descOpt, "RIFT"))
    {
        this->myDescriptor = new DescRIFT(RIFT);
        this->mydesc = NERIFT;
    }
    else if(!strcmp(descOpt, "NRIFT"))
    {
        this->myDescriptor = new DescRIFT(NRIFT);
        this->mydesc = NERIFT;
    }
    else if(!strcmp(descOpt, "FIFT"))
    {
        this->myDescriptor = new DescFIFT(FIFT);
        this->mydesc = FIFT;
    }
    else if(!strcmp(descOpt, "NFIFT"))
    {
        this->myDescriptor = new DescFIFT(NFIFT);
        this->mydesc = NFIFT;
    }
    else if(!strcmp(descOpt, "SURF"))
    {
        this->myDescriptor = new DescSURF(SURF);
        this->mydesc = SURF;
    }
    else if(!strcmp(descOpt, "NSURF"))
    {
        this->myDescriptor = new DescSURF(NSURF);
        this->mydesc = NSURF;
    }
    else if(!strcmp(descOpt, "ESURF"))
    {
        this->myDescriptor = new DescSURF(ESURF);
        this->mydesc = ESURF;
    }
    else if(!strcmp(descOpt, "NESURF"))
    {
        this->myDescriptor = new DescSURF(NESURF);
        this->mydesc = NESURF;
    }
    else if(!strcmp(descOpt, "AOD"))
    {
        this->myDescriptor = new DescSURF(AOD);
        this->mydesc = AOD;
    }
    else if(!strcmp(descOpt, "NAOD"))
    {
        this->myDescriptor = new DescSURF(NAOD);
        this->mydesc = NAOD;
    }
    else if(!strcmp(descOpt, "FIND"))
    {
        this->myDescriptor = new DescFIND(FIND);
        this->mydesc = FIND;
    }
    else if(!strcmp(descOpt, "NFIND"))
    {
        this->myDescriptor = new DescFIND(NFIND);
        this->mydesc = NFIND;
    }
    else if(!strcmp(descOpt, "ASIFT"))
    {
        this->myDescriptor = new DescASIFT(ASIFT);
        this->mydesc = ASIFT;
    }
    else if(!strcmp(descOpt, "NASIFT"))
    {
        this->myDescriptor = new DescASIFT(NASIFT);
        this->mydesc = NASIFT;
    }
    else
    {
        if(!strcmp(descOpt, ""))
        {
            cout<<"No descriptor selected!\n";
            cout<<descOpt<<endl;
        }
        else
        {
            cout<<"Wrong Option of descriptor '"<<descOpt<<"'\n";
        }
        this->myDescriptor = NULL;
    }

    cout<<"Initilaizing ............................ start";

    for(int i = 0; i < Numb_PROP; i++)
    {
        this->kp_property[i] = false;
    }

    this->kp_property[_scale_] = true;
    this->kp_property[_angle_] = true;
    this->kp_property[_flip_]  = false;

    this->CIRCLE_OUT = true;
    if(this->paras.find("circle") != this->paras.end())
    {
        strcpy(tmp,paras["circle"]);
        VString::toLower(tmp);
        if(!strcmp(tmp,"yes"))
        {
            this->CIRCLE_OUT = true;
        }
        else
        {
            this->CIRCLE_OUT = false;
        }
    }

    this->TIMEON = false;
    if(paras.find("log") != paras.end())
    {
        strcpy(this->tm_logfn, paras["log"]);
        if(VString::validatePath(this->tm_logfn))
        {
            cout<<"\nLog on processing ....................... on\n";
            this->TIMEON = true;
        }
        else
        {
            cout<<"\nLog on processing ....................... off\n";
            cout<<"Path of log file '"<<this->tm_logfn<<"' is invalid!\n";
            this->TIMEON = false;
            strcpy(this->tm_logfn, "");
        }
    }
    else
    {
        cout<<"\nLog on processing ....................... off\n";
        strcpy(this->tm_logfn, "");
        this->TIMEON = false;
    }

    if(this->paras.find("format") != this->paras.end())
    {
        strcpy(tmp,paras["format"]);
        VString::toLower(tmp);
        if(!strcmp(tmp, "vgg"))
        {
            this->out_FMRT = _vgg_fmrt;
            cout<<"Descriptor output ....................... format = vgg\n";
        }
        else if(!strcmp(tmp, "train"))
        {
            this->out_FMRT = _train_fmrt;
            cout<<"Descriptor output ....................... format = train\n";
        }
        else
        {
            cout<<"Descriptor output ....................... format = vireo\n";
        }
    }
    else
    {
        cout<<"Descriptor output ....................... format = vireo\n";
    }

    if(this->out_FMRT == _vgg_fmrt)
    {
        saveKpts = &AbstractDetector::saveKpVGG;
    }
    else
    {
        saveKpts = &AbstractDetector::saveKpVIREO;
    }

    if(this->myDescriptor != NULL)
    {
        this->myDescriptor->setOutFormat(this->kp_property, this->AFF_OUT, this->out_FMRT);
    }
    INIT = true;
    cout<<"Initializing ............................ done.\n";

    Cleaner::clearParas(paras);

    return INIT;
}

void AbstractDetector::BuildOctaves(Image * image, GBlur blur_opt, vector<vector<Image *> > &Ggoctaves,
        const float _SIGMA0, const int NUM_SCALES0)
{
    assert(image);
    int dim        = min(image->height, image->width);
    int numoctaves = int (log((double) dim) / log(2.0)) - 2;
    numoctaves     = numoctaves>MaxOctaves?MaxOctaves:numoctaves;
    Image *timage  = image->clone();
    Image *simage  = NULL;
    vector<Image *> imgScales;
    int i = 0;

    for (i = 0; i < numoctaves; i++)
    {
        imgScales = BuildScales(timage, blur_opt, _SIGMA0, NUM_SCALES0);
        Ggoctaves.push_back(imgScales);
        Filter::BlurImage(timage, _SIGMA0);
        simage = Image::halfSizeImage(timage);
        delete timage;
        timage = simage;
    }
    delete timage;
    timage = NULL;
}

vector<Image*> AbstractDetector::BuildScales(Image * image, GBlur blur_opt,
                                        const float _SIGMA0, const int NUM_SCALES0)
{
    assert(image);
    vector<Image*> Scales;

    double k = pow(2, 1.0/(float)NUM_SCALES0);
    float sigma1, sigma2, sigma;
    Image *dggImg = NULL;
    Image *tmpImg = NULL;
    int   Scales_Num;

    Scales_Num = NUM_SCALES0 + 2;
    vector<float> Gkern;
    vector<float> dkern;

    if(blur_opt == _Dxy_)
    {
        tmpImg = new Image(image->width, image->height);
    }

    for (int i =  1; i <= Scales_Num; i++)
    {
        dggImg = new Image(image->width, image->height);
        sigma1 = pow(k, i - 1) *_SIGMA0;
        sigma2 = pow(k, i) *_SIGMA0;

        switch(blur_opt)
        {
        case _Conv2:
        {
            sigma = sqrt(sigma2*sigma2 - sigma1*sigma1);
            if(Scales.size() == 0)
            {
                Filter::BlurImage(image, dggImg, sigma);
            }
            else
            {
                Filter::BlurImage(Scales[Scales.size()-1], dggImg, sigma);
            }
            break;
        }
        case _Dx_:
        {
            dkern = Filter::GaussianDxKernel1D(sigma2);
            Filter::Convolve1DWidth(dkern, image, dggImg);
            dkern.clear();
            break;
        }
        case _Dy_:
        {
            dkern = Filter::GaussianDxKernel1D(sigma2);
            Filter::Convolve1DHeight(dkern, image, dggImg);
            dkern.clear();
            break;
        }
        case _Dxx_:
        {
            dkern = Filter::GaussianDxxKernel1D(sigma2);
            Filter::Convolve1DWidth(dkern, image, dggImg);
            dkern.clear();
            break;
        }
        case _Dyy_:
        {
            dkern = Filter::GaussianDxxKernel1D(sigma2);
            Filter::Convolve1DHeight(dkern,image,dggImg);

            dkern.clear();
            break;
        }
        case _Dxy_:
        {
            dkern = Filter::GaussianDxKernel1D(sigma2);
            Filter::Convolve1DWidth(dkern, image, tmpImg);
            Filter::Convolve1DHeight(dkern, tmpImg, dggImg);
            dkern.clear();

            break;
        }
        }
        Scales.push_back(dggImg);
    }

    if(blur_opt == _Dxy_)
    {
        delete tmpImg;
    }

    return Scales;
}


bool AbstractDetector::saveKpVIREO(vector<KeyPoint*> &kplst, const int kpnum, const char *dstfn,
                                   const float resize_rate0, bool AFF_OUT)
{
    FILE *fp = fopen(dstfn, "w");
    if(fp == NULL)
    {
        cout<<"File '"<<dstfn<<"' cannot open!\n";
        return false;
    }

    vector<KeyPoint*>::iterator it;
    KeyPoint *crntPt;
    int counter = 0, tmpx = 0, tmpy = 0;
    for(it = kplst.begin(); it != kplst.end(); it++)
    {
        crntPt = *it;
        if(crntPt->KP)
        {
            counter++;
        }
    }

    if(AFF_OUT)
    {
        fprintf(fp,"%d 8\n",counter);
    }
    else
    {
        fprintf(fp,"%d 8\n",counter);
    }

    for(it = kplst.begin(); it != kplst.end(); it++)
    {
        crntPt = *it;
        if(crntPt->KP)
        {
            tmpx = (int)round(crntPt->x*resize_rate0);
            tmpy = (int)round(crntPt->y*resize_rate0);
            if(AFF_OUT)
            {
                fprintf(fp,"%d %d %f %f %f %f %f %f\n"
                        ,tmpx,tmpy, crntPt->a, crntPt->b, crntPt->c, crntPt->iscale, crntPt->ori, crntPt->funcVal);
            }
            else
            {
                fprintf(fp,"%d %d 1 0 1 %f %f %f\n"
                        ,tmpx, tmpy, crntPt->iscale, crntPt->ori, crntPt->funcVal);
            }
        }
    }
    fclose(fp);
    return true;
}


bool AbstractDetector::saveKpVGG(vector<KeyPoint*> &kplst, const int kpnum, const char *dstfn,
                                 const float resize_rate0, bool AFF_OUT)
{
    ofstream outStrm(dstfn);
    if(!outStrm.is_open())
    {
        cout<<"File '"<<dstfn<<"' cannot open!\n";
        return false;
    }

    int counter = 0, tmpx = 0, tmpy = 0;
    vector<KeyPoint*>::iterator it;
    KeyPoint *crntPt = NULL;
    for(it = kplst.begin(); it != kplst.end(); it++)
    {
        crntPt = *it;
        if(crntPt->KP)
        {
            counter++;
        }
    }

    outStrm<<"1.0"<<endl;
    outStrm<<counter<<endl;
    float e;

    for(it = kplst.begin(); it != kplst.end(); it++)
    {
        crntPt = *it;
        if(crntPt->KP)
        {
            tmpx = (int)round(crntPt->x*resize_rate0);
            tmpy = (int)round(crntPt->y*resize_rate0);
            e = 1.0/crntPt->dscale;
            e = e*e;
            float sc = crntPt->iscale*crntPt->iscale;
            outStrm<<crntPt->fx<<" "<<crntPt->fy<<" "<<crntPt->a<<" "<<crntPt->b<<" "<<crntPt->c<<endl;
        }
    }
    outStrm.close();
    return true;
}

bool AbstractDetector::topkSelect(vector<KeyPoint*> &kplst, const int topk)
{
    KeyPoint *crntPt = NULL;
    int i = 0;
    vector<KeyPoint *>::iterator  itpt;

    stable_sort(kplst.begin(), kplst.end(), KeyPoint::keypCompF);
    for(itpt = kplst.begin(); itpt != kplst.end(); itpt++)
    {
        i++;
        crntPt = *itpt;

        if(i > topk)
        {
            crntPt->KP = false;
        }
        ///cout<<crntPt->KP<<endl;
    }
    return true;
}

bool AbstractDetector::topkEqDnSelect(vector<KeyPoint*> &kplst, const int w,const int h,const int topk)
{
    float ratio = (w + 0.0f)/(h + 0.0f);
    float area  = w*h/topk;
    float h1    = sqrt(area/ratio);
    float w1    = ratio*h1;
    int wbound  = (int)floor(w / w1);

    KeyPoint *crntPt;
    int key, xkey, ykey;

    map<int,vector<KeyPoint *>* >::iterator mit;
    map<int,vector<KeyPoint *>* > grid_map;
    vector<KeyPoint*>::iterator it;
    vector<KeyPoint *> *crnt_lst;

    for(it = kplst.begin(); it != kplst.end(); it++)
    {
        crntPt = *it;
        xkey   = (int)floor(crntPt->x/w1);
        xkey   = xkey > wbound?wbound:xkey;
        ykey   = (int)floor(crntPt->y/h1);
        ykey   = ykey > wbound?wbound:ykey;
        key    = wbound*ykey + xkey;
        if(grid_map.find(key) == grid_map.end())
        {
            crnt_lst = new vector<KeyPoint*>;
            crnt_lst->push_back(crntPt);
            grid_map.insert(pair<int,vector<KeyPoint*>* >(key,crnt_lst));
        }
        else
        {
            crnt_lst = grid_map[key];
            crnt_lst->push_back(crntPt);
            grid_map[key] = crnt_lst;
        }
    }

    ///sort first
    for(mit = grid_map.begin(); mit != grid_map.end(); mit++)
    {
        crnt_lst = mit->second;
        stable_sort(crnt_lst->begin(),crnt_lst->end(),KeyPoint::keypCompF);
    }

    int n_blocks = grid_map.size();
    int kp_num = kplst.size();
    int i = 0;
    int counter = 0;

    unsigned int *cusors = new unsigned int[n_blocks];
    memset(cusors, 0, sizeof(unsigned int)*n_blocks);
    bool DONE = false;
    ///choose keypoints from block to block
    do
    {
        i = 0;
        for(mit = grid_map.begin(); mit != grid_map.end() && !DONE; mit++)
        {
            crnt_lst = mit->second;
            if((cusors[i] +1) <= crnt_lst->size())
            {
                crntPt = (*crnt_lst)[cusors[i]];
                crntPt->KP = true;
                cusors[i]  = cusors[i] + 1;
                counter++;
                if(counter == topk)
                {
                    DONE = true;
                }
            }
            i++;
        }
    }
    while(!DONE && counter < kp_num);

    delete [] cusors;
    for(mit = grid_map.begin(); mit != grid_map.end(); mit++)
    {
        crnt_lst = mit->second;
        crnt_lst->clear();
        delete crnt_lst;
    }
    grid_map.clear();
    return true;
}

bool AbstractDetector::detectKeyPoints(const char *srcdir, const char *kpdir, const char *drdir,
                                       const char *descdir, const char *dvdir, bool ONE_MISSION)
{
    char kpfn[FNLEN],  drfn[FNLEN],  dvfn[FNLEN];
    char srcfn[FNLEN], fname[FNLEN], dstfn[FNLEN];
    struct dirent *result = NULL, entry;
    bool done = false;
    char detID[10];
    DIR *dp = NULL;

    strcpy(detID, detectorsID[DETECTOR]);

    int i = 0, fnlen, suf_len = 4;

    KPDrawer *mydrawer = new KPDrawer(this->CIRCLE_OUT, this->DETECTOR);

    if(this->DETECTOR == non)
    {
        if(!VString::existDir(kpdir))
        {
            cout<<"Keypoint files are required for 'non' detector!\n";
            exit(0);
        }
    }

    vector<KeyPoint*>::iterator it;
    clock_t start = 0, finish = 0;
    float tm_span = 0;

    cout<<"Keypoints detection ..................... for '"<<srcdir<<"'\n";
    if(descdir != NULL && strcmp(descdir,""))
    {
        IOTool::creatDIR(descdir);
    }

    if(kpdir != NULL && strcmp(kpdir, ""))
    {
        IOTool::creatDIR(kpdir);
    }

    if(drdir != NULL && strcmp(drdir, ""))
    {
        IOTool::creatDIR(drdir);
    }

    if(dvdir != NULL && strcmp(dvdir, ""))
    {
        IOTool::creatDIR(dvdir);
    }

    if(ONE_MISSION)
    {
        strcpy(srcfn, srcdir);

        VString::parseFName(fname, srcdir);

        if(!strcmp(fname,""))
        {
            cout<<"The input source image is not defined\n";
            exit(0);
        }

        if(kpdir != NULL)
        {
            strcpy(kpfn, kpdir);
            strcat(kpfn, fname);
            strcat(kpfn, KP_SFIX);
        }
        else
        {
            strcpy(kpfn, "");
        }

        if(descdir == NULL)
        {
            strcpy(dstfn, "");
        }
        else
        {
            strcpy(dstfn, descdir);
            strcat(dstfn, fname);
            strcat(dstfn, KEY_SFIX);
        }

        if(dvdir == NULL)
        {
            strcpy(dvfn, "");
        }
        else
        {
            strcpy(dvfn, dvdir);
            strcat(dvfn, fname);
            strcat(dvfn, detID);
            strcat(dvfn, DV_SFIX);
        }

        if(descdir != NULL && this->myDescriptor == NULL)
        {
            cout<<" Choose proper description with '-p' option!\n";
        }
        else
        {
            done = this->KeypointBuild(srcfn, kpfn, dstfn, dvfn);

            if(drdir == NULL)
            {
                strcpy(drfn, "");
            }
            else
            {
                strcpy(drfn, drdir);
                strcat(drfn, fname);
                strcat(drfn, detID);
                strcat(drfn, JPG_SFIX);
                if(this->DETECTOR == dsurf)
                {
                    mydrawer->draw_rects(this->kps, srcfn, this->resize_rate, drfn);
                }
                else
                {
                    mydrawer->draw_shapes(this->kps, srcfn, this->resize_rate, drfn);
                }
            }
            AbstractDetector::releaseKpList(this->kps);
        }
        return true;
    }


    if((dp  = opendir(srcdir)) == NULL)
    {
        cout<<" Source directory '"<<srcdir<<"' do not exist!\n";
        return false;
    }

    if(descdir != NULL && this->myDescriptor == NULL)
    {
        cout<<" Choose proper description with '-p' option!\n";
    }

    FILE *tm_fp = NULL;

    if(this->TIMEON)
    {
        tm_fp = fopen(TMLOG, "a");
    }

    while (readdir_r(dp, &entry, &result) == 0 && result != NULL)
    {
        if(entry.d_type != DT_REG)
        {
            continue;
        }

        if((!VString::endWith(entry.d_name, "jpg"))&&(!VString::endWith(entry.d_name, "bmp"))
                &&(!VString::endWith(entry.d_name, "pgm")) && (!VString::endWith(entry.d_name, "ppm"))
                &&(!VString::endWith(entry.d_name, "png")))
            continue;

        strcpy(srcfn, srcdir);

        fnlen = strlen(entry.d_name);
        strcat(srcfn, entry.d_name);

        if(!VString::validatePath(srcfn))
        {
            cout<<"Invalid file name: '"<<srcfn<<"'!!!\n";
            continue;
        }

        strncpy(fname, entry.d_name, fnlen-suf_len);
        fname[fnlen-suf_len] = '\0';

        if(kpdir == NULL)
        {
            strcpy(kpfn,"");
        }
        else
        {
            strcpy(kpfn, kpdir);
            strcat(kpfn, fname);
            strcat(kpfn, KP_SFIX);
        }

        if(descdir == NULL)
        {
            strcpy(dstfn, "");
        }
        else
        {
            strcpy(dstfn, descdir);
            strcat(dstfn, fname);
            strcat(dstfn, KEY_SFIX);
        }

        if(dvdir == NULL)
        {
            strcpy(dvfn, "");
        }
        else
        {
            strcpy(dvfn, dvdir);
            strcat(dvfn, fname);
            strcat(dvfn, detID);
            strcat(dvfn, DV_SFIX);
        }

        if(this->TIMEON)
        {
            start = clock();
        }

        done  = this->KeypointBuild(srcfn, kpfn, dstfn, dvfn);

        if(!done)
        {
            continue;
        }

        if(this->TIMEON)
        {
            finish = clock();
            tm_span = finish-start;
            tm_span = tm_span/CLOCKS_PER_SEC;
            fprintf(tm_fp,"%s\t%f\n", fname, tm_span);
        }

        if(drdir == NULL)
        {
            strcpy(drfn, "");
        }
        else
        {
            strcpy(drfn, drdir);
            strcat(drfn, fname);
            strcat(drfn, detID);
            strcat(drfn, JPG_SFIX);
            if(this->DETECTOR == dsurf)
            {
                mydrawer->draw_rects(this->kps, srcfn, this->resize_rate, drfn);
            }
            else
            {
                mydrawer->draw_shapes(this->kps, srcfn, this->resize_rate, drfn);
            }
        }

        AbstractDetector::releaseKpList(this->kps);
        i++;
        cout<<"\r\r\r\r\t"<<i;
        //exit(0);
    }
    cout<<endl;
    closedir(dp);
    if(this->TIMEON)
    {
        fclose(tm_fp);
    }
    delete mydrawer;
    return true;
}

bool AbstractDetector::mapPCA(const float *feats, const float *pcamat, vector<float> &means,
                            vector<float> &var, float *pcaft, const int dim, const int pdim)
{
    int i = 0, k = 0, loc = 0;
    assert(pcaft != NULL);
    assert(feats != NULL);
    assert(dim >= pdim);
    memset(pcaft, 0, sizeof(float)*pdim);
    for(i = 0; i < pdim; i++)
    {
        loc = i*dim;
        for(k = 0; k < dim; k++)
        {

            pcaft[i] += pcamat[loc+k]*(feats[k] - means[k]);
        }
        pcaft[i] = pcaft[i]/var[i];
    }
    return true;
}

bool AbstractDetector::proj(const float *feats, const float *prjMat,
                            float *prjFt, const int dim, const int pdim)
{
    int i = 0, k = 0, loc = 0;
    assert(prjFt != NULL);
    assert(feats != NULL);
    assert(dim >= pdim);

    memset(prjFt, 0, sizeof(float)*pdim);
    for(i = 0; i < pdim; i++)
    {
          loc = i*dim;
        for(k = 0; k < dim; k++)
        {
            prjFt[i] += prjMat[loc+k]*feats[k];
        }
    }

    return true;
}

/**SURF's approach***/
void AbstractDetector::FindOrientation(vector<KeyPoint*> &keyps)
{
    assert(this->intImg != NULL);
    float weight = 0.f;
    int s, r, c, i, j;
    float resX[109], resY[109], Ang[109];
    const int id[13] = {6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6};
    unsigned int k, idx = 0;
    float sumX = 0.0f, sumY = 0.f, mag = 0.0f;
    float max  = 0.0f, orientation = 0.f;
    float ang1 = 0.0f,  ang2 = 0.f, ang, deltA = PI/6.0;
    vector<KeyPoint*>::iterator vit;
    KeyPoint *crnt_kpt;

    for(vit = keyps.begin(); vit != keyps.end(); vit++)
    {
        crnt_kpt = *vit;
        memset(resX, 0, sizeof(float)*109);
        memset(resY, 0, sizeof(float)*109);
        memset(Ang,  0, sizeof(float)*109);
        idx = 0;
        s = round(crnt_kpt->dscale);
        r = round(crnt_kpt->y);
        c = round(crnt_kpt->x);
        gauss25 = Filter::GaussianKernel2D(2.5*s, 6);
        for(i = -6; i <= 6; i++)
        {
            for(j = -6; j <= 6; j++)
            {
                if(i*i + j*j < 36)
                {
                    weight    = gauss25[id[i+6]][id[j+6]];
                    /**this version works well also**********
                    resX[idx] = weight * this->getDx(c+j, r+i, 2*s, 4*s);
                    resY[idx] = weight * this->getDy(c+j, r+i, 2*s, 4*s);
                    ***/
                    resX[idx] = weight * this->getDx(c+j*s, r+i*s, 2*s, 4*s);
                    resY[idx] = weight * this->getDy(c+j*s, r+i*s, 2*s, 4*s);
                    Ang[idx]  = getAngle(resX[idx], resY[idx]);
                    idx++;
                }
            }
        }
        Cleaner::clear2DArray(gauss25);

        sumX = sumY = mag = 0.0f;
        max  = orientation = 0.0f;
        ang1 = ang2 = ang = 0;

        for(ang1 = 0; ang1 < PI2;  ang1 += 0.15f)
        {
            ang2 = ang1 + deltA;
            while(ang2 > PI2)
                ang2 = ang2 - PI2;

            sumX = sumY = 0.f;
            for(k = 0; k < 109; ++k)
            {
                ang = Ang[k];
                if(ang1 < ang2 && ang1 < ang && ang < ang2)
                {
                    sumX += resX[k];
                    sumY += resY[k];
                }
                else if (ang2 < ang1 && ((ang > 0 && ang < ang2)
                                         || (ang > ang1 && ang < PI2)))
                {
                    sumX += resX[k];
                    sumY += resY[k];
                }
            }
            mag = sumX*sumX + sumY*sumY;
            if (mag > max)
            {
                max = mag;
                orientation = getAngle(sumX, sumY);
            }
        }

        while(orientation < 0.0)
            orientation += PI2;

        while(orientation >= PI2)
            orientation -= PI2;

        crnt_kpt->ori = orientation;
    }//for(vit)
    return ;
}

int AbstractDetector::findFlip(const float thetas[], const int dbin, const int nbin, float &rate)
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

    rate = righta > lefta? (lefta/righta):(righta/lefta);
    if(righta > lefta)
    {
        return  1;
    }
    else
    {
        return -1;
    }
}

IntImage *AbstractDetector::buildIntImage(Image *img0)
{
    unsigned int h = img0->height, w = img0->width;
    unsigned int i = 0, j, loc = 0;
    IntImage *nwintImg = new IntImage(w, h);
    double val = 0;
    int x, y;

    for(i = 0; i < h; i++)
    {
        y = i;
        for(j = 0; j < w; j++)
        {
            x = j;
            loc = i*w + j;     ///[x, y]
            val = nwintImg->getPixel(x, y-1) + nwintImg->getPixel(x-1, y)
                  + img0->pix[loc] - nwintImg->getPixel(x-1, y-1);
            nwintImg->setPixel(x, y, val);
        }
    }

    return nwintImg;
}


float AbstractDetector::getAngle(float dx, float dy)
{
    if(dx >= 0 && dy >= 0)
        return atan(dy/dx);

    if(dx < 0 && dy >= 0)
        return PI - atan(-dy/dx);

    if(dx < 0 && dy < 0)
        return PI + atan(dy/dx);

    if(dx >= 0 && dy < 0)
        return 2*PI - atan(-dy/dx);

    return 0;
}

float AbstractDetector::getDx(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral(x0, (y0 - rd0), rd0, ts) - this->intImg->boxIntegral(x0-rd0, y0-rd0, rd0, ts);
    return val;
}

float AbstractDetector::getDy(const int x0, const int y0, const int rd0, const int ts)
{
    float val = this->intImg->boxIntegral(x0-rd0, y0, ts, rd0) - this->intImg->boxIntegral(x0-rd0, y0-rd0, ts, rd0);
    return val;
}


bool AbstractDetector::releaseKpList(vector<KeyPoint *> &keyps)
{
    vector<KeyPoint *>::iterator it;
    KeyPoint *crntPt;
    for(it = keyps.begin(); it != keyps.end(); it++)
    {
        crntPt = *it;
        delete crntPt;
    }
    keyps.clear();
    return true;
}


bool AbstractDetector::releaseScales(vector<Image*> &scales)
{
    vector<Image *>::iterator it;
    Image *tempim;
    for(it=scales.begin(); it!=scales.end(); it++)
    {
        tempim = *it;
        delete tempim;
    }
    scales.clear();
    return true;
}


AbstractDetector::~AbstractDetector()
{
    vector<vector<float> >::iterator it;
    for(it = gauss25.begin(); it != gauss25.end(); it++)
    {
        vector<float> &crnt_vect = *it;
        crnt_vect.clear();
    }
    gauss25.clear();

    if(this->crntimg != NULL)
    {
        delete this->crntimg;
        this->crntimg = NULL;
    }

    if(this->myDescriptor != NULL)
    {
       delete this->myDescriptor;
       this->myDescriptor = NULL;
    }


}

