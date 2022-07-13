#include <iostream>
#include <cstdlib>

#include "kpdrawer.h"
#include "vstring.h"
#include "pallete.h"
#include "ioimage.h"
#include "cimage.h"
#include "filter.h"
#include "image.h"
#include "vmath.h"
#include "haar.h"

#include "abstractdetector.h"
#include "nondetector.h"
#include "hessian.h"
#include "hesslap.h"
#include "hessnaff.h"
#include "harris.h"
#include "harlap.h"
#include "hesaff.h"
#include "dsurf.h"
///#include "msurf.h"
#include "dense.h"
#include "dog.h"
#include "log.h"

#include "descsift.h"
#include "desccm.h"
#include "thinner.h"
#include "canny.h"

#include "kernel.h"

#include <cstdlib>

#include <cassert>
#include <string>
#include <ctime>
#include <map>

using namespace std;

/***
*@author  Wan-Lei Zhao
*@version 1.38
*@date    18-May-2018
*
*
*All rights reserved by Wan-Lei Zhao
*
*
*
**/

void usage()
{
    const char *version = "1.40 stable";
    cout<<" Usage:\n\n";
    cout<<" lip-vireo [-img fn |-dir dir] [-d detector] [-kpdir kpdir] [-drdir imgdir] [-p descriptor] [-dsdir descdir] [-dvdir pviewdir] [-c conf]\n\n";
    cout<<" -img fn |-dir dir\n";
    cout<<"\tpgm/ppm/png/bmp/jpg image or folder contains pgm/ppm/png/bmp/jpg image(s)\n\n";
    cout<<" -d hess|dog|harr|log|harlap|hesslap|hesaff|dense|dsurf|non\n";
    cout<<"\thess\tBased on Hessian matrix\n";
    cout<<"\thesaff\tHessian-Affine detector\n";
    cout<<"\thesslap Based on Hessian-Laplacian function\n";
    cout<<"\tdog\tBased on Difference of Gaussian pyramid\n";
    cout<<"\tlog\tBased on Laplacian of Gaussian function\n";
    cout<<"\tharr\tBased on Harris function (multiple scales)\n";
    cout<<"\tharlap\tBased on Harris-Laplacian function\n";
    cout<<"\tdsurf\tFast Hessian detector\n";

    cout<<"\tdense\tDense sampling on multiple octaves and multiple scales\n";
    cout<<"\tnon\tKeypoints are supplied by the user\n\n";
    cout<<" -p SIFT|LJET|ERIFT|FIFT|SPIN|RIFT|SURF|ESURF|AOD|FIND\n";
    cout<<"\tSIFT\tScale Invariant Feature Transform\n";
    cout<<"\tRIFT\tRotation Invariant Feature Transform\n";
    cout<<"\tERIFT\tEnhanced Rotation Invariant Feature Transform\n";
    cout<<"\tFIFT\tFlip invariant SIFT\n";
    cout<<"\tSPIN\tSpin images \n";
    cout<<"\tLJet\tSteerable filter (also known as 'local jet')\n";
    cout<<"\tSURF\tSpeed-Up Robust Features (based on Haar wavelet)\n";
    cout<<"\tESURF\tExtended SURF (128 dimensions)\n";
    cout<<"\tAoD\tAggregation on Derivatives (Dx, Dy, |Dx|, |Dy|)\n";
    cout<<"\tFIND\tFlip INvariant Descriptor proposed by X.-J. Guo\n\n";
    cout<<" -kpdir kpdir\n";
    cout<<"\tDirectory where to save keypoint file(s), optional except for '-d non'\n\n";
    cout<<" -drdir imgdir\n";
    cout<<"\tDirectory where to save image(s) with keypoints plotted (optional)\n\n";
    cout<<" -dsdir descdir\n";
    cout<<"\tDirectory where to save descriptor file(s) (optional)\n\n";
    cout<<" -dvdir pviewdir\n";
    cout<<"\tDirectory where to save image(s) of keypoint patches (optional)\n\n";
    cout<<" -c conf\n";
    cout<<"\tConfigure file for parameters\n\n";

    cout<<"\tAuthor:   Wan-Lei Zhao, wlzhao@xmu.edu.cn\n";
    cout<<"\tVersion:  "<<version<<" (2008-2018)\n";
}

void license()
{
    cout<<"\t\t\t  Copyright\n\n";
    cout<<" 1. All rights are reserved by the author and the patent holder.\n";
    cout<<" Permission to use, copy, and distribute this software and its do-\n";
    cout<<" cumentation is hereby granted free of charge, provided that (1) it\n";
    cout<<" is not a component of a commercial product, and (2) this notice\n";
    cout<<" appears in all copies of the software and related documentation;\n";
    cout<<" 2. DoG+SIFT has been patented by David G. Lowe in US. The patent\n";
    cout<<" holder is University of British Columbia, Canada;\n";
    cout<<" 3. Copyright holder for DSURF+SURF is Tinne Tuytelaars et al.;\n";
    cout<<" 4. Copyright holder for FIND is Xiao-Jie Guo et al.\n";
    cout<<"\n\n";
}

void test1()
{
    return ;
}

void test()
{
    vector<vector<float> > kern;
    ///Harris::test();
    ///HarLap::test();
    ///HessLap::test();
    ///HesAff::test();
    ///Hessian::test();
    ///DoG::test();
    ///LoG::test();
    ///Dense::test();
    ///DSURF::test();

    //IOImage::test();
    //IOImage::testbmp();
    //Kernel::test();
    //Haar::test();
    ///Filter::GaussianKernel2D(1.4, kern);
    ///Filter::printMat(kern);
    //CImage::test();
    //Image::test();
    ///VMath::test();
    //Pallete::test();

    //DescCM::test();
    //DescSIFT::test();

    //Canny::test();
    //CannyCorner::test();
    //Thinner::test();

    //KPDrawer::test();
    //VString::test();
    //license();

    ///test1();

    ///cout<<sizeof(int)<<endl;
    //cout<<sizeof(unsigned long long)<<endl;
    //usage();
    ///cout<<floor(-0.1)<<endl;
}

int main(int argc,char *argv[])
{
    /**
    test();
    return 0;
   /**/

    const char *args[9] = {"-img", "-dir", "-d", "-kpdir", "-p","-dsdir",
                           "-c",   "-drdir", "-dvdir"};

    const char *detectors[13] = {"hess", "dog", "harr", "non", "log", "harlap", "hesslap",
                                 "mser", "dense", "dsurf", "msurf", "corner", "hesaff"};

    int  i = 0, choice = 0, charslen = 0;
    char src_path[2048] = "";
    char dt_option[64]  = "";
    char *kp_path       = NULL;
    char *dr_path       = NULL;
    char *dv_path       = NULL;
    char conf[2048]     = "";
    char *dst_path      = NULL;
    char descOpt[64]    = "";

    map<string, int> argsmap;
    map<string, int> detectormap;

    if(argc < 9)
    {
        license();
        usage();
        return 1;
    }

    for(i = 0; i < 9; i++)
    {
        argsmap.insert(pair<string, int>(args[i], i));
    }

    for(i = 0; i < 13; i++)
    {
        detectormap.insert(pair<string, int>(detectors[i], i));
    }

    bool SUCCESS    = true;
    bool ONE_MISSON = true;
    bool SET_SRC    = false;

    for(i = 1; i < argc; i = i+2)
    {
        if(argsmap.find(argv[i]) != argsmap.end())
        {
            choice = argsmap[argv[i]];
            switch(choice)
            {
            case 0:
            {
                strcpy(src_path,argv[i+1]);
                if(SET_SRC)
                {
                    cout<<"argument '-dir', '-img' is allowed to select either one at once!\n";
                    SUCCESS = false;
                }
                else
                {
                    SET_SRC = true;
                }
                ONE_MISSON = true;

                break;
            }
            case 1:
            {
                strcpy(src_path, argv[i+1]);
                if(SET_SRC)
                {
                    cout<<"argument '-dir', '-img' is allowed to select either one at once!\n";
                    SUCCESS = false;
                }
                else
                {
                    SET_SRC = true;
                }
                ONE_MISSON = false;
                break;
            }
            case 2:
            {
                strcpy(dt_option, argv[i+1]);
                VString::toLower(dt_option);

                if(detectormap.find(dt_option) == detectormap.end())
                {
                    cout<<"Unknow detector option '"<<argv[i+1]<<"'\n";
                }
                break;
            }
            case 3:
            {
                charslen = strlen(argv[i+1]) + 1;
                kp_path = new char[charslen];
                strcpy(kp_path, argv[i+1]);
                break;
            }
            case 4:
            {
                strcpy(descOpt, argv[i+1]);
                break;
            }
            case 5:
            {
                charslen = strlen(argv[i+1]) + 1;
                dst_path = new char[charslen];
                strcpy(dst_path, argv[i+1]);
                break;
            }
            case 6:
            {
                strcpy(conf,argv[i+1]);
                break;
            }
            case 7:
            {
                charslen = strlen(argv[i+1]) + 1;
                dr_path = new char[charslen];
                strcpy(dr_path, argv[i+1]);
                break;
            }
            case 8:
            {
                charslen = strlen(argv[i+1]) + 1;
                dv_path = new char[charslen];
                strcpy(dv_path, argv[i+1]);
                break;
            }
         };
      }
      else
      {
            cout<<"Option '"<<argv[i]<<"' is not valid!\n";
            cout<<"valid options are '-img|-dir', '-d', '-p', '-kpdir', '-drdir', '-dsdir, '-c'\n";
      } ///end if-else
   }///for-loop

    if(!SUCCESS)
        return 0;

    IDetector *mydetector;
    VString::toLower(dt_option);
    VString::toUpper(descOpt);

    if(!strcmp(dt_option, detectors[0]))
    {
        mydetector = new Hessian();
    }
    else if(!strcmp(dt_option, detectors[1]))
    {
        mydetector = new DoG();
    }
    else  if(!strcmp(dt_option, detectors[2]))
    {
        mydetector = new Harris();
    }else if(!strcmp(dt_option, detectors[3]))
    {
        if(!strcmp(kp_path,""))
        {
            cout<<"Folder that contains keypoint files must be supplied for 'non' detector!\n";
            exit(1);
        }
        mydetector = new NonDetector();
    }
    else if(!strcmp(dt_option, detectors[4]))
    {
        mydetector = new LoG();
    }
    else if(!strcmp(dt_option, detectors[5]))
    {
        mydetector = new HarLap();
    }
    else if(!strcmp(dt_option, detectors[6]))
    {
        mydetector = new HessLap();
    }
    else if(!strcmp(dt_option, detectors[8]))
    {
        mydetector = new Dense();
    }
    else if(!strcmp(dt_option, detectors[9]))
    {
        mydetector = new DSURF();
    }
    else
    {
        exit(1);
    }

    if(!strcmp(conf, ""))
    {
        cout<<"Configure file should be set!\n";
        return 0;
    }

    mydetector->Init(conf, descOpt);

    mydetector->detectKeyPoints(src_path, kp_path, dr_path, dst_path, dv_path, ONE_MISSON);

    if(kp_path != NULL)
        delete [] kp_path;

    if(dr_path != NULL)
        delete [] dr_path;

    if(dst_path != NULL)
        delete [] dst_path;

    return 0;
}
