#ifndef IDETECTOR_H
#define IDETECTOR_H

#ifndef FNLEN
#define FNLEN 1024
#endif

#ifndef DIRLEN
#define DIRLEN 1024
#endif

#define KEY_SFIX ".pkeys"
#define KP_SFIX  ".keys"
#define JPG_SFIX "_dr.jpg"
#define DV_SFIX  "_dv.jpg"

#define TMLOG    "./time_cost.txt"

/**
Interface for keypoint detectors
**/
enum SelOPT {TOPK = 0, DENS = 1, THRSH = 2, NONS = 3};
enum DESC_FMT {_vgg_fmrt = 0, _vireo_fmrt = 1, _train_fmrt = 2};
enum Detector {hessian = 0, dog, harris, non, _log, harlap, hesslap,
               mser, dense, dsurf, msurf, corner, hessaff, hesaff};

enum DESC{
  NSIFT = 0,  SIFT = 1,    NLJET = 2,  LJET = 3,   NCM = 4,   CM = 5,     YKPCA = 6, NSPIN = 10,
  SPIN = 11,  NERIFT = 12, ERIFT = 13, NFIFT = 14, FIFT = 15, NGLOH = 16, GLOH = 17, NPCASIFT = 18,
  PCASIFT = 19, NRIFT = 20, RIFT = 21, NSURF = 22, SURF = 23,  NESURF = 24, ESURF = 25,
  NAOD = 26, AOD = 27, PCAPATCH = 28,  NASIFT = 30, ASIFT = 31, NGBLUR = 32, GBLUR = 33,
  NFIND = 34, FIND = 35, NPVIEW = 36, PVIEW = 37,  NCVSIFT = 38, CVSIFT = 39, CTSIFT = 40,
  NCTSIFT = 41, NCNN = 42, CNN = 43};

enum KP_FEAT {_scale_ = 0, _angle_ = 1, _flip_ = 2};


class IDetector
{
    protected:
        float resize_rate;
    public:
        static const int Numb_PROP = 5;
        bool kp_property[Numb_PROP];
        Detector DETECTOR;
        DESC     mydesc;
        bool AFF_OUT, CIRCLE_OUT;
        DESC_FMT out_FMRT;
        bool TIMEON;
    public:
        virtual bool Init(const char *scriptfn, const char *descopt) = 0;
        virtual bool paramsCheck() = 0;
        virtual bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn) = 0;
        virtual bool detectKeyPoints(const char *srcdir, const char *kpdir, const char *drdir,
                                     const char *destdir,const char *dvdir, bool ONE_MISSION) = 0;
        IDetector(){}
        virtual ~IDetector(){}
};

#endif
