#ifndef DENSE_H
#define DENSE_H

#include "abstractdetector.h"


/****
step=5 (only rotate, do not square-root)
dense+nrsift  Holidays     VLAD*       0.697
dense+nrsift  Ox5k         VLAD*       0.313

step=5 (rotate, square-root)
dense+nrrsift  Holidays     VLAD*       0.772
dense+nrrsift  Ox5k         VLAD*       0.381

***/

class Dense: public AbstractDetector
{
    public:
        Dense();
        virtual ~Dense(){};

        int numb, step;
        int BORDER;
        float sigma0;

     private:
        static const int MaxOctaves;
        static const int SCALES;
        static const float _SIGMA;
        static const float INITSIGMA;

        static const float k0;
        static const int THRESH;
        static const float mag;
        static const int DEGREE;
        static const int Delta;

        static const int NumOrient;
        static const int DEGPERBIN;
        static const float NwKpThresh;

        void   BuildOctaves(Image * image, GBlur blur_opt, vector<vector<Image *> > &Ggoctaves, const float _SIGMA0);
        vector<Image*> BuildGaussianScales(Image * image, GBlur blur_opt, const int noctave,const float _SIGMA0);
        vector<KeyPoint *> FindOrientation(vector<KeyPoint *> & kps, vector<vector<Image *> > & GOctaves);

     public:

        int keypDetect(const float sigma,const int level,const int scale);
        vector<KeyPoint *> generatePeaksScales(const int octave, vector<Image *> & GScales);
        void generatePeaksOctaves(vector<vector<Image *> > & GOctaves);

        bool paramsCheck();
        bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char* dvfn);
        void writeKeypoint(const char*fn);

    static void test();
};

#endif
