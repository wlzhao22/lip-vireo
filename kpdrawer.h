#ifndef KPDRAWER_H
#define KPDRAWER_H

#include "idetector.h"
#include "keypoint.h"
#include "cimage.h"

#include <vector>

class KPDrawer
{
    private:
        static const unsigned int _cwith, _cheight, _awidth, _aheight;
        static const float PI0;
        float WHITE[3];
        float RED[3];
        float PURPL[3];
        float YELLW[3];
        bool draw_circle;
        Detector myoption;

    protected:

        void draw_ellipse(CImage *img, const unsigned int x0, const unsigned int y0, const float a,
                           const float b, const float c, const float sc);
        void draw_cross(CImage *img, const unsigned int x0, const unsigned int y0);
        void draw_arrow(CImage *img, const float theta0, const unsigned int x0,
                        const unsigned int y0);

        void lineto(CImage *Img, const int x0, const int y0, const int x1, const int y1);
        void lineto(CImage *Img, const int x0, const int y0,
                    const int x1, const int y1, const float *color);

    public:
        KPDrawer(const bool CIRCLE, const Detector det_option);

        void draw_rects(vector <KeyPoint *> kps, const char *srcimgfn, const float scale_rate0, const char *dstimgfn);
        void draw_shapes(vector <KeyPoint *> kps, const char *srcimgfn, const float scale_rate0, const char *dstimgfn);

        static void load_keys(vector<KeyPoint*> &kps, const char *srcfn);
        static void load_siftgeo(vector<KeyPoint*> &kps, const char *srcfn);
        static void clear_kps(vector<KeyPoint*> &kps);

        static void test();
        virtual ~KPDrawer(){}
};

#endif // KPDRAWER_H
