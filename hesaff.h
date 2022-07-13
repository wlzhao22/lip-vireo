#ifndef HESAFF_H
#define HESAFF_H


#include <algorithm>
#include <fstream>
#include <cmath>
#include "abstractdetector.h"
#include "intimage.h"

class HesAff : public AbstractDetector {
public:
    HesAff();

    virtual ~HesAff();

private:

    int numb;
    float thresh;
    //float sigma;
    static const int MaxOctaves;
    static const int SCALES;
    static const float _SIGMA;
    static const bool INTERP_KEYS;
    static const float INITSIGMA;
    static bool ZOOM_OUT;

    static bool LOGON;
    static bool useD1L;             // 是否使用 D1L>0作为辅助阈值
    static bool useThresh;          // 函数值是否要大于阈值
    static const float D1L_ALPHA;         // <Image Matching Using Generalized Scale-Space Interest Points> Tony Lindeberg

    //for detector
private:

    static const int BORDER;        // keypoints必须在图片的BORDER以内
    static const int THRESH;        // Hessian行列式的阈值
    static const float mag;         // integration scale = mag*derivative scale

    // FindOrientByGrad 中 生成直方图的参数
    static const int DEGREE; // 方向直方图中每一柱的角度范围, 共有360/DEGREE柱
    static const int NOrient;
    static const int DEGPERBIN;
    static const float NwKpThresh;

    // adapt affine 用到的参数
    static const int maxAdaptIter; // 允许的最大迭代次数
    static const float EPSILON;    // 浮点数最小精度
    vector<Image *> octave;
    vector<vector<Image*> > pyramid;
    vector<vector<Image*> > dxxPyramid, dyyPyramid, dxyPyramid;

    Image *CreateInitialImage(Image *crntImg);
    void BuildGaussianPyramid(Image *initImg, vector<vector<Image*> > &pyramid, const float sigma0, const int numScales);
    void BuildHessianPyramid(vector<vector<Image*> > &gaussPyramid, vector<vector<Image*> > &hessPyramid, vector<vector<Image*> > &D1LPPyramid);
    void BuildHessianImage(Image *src, float sigma, Image *&hessImg, Image *&D1LImg);
    ///void BuildHessianImage(Image *src, float sigma, Image *&hessImg, Image *&D1LImg, Image *&DxxImg, Image *&DxyImg, Image *&DyyImg);
    void FindKeypoints(const int octIndex, vector<Image *> &hessOctave, vector<Image *> &D1Loctave, vector<KeyPoint *> &peaks);
    KeyPoint *FindKeypointNearby(Image *hessImg, Image *D1LImg, const int cx, const int cy, const int size);
    bool  InterpKey(int x, int y, int s, vector<Image *> &LoGImages, float *fx, float *fy, float *fs, float *dogVal);
    float InterpKeyStep(int x, int y, int s, vector<Image *> &DI, float *dx, float *dy, float *ds);

    vector<KeyPoint*> FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GaussianOctaves);
    void affineAdapt(vector<KeyPoint *> &peaks,   vector<vector<Image*> > &gaussPyramid, vector<vector<Image*> > &dxPyramid,  vector<vector<Image*> > &dyPyramid);
    void affineAdapt2D(vector<KeyPoint *> &peaks, vector<vector<Image*> > &gaussPyramid);

    Image *normalizeWindow(Image *window, float sqrtU[2][2], int cx, int cy);
    void buildDxPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &DxPyramid);
    void buildDyPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &DyPyramid);
    void buildDxImage(Image *src, Image *dxImg);
    void buildDyImage(Image *src, Image *dyImg);

    void cal2ndMomentMat(Image *window, const int x, const int y, const float dscale, const float iscale, float u[2][2]);

    void saveKeypoints(const string dstFn, vector<KeyPoint *> kps)
    {
        vector<KeyPoint*>::iterator it;
        KeyPoint* kp = NULL;
        stable_sort(kps.begin(), kps.end(), KeyPoint::keypCompF);
        FILE *fp = fopen(dstFn.c_str(), "w");
        for(it = kps.begin(); it != kps.end(); it++)
        {
            kp = *it;
            std::fprintf(fp, "%d %d %.4f %.4f %.4f %.4f %.4f\n",
                         kp->x, kp->y, kp->dscale, kp->a, kp->b, kp->c, kp->funcVal);
        }
        fclose(fp);
    }

public:

    // (cx,cy)是否是(x,y,s)空间上的局部极值
    static bool isLocalExtrema(const Image *prevImg, const Image *img, const Image *nextImg, const int cx, const int cy);
    // (cx,cy)是否是(x,y)空间上的局部极值
    static bool isLocalExtrema(const Image *img, const int cx, const int cy);

    // 2*2矩阵的行列式
    static inline float det2D(float m[2][2]) {
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    }

    // Scale-Space Behaviour of Local Extrema and Blobs,equation(7)
    // 根据隐函数定理求trajectory关于sigma的导数
    // 但是论文中有一项错了,-LxyLxxx应该改一个变成-LxyLxyy
    static bool getdRds(const float block[5][5], const float sigma, float &dxds, float &dyds);

    static float getdHds(const float block[5][5], float sigma);
    // 2*2矩阵的特征值
    static inline bool Eigenvalue2D(float m[2][2], float &eig1, float &eig2) {
        // 解方程（a-eig1)(d-eig2)-bc=0
        float a = 1, b = -(m[0][0] + m[1][1]), c = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        float delta = b * b - 4 * a * c;

        if (delta < -EPSILON)
            return false;
        else {
            if (delta < 0) delta = -delta;
            float sqrtDelta = sqrt(delta);
            eig1 = (-b - sqrtDelta) / (2 * a);
            eig2 = (-b + sqrtDelta) / (2 * a);
            return true;
        }
    }

    // 2*2矩阵的逆
    static inline bool Inverse2D(float m[2][2], float invm[2][2]) {
        float det = det2D(m);
        if (fabs(det) > EPSILON) {
            float divDet = 1 / det;
            invm[0][0] = divDet * m[1][1];
            invm[0][1] = -divDet * m[0][1];
            invm[1][0] = -divDet * m[1][0];
            invm[1][1] = divDet * m[0][0];
            return true;
        } else
            return false;
    }

    static inline bool Inverse2D(float m[2][2]) {
        float det = det2D(m);
        if (fabs(det) > EPSILON) {
            float divDet = 1 / det;
            float a = m[0][0];
            m[0][0] = divDet * m[1][1];
            m[0][1] = -divDet * m[0][1];
            m[1][0] = -divDet * m[1][0];
            m[1][1] = divDet * a;
            return true;
        } else
            return false;
    }

    // 2*2矩阵的平方根, 只计算出实根
    // https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
    static inline bool SqrtMatrix2D(float m[2][2], float mr[2][2]) {
        float det = det2D(m);
        if (det > EPSILON) {
            det = sqrt(det);
            float t = m[0][0] + m[1][1] + 2 * det;
            if (t > EPSILON) {
                t = 1 / sqrt(t);
                mr[0][0] = t * (m[0][0] + det);
                mr[0][1] = t * m[0][1];
                mr[1][0] = t * m[1][0];
                mr[1][1] = t * (m[1][1] + det);
                return true;
            }
        }
        return false;
    }

    // 将initImg的[x1~x2, y1~y2]部分切割出来并返回
    static inline Image *ImageCut(Image *initImg, int x1, int x2, int y1, int y2) {
        int width = x2 - x1 + 1, height = y2 - y1 + 1;
        Image *newImg = new Image(width, height);
        for (int i = 0; i < width; i++)
            for (int j = 0; j < height; j++)
                newImg->setPixel(i, j, initImg->getPixel(x1 + i, y1 + j));
        return newImg;
    }

    // 返回x的符号 正数1 负数-1 0
    static inline int Sign(const float x) { return (fabs(x) < EPSILON ? 0 : (x < -EPSILON ? -1 : 1)); }

    // 计算window(cx,cy,sigma)的det(H)对sigma的导数
    static float getdHds(const Image *window, const int cx, const int cy, const float sigma);

    // 计算normalized det(Hessian(cx, cy, sigma))
    static float getDetH(const Image *window, const int cx, const int cy, const float sigma);

    // m1矩阵右乘m2
    inline static void MultiplyMatrix(float m1[2][2], float m2[2][2]) {
        float a = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0];
        float b = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1];
        float c = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0];
        float d = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1];
        m1[0][0] = a;
        m1[0][1] = b;
        m1[1][0] = c;
        m1[1][1] = d;
    }

    // m/sqrt(det(m))
    inline static void Normalize(float m[2][2]) {
        float v = 1 / sqrt(det2D(m));
        m[0][0] *= v;
        m[0][1] *= v;
        m[1][0] *= v;
        m[1][1] *= v;
    }


public:
    bool paramsCheck();
    bool KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn);
    void writeKeypoint(const char *fn);
    void extractKeyp(string imgPath);
    static void test();

};

#endif // HESAFF_H
