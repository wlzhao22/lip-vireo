#include "hessaff.h"

/**constant initialization**/
#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <limits>
#include <fstream>

/**constant initialization**/

const int HessAff::MaxOctaves = 8;
const int HessAff::SCALES = 6;
const float HessAff::_SIGMA = 1.4;
const float HessAff::INITSIGMA = 0.5;
const bool  HessAff::INTERP_KEYS = false; // TODO true?

bool HessAff::ZOOM_OUT = false;  // 是否被放大两倍

const float HessAff::mag = 3.5; //optimal 8
const int HessAff::BORDER = 5;
const int HessAff::THRESH = 300;
const float HessAff::D1L_ALPHA = 0.05; // TODO D1L_ALPHA

const int HessAff::DEGREE = 10;
const int HessAff::NOrient = 36;
const int HessAff::DEGPERBIN = (360 / NOrient);
const float HessAff::NwKpThresh = 0.8;

const int HessAff::maxAdaptIter = 30;
const float HessAff::EPSILON = 1e-6;
bool HessAff::useD1L = true;
bool HessAff::useThresh = true;

#define _dx(src, x, y)  ((src->getPixel(x + 1, y) - src->getPixel(x - 1, y)) / 2);
#define _dy(src, x, y)  ((src->getPixel(x, y + 1) - src->getPixel(x, y - 1)) / 2);
#define _dxx(src, x, y) (src->getPixel(x + 1, y) + src->getPixel(x - 1, y) - 2 * src->getPixel(x, y));
#define _dyy(src, x, y) (src->getPixel(x, y + 1) + src->getPixel(x, y - 1) - 2 * src->getPixel(x, y));
#define _dxy(src, x, y) ((src->getPixel(x + 1, y + 1) + src->getPixel(x - 1, y - 1)- src->getPixel(x - 1, y + 1) - src->getPixel(x + 1, y - 1)) / 4);


//TODO 输出迭代过程
bool HessAff::LOGON = false;

HessAff::HessAff() {
    cout << "Detector ................................ HessAff\n";
    this->DETECTOR = hessaff;
    this->intImg = NULL;
    this->sel_option = THRSH;
}

bool HessAff::paramsCheck() {
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    if (this->paras.find(argv[1]) != this->paras.end()) {
        this->thresh = (float) atof(paras["thresh"]);
    } else {
        this->thresh = THRESH;
    }
    cout << "Thresh .................................. " << this->thresh << endl;

    if (this->paras.find(argv[2]) != this->paras.end()) {
        this->fix_kp_numb = atoi(paras["topk"]);
        this->sel_option = TOPK;
        cout << "Topk .................................... " << this->fix_kp_numb << endl;
    }

    if (this->paras.find(argv[3]) != this->paras.end()) {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option = DENS;
        cout << "Topk .................................... " << this->fix_kp_numb << endl;
    }

    return true;
}

Image *HessAff::CreateInitialImage(Image *crntImg) {
    Image *tmpImg;
    int minsz = min(crntImg->width, crntImg->height);
    vector<Image *> initImgs;

    if (minsz <= 200) {
        tmpImg = Image::doubleSizeImage(crntImg);
        this->resize_rate = 0.5;
        this->RESIZE = true;
        ZOOM_OUT = true;
    } else {
        this->resize_rate = 1.0f;
        ZOOM_OUT = false;
        // 如果size超过限制, 则每次缩小一半, 直到小于最大允许的尺寸
        tmpImg = crntImg->clone();
        while (minsz > imgSzBound) {
            tmpImg = Image::halfSizeImage(tmpImg);
            minsz = min(tmpImg->width, tmpImg->height);
            this->resize_rate = this->resize_rate * 2;
            this->RESIZE = true;
        }
    }

    // 将初始图片的尺度设置为_SIGMA
    float dsigma;
    if (ZOOM_OUT)
        dsigma = sqrt(_SIGMA * _SIGMA - 4 * INITSIGMA * INITSIGMA);
    else
        dsigma = sqrt(_SIGMA * _SIGMA - INITSIGMA * INITSIGMA);

    Filter::BlurImage(tmpImg, dsigma);

    return tmpImg;
}

float HessAff::InterpKeyStep(int x, int y, int s, vector<Image *> &DI, float *dx, float *dy, float *ds)
{
    unsigned int i = 0;
    float Dp[3] = {0}; // first derivative of D with respect to x, y, s
    float Dpp[3][3] = {{0}}; // HessAff of D
    Dp[0] = (float) ((DI[s]->getPixel(x + 1, y) - DI[s]->getPixel(x - 1, y)) / 2.0); // Dx
    Dp[1] = (float) ((DI[s]->getPixel(x, y + 1) - DI[s]->getPixel(x, y - 1)) / 2.0); // Dy
    Dp[2] = (float) ((DI[s + 1]->getPixel(x, y) - DI[s - 1]->getPixel(x, y)) / 2.0); // Ds

    // Dxx
    Dpp[0][0] = (float) (DI[s]->getPixel(x + 1, y) + DI[s]->getPixel(x - 1, y)
                         - 2.0 * DI[s]->getPixel(x, y));

    // Dyy
    Dpp[1][1] = (float) (DI[s]->getPixel(x, y + 1) + DI[s]->getPixel(x, y - 1)
                         - 2.0 * DI[s]->getPixel(x, y));

    // Dzz
    Dpp[2][2] = (float) (DI[s + 1]->getPixel(x, y) + DI[s - 1]->getPixel(x, y)
                         - 2.0 * DI[s]->getPixel(x, y));


    // Dxy = Dyx
    Dpp[0][1] = Dpp[1][0] = (float) ((DI[s]->getPixel(x + 1, y + 1) - DI[s]->getPixel(x - 1, y + 1)
                                      - DI[s]->getPixel(x + 1, y - 1) + DI[s]->getPixel(x - 1, y - 1)) / 4.0);

    // Dxs = Dsx
    Dpp[0][2] = Dpp[2][0] = (float) ((DI[s + 1]->getPixel(x + 1, y) - DI[s + 1]->getPixel(x - 1, y)
                                      - DI[s - 1]->getPixel(x + 1, y) + DI[s - 1]->getPixel(x - 1, y)) / 4.0);

    // Dys = Dsy
    Dpp[1][2] = Dpp[2][1] = (float) ((DI[s + 1]->getPixel(x, y + 1) - DI[s + 1]->getPixel(x, y - 1)
                                      - DI[s - 1]->getPixel(x, y + 1) + DI[s - 1]->getPixel(x, y - 1)) / 4.0);

    float invDpp[3][3];

    VMath::mInv33(Dpp, invDpp); // Jacobian matrix

    // Solve for delta positions
    *dx = 0;
    for (i = 0; i < 3; i++)
        *dx -= invDpp[0][i] * Dp[i];

    *dy = 0;
    for (i = 0; i < 3; i++)
        *dy -= invDpp[1][i] * Dp[i];

    *ds = 0;
    for (i = 0; i < 3; i++)
        *ds -= invDpp[2][i] * Dp[i];

    //printf("Interp: %f %f %f\n", *dx, *dy, *ds);

    float val = DI[s]->getPixel(x, y);

    // 泰勒展开 Dp(X+∆X) = Dp(X) + Dp'*∆X + 0.5*∆X^T*Dp''*∆X
    // val = Dp(X+∆X) = Dp(X) + Dp'(X)*∆X
    val += 0.5 * (Dp[0] * *ds + Dp[1] * *dy + Dp[2] * *ds);

    return fabs(val);
}

bool
HessAff::InterpKey(int x, int y, int s, vector<Image *> &LoGImages, float *fx, float *fy, float *fs, float *dogVal) {
    bool addkey = true;
    int moves_left = 5;
    int tx = x;
    int ty = y;
    int ts = s;

    float dx, dy, ds, val;
    bool updated;
    float contrast_thresh = (float) (0.8 * this->thresh);

    do {
        moves_left--;
        updated = false;

        val = InterpKeyStep(tx, ty, ts, LoGImages, &dx, &dy, &ds);
        if (useThresh && val < contrast_thresh) {
            addkey = false;
            continue;
        }


        if (dx > 0.6 && tx < LoGImages[0]->width - 3) {
            tx++;
            updated = true;
        } else if (dx < -0.6 && tx > 3) {
            tx--;
            updated = true;
        }


        if (dy > 0.6 && ty < LoGImages[0]->height - 3) {
            ty++;
            updated = true;
        } else if (dy < -0.6 && ty > 3) {
            ty--;
            updated = true;
        }

    } while (moves_left > 0 && updated);
    *dogVal = val;

    if (addkey && fabs(dx) < 1.5 && fabs(dy) < 1.5 && fabs(ds) < 1.5) {
        *fx = tx + dx;
        *fy = ty + dy;
        *fs = ts + ds;
        return true;
    }

    return false;
}


void
HessAff::BuildGaussianPyramid(Image *initImg, vector<vector<Image*> > &gaussianPyramid, const float sigma0, const int numScales)
{

    // octave的总数 log2(size)-2
    int numOctaves = int(log((float) min(initImg->height, initImg->width)) / log(2.0f)) - 2;
    numOctaves = min(numOctaves, MaxOctaves);
    unsigned int j = 0, o = 0, layer = 0;

    // 每个octave中的层数
    const int numLayers = numScales + 2;

    // sigma[o,i] = sigma[o, i-1]*k
    float k = pow(2.0, 1.0 / numScales);

    // dsigma[i] = sqrt(sigma[i]^2-sigma[i-1]^2)
    vector<float> dsigma;
    dsigma.push_back(sigma0);
    float sigTotal, sigPrev = sigma0;
    for (j = 1; j < numLayers; j++)
    {
        sigTotal = sigPrev * k;
        dsigma.push_back(sqrt(sigTotal * sigTotal - sigPrev * sigPrev));
        sigPrev = sigTotal;
    }

    for (o = 0; o < numOctaves; o++)
    {
        // 每个octave底层的图片, 尺度为 sigma0 * 2^o, o=[0..numOctaves-1]
        Image *base = NULL;

        if (o == 0)
        {
            base = initImg;
        }else
        {
            // sigma0*2^o = sigma0*2^(o-1)*k^numLayers = sigma0*2^(o-1)*2
            base = gaussianPyramid.back()[numScales];
            base = Image::halfSizeImage(base);
        }

        vector<Image*> newOct;
        newOct.push_back(base);
        for (layer = 1; layer < numLayers; layer++)
        {
            Image *newImg = new Image(base->width, base->height);
            Filter::BlurImage(newOct.back(), newImg, dsigma[layer]);
            newOct.push_back(newImg);
        }

        gaussianPyramid.push_back(newOct);
    }
}

void HessAff::BuildHessianImage(Image *src, float sigma, Image *&hessImg, Image *&D1LImg)
{
    int width = src->width, height = src->height, x = 0, y = 0;
    hessImg = new Image(width, height);
    D1LImg = new Image(width, height);
    for(x = 0; x < width; x++)
        for(y = 0; y < height; y++)
        {
            float dxx = _dxx(src, x, y);
            float dyy = _dyy(src, x, y);
            float dxy = _dxy(src, x, y);
            float s2 = sigma * sigma;
            float trace = s2 * (dxx + dyy);
            float det = s2 * s2 * (dxx * dyy - dxy * dxy);
            float d1L = det - D1L_ALPHA * trace * trace;
            hessImg->setPixel(x, y, det);
            D1LImg->setPixel(x, y, d1L);
        }
}

void HessAff::BuildHessianPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &HessPyramid,
                                  vector<vector<Image*> > &D1LPyramid)
{
    assert(gaussianPyramid.size() > 0);
    unsigned int i = 0, layer = 0;
    int numLayers = gaussianPyramid.front().size();
    assert(numLayers > 0);

    float sigma;
    Image *crntImg, *hessImg = NULL, *D1LImg = NULL;

    for (i = 0; i < gaussianPyramid.size(); i++)
    {
        vector<Image*> hessOctave, D1LOctave;
        for (layer = 0; layer < numLayers; layer++)
        {
            crntImg = gaussianPyramid[i][layer];
            sigma = _SIGMA * pow(2.0, (float) layer / SCALES);
            BuildHessianImage(crntImg, sigma, hessImg, D1LImg);
            hessOctave.push_back(hessImg);
            D1LOctave.push_back(D1LImg);
        }
        HessPyramid.push_back(hessOctave);
        D1LPyramid.push_back(D1LOctave);
    }
}

bool HessAff::isLocalExtrema(const Image *img, const int x, const int y)
{
    int i = 0, j = 0;
    float val = img->getPixel(x, y);
    bool maximal = true, minimal = true;
    for(i = x - 1; i <= x + 1; i++)
        for(j = y - 1; j <= y + 1; j++)
        {
            if (i == x && j == y)
                continue;
            if (maximal && img->getPixel(i, j) >= val)
                maximal = false;
            if (minimal && img->getPixel(i, j) <= val)
                minimal = false;
            if (not maximal && not minimal)
                return false;
        }

    return true;
}

bool HessAff::isLocalExtrema(const Image *prevImg, const Image *img, const Image *nextImg, const int x, const int y) {
    float val = img->getPixel(x, y);
    const Image *imgs[3] = {prevImg, img, nextImg};
    bool maximal = true, minimal = true;
    int i = 0, j = 0, s = 0;
    for(i = x - 1; i <= x + 1; i++)
        for(j = y - 1; j <= y + 1; j++)
            for (s = 0; s <= 2; s++)
            {
                if (i == x && j == y && s == 1)
                    continue;
                if (maximal && imgs[s]->getPixel(i, j) >= val)
                    maximal = false;
                if (minimal && imgs[s]->getPixel(i, j) <= val)
                    minimal = false;
                if (not maximal && not minimal)
                    return false;
            }
    return true;
}

void HessAff::FindKeypoints(const int octIndex, vector<Image*> &hessOctave, vector<Image*> &D1Loctave, vector<KeyPoint *> &peaks)
{
    Image *kpfound = new Image(hessOctave[0]->width, hessOctave[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0;
    unsigned int s = 0;
    int x = 0, y = 0;

    for (s = 1; s < hessOctave.size() - 1; s++) {

        for (y = BORDER; y < (hessOctave[0]->height - BORDER); y++) {
            for (x = BORDER; x < (hessOctave[0]->width - BORDER); x++) {
                // 除去边缘
                if (useD1L && D1Loctave[s]->getPixel(x, y) <= 0)
                    continue;

                funcVal = hessOctave[s]->getPixel(x, y);
                if (useThresh && funcVal < THRESH)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                // 也可以判断Hessian(det(H))(x,y,s)必须为正定或者负定, 即特征值全部同号, 排除三维空间中的saddle
                if (!isLocalExtrema(hessOctave[s - 1], hessOctave[s], hessOctave[s + 1], x, y))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                if (INTERP_KEYS) {
                    if (!InterpKey(x, y, s, hessOctave, &fx, &fy, &fs, &funcVal))
                        continue;
                }
                // 舍弃fs<=0的特征点, 因为这个点的最优scale很可能在上个octave已经计算并添加了
                if (fs <= 0) continue;

                KeyPoint *peak = new KeyPoint();

                peak->x = (int) round(fx * pow(2.0, octIndex));
                peak->y = (int) round(fy * pow(2.0, octIndex));
                peak->dscale = HessAff::_SIGMA * pow(2.0, octIndex + fs / (float) SCALES);
                peak->octSigma = HessAff::_SIGMA * pow(2.0, fs / (float) SCALES);
                peak->iscale = peak->dscale * HessAff::mag;
                peak->funcVal = funcVal;
                peak->ori = 0;
                peak->scale = s;
                peak->fscale = fs;
                peak->gscale = octIndex + fs / HessAff::SCALES;
                peak->sx = fx;
                peak->sy = fy;
                peak->octIndex = octIndex;
                peak->KP = true;

                leveli_kps.push_back(peak);
                kpfound->setPixel((int) (fx + 0.5), (int) (fy + 0.5), 1);

            }
        }
        peaks.insert(peaks.begin(), leveli_kps.begin(), leveli_kps.end());
        leveli_kps.clear();
    }

    delete kpfound;
}

KeyPoint *HessAff::FindKeypointNearby(Image *hessImg, Image *D1LImg, const int cx, const int cy, const int size)
{
    int radius = size / 2;
    int minDistance = size * size * 2, minX, minY;
    bool found = false;
    int i = 0, j = 0;
    for(i = cx - radius; i < cx + radius; i++)
        for(j = cy - radius; j < cy + radius; j++)
        {
            if (useD1L && D1LImg->getPixel(i, j) <= 0)
                continue;
            if (!isLocalExtrema(hessImg, i, j))
                continue;
            int distance = (i - cx) * (i - cx) + (j - cy) * (j - cy);
            if (distance < minDistance) {
                found = true;
                minDistance = distance;
                minX = i;
                minY = j;
            }
        }
    if (found) {
        KeyPoint *peak = new KeyPoint();
        peak->x = minX;
        peak->y = minY;
        peak->funcVal = hessImg->getPixel(minX, minY);
        return peak;
    } else
        return NULL;
}

Image *HessAff::normalizeWindow(Image *window, float sqrtU[2][2], int cx, int cy)
{
    int width = window->width, height = window->height;
    Image *normWindow = new Image(width, height);
    int i, j = 0;
    for (i = 0; i < width; i++)
        for (j = 0; j < height; j++)
        { // (i,j)为正则化窗口中的坐标
            float dx = i - cx, dy = j - cy;
            float cordX = sqrtU[0][0] * dx + sqrtU[0][1] * dy + cx; // (cordX,cordY)为(i,j)对应的window(未被正则化)中的坐标
            float cordY = sqrtU[1][0] * dx + sqrtU[1][1] * dy + cy;
            // 用(cordX,cordY)附近的四个像素点插值
            float xfloor = floor(cordX), xplus = cordX - xfloor;
            float yfloor = floor(cordY), yplus = cordY - yfloor;
            float pixel = fabs((1 - xplus) * (1 - yplus)) * window->getPixel(xfloor, yfloor)
                          + fabs(xplus * (1 - yplus)) * window->getPixel(xfloor + 1, yfloor)
                          + fabs((1 - xplus) * yplus) * window->getPixel(xfloor, yfloor + 1)
                          + fabs(xplus * yplus) * window->getPixel(xfloor + 1, yfloor + 1);
            normWindow->setPixel(i, j, pixel);
        }
    return normWindow;
}

void
HessAff::cal2ndMomentMat(Image *window, const int x, const int y, const float dscale, const float iscale, float u[2][2]) {
    int width = window->width, height = window->height;
    int i, j = 0;
    Image *Lx2 = new Image(width, height);
    Image *Ly2 = new Image(width, height);
    Image *LxLy = new Image(width, height);

    for (i = 0; i < width; i++)
        for (j = 0; j < height; j++) {
            float dx = dscale * _dx(window, i, j);
            float dy = dscale * _dy(window, i, j);
            Lx2->setPixel(i, j, dx * dx);
            Ly2->setPixel(i, j, dy * dy);
            LxLy->setPixel(i, j, dx * dy);
        }
    Filter::BlurImage(Lx2, iscale);
    Filter::BlurImage(Ly2, iscale);
    Filter::BlurImage(LxLy, iscale);

    u[0][0] = Lx2->getPixel(x, y);
    u[0][1] = LxLy->getPixel(x, y);
    u[1][0] = u[0][1];
    u[1][1] = Ly2->getPixel(x, y);

    delete (Lx2);
    delete (Ly2);
    delete (LxLy);
}

int HessAff::adaptAffine(Image *initImg, KeyPoint *keypoint)
{
    // eg. U = uk*...u3*u2*u1
    float U[2][2] = {{1, 0},
                     {0, 1}};
    float sqrtU[2][2];

    float scaleStep = 1.0 / SCALES / 3;
    const float MAX_SIGMA = HessAff::_SIGMA * pow(2.0, 1.0 + 1.0 / HessAff::SCALES);
    const float MIN_SIGMA = _SIGMA;
    const int width = initImg->width, height = initImg->height;
    const float Qmin = 1.05;

    // affine adaption是否收敛; scale selection是否收敛; keypoint是否发散
    bool affineConverge = false, scaleConverge = false, diverge = false;

    // 当前迭代中, keypoint的坐标(crntX,crntY,sigma) (initImg的尺度为_SIGMA)
    int crntX = (int) round(keypoint->sx), crntY = (int) round(keypoint->sy);
    float sigma = keypoint->octSigma;
    float H = keypoint->funcVal, dHds;   // det(Hessian(cx,cy,sigma)); ∂det(H)/∂sigma

    const int minimalWindowSize = 9; // minimalWindowSize*2+1
    int brx;                         // half of the window size
    int cx, cy;                      // center of the window

    int iter = 0;    // iteration times
    int v, vold = 0; // optimization direction
    //cout<<"bug 1.0\n";

    // both scale and affine converge
    while (!(scaleConverge && affineConverge))
    {
        iter++;
        if (iter > HessAff::maxAdaptIter)
        {
            diverge = true;
            break;
        }

        Inverse2D(U);           // U = u1^-1*u2^-1...
        Normalize(U);

        if (!SqrtMatrix2D(U, sqrtU))
        { // sqrtU = u1^-0.5*u2^-0.5...
            diverge = true;
            break;
        }

        // ratio of two eigenvalues
        float a, b;
        Eigenvalue2D(U, a, b);
        float Q = (abs(max(a, b))) / (abs(min(a, b)));

        //cropping window from original image
        brx = max(minimalWindowSize, (int) round(Q * 4 * sigma));
        Image *window = HessAff::ImageCut(initImg, max(0, crntX - brx), min(width - 1, crntX + brx),
                                          max(0, crntY - brx), min(height - 1, crntY + brx));
        cx = min(brx, crntX);
        cy = min(brx, crntY);

        // shrink the window size by window  = U^0.5*window
        window = HessAff::normalizeWindow(window, sqrtU, cx, cy);
        brx = max(minimalWindowSize, (int) round(4 * sigma));
        window = HessAff::ImageCut(window, max(0, cx - brx), min(window->width - 1, cx + brx),
                                   max(0, cy - brx), min(window->height - 1, cy + brx));
        cx = min(brx, crntX);
        cy = min(brx, crntY);

        // the window size is sigma
        Filter::BlurImage(window, sqrt(sigma * sigma - _SIGMA * _SIGMA));

        // get the direction of ∂H/∂s
        dHds = getdHds(window, cx, cy, sigma);
        H = getDetH(window, cx, cy, sigma);
        float ds;    // rate of variation of sigma

        // window(x,y,sigma),(cx,cy), (x,y) is the extreme point, grad=0
        if (fabs(dHds) > EPSILON)
        {
            v = HessAff::Sign(H) * HessAff::Sign(dHds);
            if (v * vold < 0)  // change in the direction of scale optimization
                scaleStep = scaleStep / 2;
            vold = v;

            // check if the new scale is better otherwise reduce the step,
            float newSigma = sigma * pow(2.0, scaleStep * v);

            if (newSigma > MAX_SIGMA || newSigma < MIN_SIGMA)
            {
                diverge = true;
                break;
            }

            Image *newWindow = new Image(window->width, window->height);
            Filter::BlurImage(window, newWindow, sqrt(newSigma * newSigma - _SIGMA * _SIGMA));
            double newH = getDetH(newWindow, cx, cy, newSigma);
            if (H * H < newH * newH)
            {
                window = newWindow;
                ds = newSigma - sigma;
                sigma = newSigma;
                scaleConverge = false;
            } else
            {
                scaleConverge = true;
                ds = 0;
            }
        } else
        {
            scaleConverge = true;
            ds = 0;
        }


        if (not(scaleConverge && affineConverge))
        {
            Image *hessImg = NULL, *D1LImg = NULL;
            BuildHessianImage(window, sigma, hessImg, D1LImg);

            int newCx, newCy;// under new scale, (cx,cy) is the closest point to extreme point

            // following feature trajectory to calculate dx/ds dy/ds
            // <Scale-Space Behaviour of Local Extrema and Blobs>-Lindeberg
            float dxds, dyds, dx, dy;
            if (ds != 0)
            {
                getdRds(window, cx, cy, sigma, dxds, dyds);
                dx = ds * dxds;
                dy = ds * dyds;
                newCx = (int) round(cx + ds * dxds);
                newCy = (int) round(cy + ds * dyds);
            } else
            {
                // when there is no scale transformation, det(H) attains loal maxima at (cx, cy) according to affine covariance
                newCx = cx;
                newCy = cy;
                dx = dy = 0;
            }
            // TODO 多小? 会影响到收敛率
            // 理论上来说(newCx, newCy)是(x,y)空间中的极值点了,但是因为像素是离散的, 还是应该在小窗口中找最近的极值点
            // 窗口越小，速度越快，但是收敛率也越低
            // 直观上来说, 窗口的大小应该跟三维空间中距离变化大小成正比

            int windowSize;
            //windowSize = 10;
//            windowSize = 5;
//            windowSize = (int) (sqrt(dx * dx + dy * dy) * 2 + 2);
//            windowSize = (int) (sqrt(dx * dx + dy * dy) * 2 + 4);
//            windowSize = (int) (sqrt(dx * dx + dy * dy) * 3 + 4);
//            windowSize = (int) (sqrt(dx * dx + dy * dy) * 4 + 4);
//            windowSize = (int) (sqrt(dx * dx + dy * dy) * 4 + 6);
            windowSize = (int) (2 * sqrt(dx * dx + dy * dy + ds * ds) + 4);
//            windowSize = (int) (3 * sqrt(dx * dx + dy * dy + ds * ds) + 4);
//            windowSize = (int) (3 * sqrt(dx * dx + dy * dy + ds * ds) + 6);


            assert(windowSize < min(window->width, window->height));

            KeyPoint *nearbyKeypoint = FindKeypointNearby(hessImg, D1LImg, newCx, newCy, windowSize);

            // theoritically speaking, nearbyKeypoint is (newCx, newCy)
            if (nearbyKeypoint != NULL)
            {
                H = nearbyKeypoint->funcVal;
                newCx = nearbyKeypoint->x;
                newCy = nearbyKeypoint->y;
                if (useThresh && H < THRESH)
                {
                    diverge = true;
                    break;
                }

                float dx = newCx - cx, dy = newCy - cy;
                crntX = (int) round(sqrtU[0][0] * dx + sqrtU[0][1] * dy + cx) + max(0, crntX - brx);
                crntY = (int) round(sqrtU[1][0] * dx + sqrtU[1][1] * dy + cy) + max(0, crntY - brx);

                float u[2][2]; // second-moment matrix
                cal2ndMomentMat(window, newCx, newCy, sigma, sigma * HessAff::mag, u);

                if (!Inverse2D(u)) { // u = u^-1
                    diverge = true;
                    break;
                }
                Normalize(u);

                MultiplyMatrix(U, u); // Unew = U*u = u1^-1*u2^-1...ui^-1
                Normalize(U);

                // check whether max(eig)/min(eig)
                float ua, ub;
                Eigenvalue2D(u, ua, ub);

                if (max(ua, ub) / min(ua, ub) < Qmin)
                    affineConverge = true;
                else
                    affineConverge = false;
                Inverse2D(U); // U = ui*ui-1*....u3*u2*u1
                Normalize(U);
            } else
            {
                diverge = true; // failed to find the extreme point
                break;
            }
        }

    }

    if (!(scaleConverge && affineConverge))
        diverge = true;

    if (!diverge)
    {
        keypoint->sx = crntX;
        keypoint->sy = crntY;
        keypoint->x = (int) round(crntX * pow(2.0, keypoint->octIndex));
        keypoint->y = (int) round(crntY * pow(2.0, keypoint->octIndex));
        keypoint->octSigma = sigma;
        keypoint->dscale = sigma * pow(2.0, keypoint->octIndex);
        keypoint->funcVal = H;
        keypoint->u1 = U[0][0];
        keypoint->u2 = U[0][1];
        keypoint->u3 = U[1][0];
        keypoint->u4 = U[1][1];
    }
    else
      iter = 0;

    return iter;
}

float HessAff::getdHds(const Image *window, const int cx, const int cy, const float sigma) {

    float t = sigma * sigma;
    float dxx, dyy, dxy, dxxxx, dyyyy, dxxyy, dxxxy, dxyyy;
    float xn = window->getPixel(cx - 1, cy); // x negative 1
    float xp = window->getPixel(cx + 1, cy); // x positive 1
    float xy = window->getPixel(cx, cy);
    float yn = window->getPixel(cx, cy - 1);
    float yp = window->getPixel(cx, cy + 1);
    float xpyp = window->getPixel(cx + 1, cy + 1);
    float xnyn = window->getPixel(cx - 1, cy - 1);
    float xnyp = window->getPixel(cx - 1, cy + 1);
    float xpyn = window->getPixel(cx + 1, cy - 1);

    dxx = xn + xp - 2 * xy;
    dyy = yn + yp - 2 * xy;
    dxy = (xpyp + xnyn - xnyp - xpyn) / 4;
    dxxxx = window->getPixel(cx - 2, cy) + window->getPixel(cx + 2, cy) - 4 * (xn + xp) + 6 * xy;
    dyyyy = window->getPixel(cx, cy - 2) + window->getPixel(cx, cy + 2) - 4 * (yn + yp) + 6 * xy;
    dxxyy = xpyp + xpyn + xnyp + xnyn - 2 * (xn + xp + yn + yp) + 4 * xy;
    dxxxy = (window->getPixel(cx - 2, cy - 1) + window->getPixel(cx + 2, cy + 1)
             - window->getPixel(cx - 2, cy + 1) - window->getPixel(cx + 2, cy - 1)) / 4
            + (xnyp + xpyn - xnyn - xpyp) / 2;
    dxyyy = (window->getPixel(cx - 1, cy - 2) + window->getPixel(cx + 1, cy + 2)
             - window->getPixel(cx - 1, cy + 2) - window->getPixel(cx + 1, cy - 2)) / 4
            + (xnyp + xpyn - xnyn - xpyp) / 2;
    float dhs = ((dxxxx + dxxyy) * dyy + (dxxyy + dyyyy) * dxx) / 2 - dxy * (dxxxy + dxyyy);
    return 2 * t * (dxx * dyy - dxy * dxy) + t * t * dhs;
}

float HessAff::getDetH(const Image *window, const int cx, int cy, float sigma)
{
    float dxx = _dxx(window, cx, cy);
    float dyy = _dyy(window, cx, cy);
    float dxy = _dxy(window, cx, cy);
    return (float) pow(sigma, 4.0) * (dxx * dyy - dxy * dxy);
}

bool HessAff::getdRds(const Image *window, const int cx, const int cy, const float sigma, float &dxds, float &dyds)
{
    float dxx = _dxx(window, cx, cy);
    float dyy = _dyy(window, cx, cy);
    float dxy = _dxy(window, cx, cy);
    float det = dxx * dyy - dxy * dxy;

    if (fabs(det) < EPSILON)
        return false;

    // dRds = dRdt * dtds = 2s*dRdt
    float xn = window->getPixel(cx - 1, cy); // x negative 1
    float xp = window->getPixel(cx + 1, cy); // x positive 1
    float yn = window->getPixel(cx, cy - 1);
    float yp = window->getPixel(cx, cy + 1);
    float xpyp = window->getPixel(cx + 1, cy + 1);
    float xnyn = window->getPixel(cx - 1, cy - 1);
    float xnyp = window->getPixel(cx - 1, cy + 1);
    float xpyn = window->getPixel(cx + 1, cy - 1);
    float v = -sigma / det;

    float dxxx, dxxy, dxyy, dyyy;
    float xnn = window->getPixel(cx - 2, cy);
    float xpp = window->getPixel(cx + 2, cy);
    float ynn = window->getPixel(cx, cy - 2);
    float ypp = window->getPixel(cx, cy + 2);
    dxxx = (xpp - 2 * xp + 2 * xn - xnn) / 2;
    dyyy = (ypp - 2 * yp + 2 * yn - ynn) / 2;
    dxyy = (2 * (xn - xp) + xpyn + xpyp - xnyn - xnyp) / 2;
    dxxy = (2 * (yn - yp) + xnyp + xpyp - xnyn - xpyn) / 2;

    dxds = v * (dyy * dxxx + dyy * dxyy - dxy * dxxy - dxy * dyyy);
    dyds = v * (dxx * dyyy + dxx * dxxy - dxy * dxxx - dxy * dxyy);

    return true;
}

bool HessAff::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    clock_t start = clock();
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);

    Image rawImage(fn);

    if (!rawImage.isActive())
        return false;

    // ZOOM_OUT;
    Image *initImg = CreateInitialImage(&rawImage);
    this->crntimg = initImg;

    //Gaussian Pyramid
    vector<vector<Image*> > gaussianPyramid;
    BuildGaussianPyramid(initImg, gaussianPyramid, HessAff::_SIGMA, HessAff::SCALES);

    //Hessian Pyramid
    vector<vector<Image*> > hessPyramid, D1LPyramid;
    BuildHessianPyramid(gaussianPyramid, hessPyramid, D1LPyramid);

    //detect initial points' locations
    vector<KeyPoint *> initKeypoints, peaks;
    vector<KeyPoint *>::iterator vit;
    KeyPoint *crntPt = NULL;
    for (int ioctave = 0; ioctave < hessPyramid.size(); ioctave++)
    {
        FindKeypoints(ioctave, hessPyramid[ioctave], D1LPyramid[ioctave], peaks);
        initKeypoints.insert(initKeypoints.end(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    clock_t finish = clock();
    cout << "Time costs: " << (float) (finish - start) / CLOCKS_PER_SEC << endl;
    start = clock();

    cout << initKeypoints.size() << endl;

    int index = 0, sumIter = 0;
    for (vit = initKeypoints.begin(); vit != initKeypoints.end(); vit++)
    {
        index++;
        crntPt = *vit;
        cout << "adapting: " << index << "/" << initKeypoints.size() << "\n";
        int iter = adaptAffine(gaussianPyramid[crntPt->octIndex][0], crntPt);
        if (iter != 0)
        {
            this->kps.push_back(crntPt);
            sumIter += iter;
        }
    }
    printf("Num. of iteration: %.4f\n", sumIter / (float) this->kps.size());

    finish = clock();
    cout << "adaptation time:" << (float) (finish - start) / CLOCKS_PER_SEC << endl;

    cout << this->kps.size() << endl;
    cout << "convergence rate: " << ((float) this->kps.size()) / initKeypoints.size() << endl;

//    SaveKeypoints(this->crntimg, fn, this->kps);
    initKeypoints.clear();

    /**
    vector<KeyPoint *> extra_peaks = this->FindOrientByGrad(kps, gaussianPyramid);

    if (extra_peaks.size() > 0)
    {
        this->kps.insert(kps.begin(), extra_peaks.begin(), extra_peaks.end());
        extra_peaks.clear();
    }
    **/

    switch (this->sel_option)
    {
        case 0: {
            AbstractDetector::topkSelect(kps, this->fix_kp_numb);
            break;
        }
        case 1: {
            break;
        }
        default: {
            //AbstractDetector::topkSelect(kps,this->fix_kp_numb);
            break;
        }
    }
    cout<<"I am here\n";
// TODO
    if (strcmp(dstfn, ""))
    {
        ///(this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
        this->saveKeypoints(dstfn, kps);
    }

    //cout<<"i am out\n";

    if (strcmp(descfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }

    if (strcmp(dvfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }
//    Cleaner::releaseOctaves(hessPyramid);
//    Cleaner::releaseOctaves(gaussianPyramid);
    cout<<"i am done\n";

    delete this->intImg;
    delete this->crntimg;
    this->intImg = NULL;
    this->crntimg = NULL;
    return true;
}

void HessAff::extractKeyp(string imgPath)
{
    int i = 0, j = 0;
    cout<<imgPath<<endl;
    Image *testImg = new Image(&imgPath[0]);
    vector<vector<Image*> > gaussianPyramid;
    BuildGaussianPyramid(testImg, gaussianPyramid, HessAff::_SIGMA, HessAff::SCALES);
    vector<vector<Image*> > hessPyramid, D1LPyramid;
    BuildHessianPyramid(gaussianPyramid, hessPyramid, D1LPyramid);
//
//    Image *dv2 = Image::halfSizeImage(testImg);
//
//    string filename("exm1div2.png");
//    dv2->save(&filename[0]);
    cout<<"bug 1.0\n";

    vector<KeyPoint *> initKeypoints;
    vector<KeyPoint *> peaks;
    for (int ioctave = 0; ioctave < hessPyramid.size(); ioctave++)
    {
        FindKeypoints(ioctave, hessPyramid[ioctave], D1LPyramid[ioctave], peaks);
        initKeypoints.insert(initKeypoints.end(), peaks.begin(), peaks.end());
        peaks.clear();
    }
    ofstream output("/home/wlzhao/detH.txt");
    cout<<"bug 1.1\n";

    for (int k = 0; k < initKeypoints.size(); k++)
    {
        KeyPoint *sample = initKeypoints[k];
        int x = sample->x, y = sample->y;
        output << "cord: (" << x << "," << y << ")\n";
        output << "point: " << sample->dscale << "    " << sample->funcVal << endl;

        for (i = 0; i < hessPyramid.size(); i++)
        {
            for (j = 1; j < hessPyramid[i].size() - 1; j++)
            {
                output << _SIGMA * pow(2.0, i + (float) j / SCALES)
                       << "   " << hessPyramid[i][j]->getPixel(x, y) << "\n";
            }
            x = x / 2;
            y = y / 2;
        }
        output << "\n";
    }
    cout<<"bug 1.2\n";
    output.close();

}

vector<KeyPoint *> HessAff::FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > &GOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint *> newkps;
    KeyPoint *newkp = NULL;
    int indexa = 0, indexb = 0, indexc = 0, j = 0, c = 0, index = 0, x = 0, y = 0;
    float sigma = 0.0f, m = 0.0f, theta = 0.0f, degree = 0.0f;
    float maxval = 0.0f, maxp = 0.0f, weight = 0, Wi = 1.0f;
    float thetaa = 0.0f, thetab = 0.0f, thetac = 0.0f;
    float thetas[NOrient] = {0};
    unsigned int i = 0;
    bool valid = false;

    for (i = 0; i < kps.size(); i++) {
        sigma = 1.5 * pow(2.0, (kps[i]->fscale) / (float) SCALES) * HessAff::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size() / 2;
        Wi = 3.0f * sigma;
        Wi = Wi > 1.0f ? Wi : 1.0f;
        Wi = Wi * Wi;

        memset(thetas, 0, sizeof(float) * NOrient);
        for (y = -c; y <= c; y++) {
            for (x = -c; x <= c; x++) {
                if ((x * x + y * y) > Wi)
                    continue;

                valid = Filter::GetPixOrientation((int) (kps[i]->sx + x + 0.5), (int) (kps[i]->sy + y + 0.5),
                                                  GOctaves[kps[i]->octIndex][kps[i]->scale], m, theta);

                if (valid) {
                    degree = theta / PI * 180.0 + 180.0;
                    index = ((int) (degree / DEGREE));
                    index = index % NOrient;
                    weight = m * gmat[y + c][x + c];
                    thetas[index] += weight;
                }
            }
        }
        vector<vector<float> >::iterator
                it;
        vector<float> crntvect;
        for (it = gmat.begin(); it != gmat.end(); ++it) {
            crntvect = *it;
            crntvect.clear();
        }
        gmat.erase(gmat.begin(), gmat.end());

        for (j = 0; j < 6; j++) {
            VMath::SmoothHistogram(thetas, NOrient);
        }

        maxval = VMath::maxVec(thetas, NOrient);

        for (j = 0; j < NOrient; j++) {
            if (thetas[j] < maxval * NwKpThresh)
                continue;

            indexa = VMath::mod(j - 1, NOrient);
            indexb = j;
            indexc = VMath::mod(j + 1, NOrient);
            thetaa = thetas[indexa];
            thetab = thetas[indexb];
            thetac = thetas[indexc];

            if (!(thetab > thetaa && thetab > thetac))
                continue;

            maxp = VMath::getMaxParabola(-1, thetaa, 0, thetab, 1, thetac);

            if (thetas[j] == maxval) {
                kps[i]->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;
            } else {
                newkp = new KeyPoint();
                memcpy(newkp, kps[i], sizeof(KeyPoint));
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;
                newkps.push_back(newkp);
            }
        }
    }
    return newkps;
}

void HessAff::test()
{
     const char *srcFn = "/home/wlzhao/datasets/vgg/graf/graf1.jpg";
     const char *dstFn = "/home/wlzhao/datasets/vgg/graf1.keys";
     HessAff *myhes = new HessAff();
     ///myhes->extractKeyp(srcFn);
     myhes->KeypointBuild(srcFn, dstFn, "", "");

}

