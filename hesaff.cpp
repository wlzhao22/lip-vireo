#include "hesaff.h"

#include "cleaner.h"
#include "filter.h"
#include "board.h"
#include "vmath.h"


const int HesAff::MaxOctaves  = 10;
const int HesAff::SCALES      = 7;
const float HesAff::_SIGMA    = 1.4;
const float HesAff::INITSIGMA = 0.5;
const bool  HesAff::INTERP_KEYS = false; // TODO true?

bool HesAff::ZOOM_OUT         = false;
const float HesAff::mag       = 4.0; //optimal 8
const int HesAff::BORDER      = 5;
const int HesAff::THRESH      = 100;
const float HesAff::D1L_ALPHA = 0.05; // TODO D1L_ALPHA

const int HesAff::DEGREE      = 10;
const int HesAff::NOrient     = 36;
const int HesAff::DEGPERBIN   = (360 / NOrient);
const float HesAff::NwKpThresh = 0.8;

const int HesAff::maxAdaptIter = 30;
const float HesAff::EPSILON    = 1e-6;
bool HesAff::useD1L            = false;
bool HesAff::useThresh         = true;

#ifndef _INTMAT_
#define _INTMAT_
const float IntMat[9][9] = {
{0.000023, 0.000138, 0.000495, 0.001064, 0.001373, 0.001064, 0.000495, 0.000138, 0.000023},
{0.000138, 0.000825, 0.002953, 0.006347, 0.008191, 0.006347, 0.002953, 0.000825, 0.000138},
{0.000495, 0.002953, 0.010572, 0.022725, 0.029329, 0.022725, 0.010572, 0.002953, 0.000495},
{0.001064, 0.006347, 0.022725, 0.048852, 0.063048, 0.048852, 0.022725, 0.006347, 0.001064},
{0.001373, 0.008191, 0.029329, 0.063048, 0.081369, 0.063048, 0.029329, 0.008191, 0.001373},
{0.001064, 0.006347, 0.022725, 0.048852, 0.063048, 0.048852, 0.022725, 0.006347, 0.001064},
{0.000495, 0.002953, 0.010572, 0.022725, 0.029329, 0.022725, 0.010572, 0.002953, 0.000495},
{0.000138, 0.000825, 0.002953, 0.006347, 0.008191, 0.006347, 0.002953, 0.000825, 0.000138},
{0.000023, 0.000138, 0.000495, 0.001064, 0.001373, 0.001064, 0.000495, 0.000138, 0.000023}};
#endif // _INTMAT_
HesAff::HesAff()
{
    cout << "Detector ................................ HesAff\n";
    this->DETECTOR   = hesaff;
    this->intImg     = NULL;
    this->sel_option = THRSH;
    this->AFF_OUT    = true;
}

bool HesAff::paramsCheck()
{
    const char *argv[] = {"sigma", "thresh", "topk", "dens"};

    if (this->paras.find(argv[1]) != this->paras.end())
    {
        this->thresh = (float) atof(paras["thresh"]);
    }
    else
    {
        this->thresh = THRESH;
    }
    cout << "Thresh .................................. " << this->thresh << endl;

    if (this->paras.find(argv[2]) != this->paras.end())
    {
        this->fix_kp_numb = atoi(paras["topk"]);
        this->sel_option = TOPK;
        cout << "Topk .................................... " << this->fix_kp_numb << endl;
    }

    if (this->paras.find(argv[3]) != this->paras.end())
    {
        this->fix_kp_numb = atoi(paras["dens"]);
        this->sel_option = DENS;
        cout << "Topk .................................... " << this->fix_kp_numb << endl;
    }

    return true;
}

Image *HesAff::CreateInitialImage(Image *crntImg)
{
    Image *tmpImg;
    int minsz = min(crntImg->width, crntImg->height);
    vector<Image *> initImgs;

    if (minsz <= 200)
    {
        tmpImg = Image::doubleSizeImage(crntImg);
        this->resize_rate = 0.5;
        this->RESIZE = true;
        ZOOM_OUT = true;
    }
    else
    {
        this->resize_rate = 1.0f;
        ZOOM_OUT = false;
        tmpImg = crntImg->clone();
        while (minsz > imgSzBound)
        {
            tmpImg = Image::halfSizeImage(tmpImg);
            minsz = min(tmpImg->width, tmpImg->height);
            this->resize_rate = this->resize_rate * 2;
            this->RESIZE = true;
        }
    }

    float dsigma;
    if (ZOOM_OUT)
        dsigma = sqrt(_SIGMA * _SIGMA - 4 * INITSIGMA * INITSIGMA);
    else
        dsigma = sqrt(_SIGMA * _SIGMA - INITSIGMA * INITSIGMA);

    Filter::BlurImage(tmpImg, dsigma);

    return tmpImg;
}


void HesAff::BuildGaussianPyramid(Image *initImg, vector<vector<Image*> > &gaussianPyramid, const float sigma0, const int numScales)
{
    // octave的总数 log2(size)-2
    unsigned numOctaves = int(log((float) min(initImg->height, initImg->width)) / log(2.0f)) - 2;
    numOctaves = min((int)numOctaves, MaxOctaves);
    unsigned int j = 0, o = 0, layer = 0;
    const unsigned numLayers = numScales + 2;
    float k = pow(2.0, 1.0 / numScales);
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
        Image *base = NULL;
        if (o == 0)
        {
            base = initImg->clone();
        }
        else
        {
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

void HesAff::BuildHessianImage(Image *src, float sigma, Image *&hessImg, Image *&D1LImg)
{
    float dxx = 0, dyy = 0, dxy = 0, s2 = 0, trace = 0, det = 0, d1L = 0;
    int width = src->width, height = src->height, x = 0, y = 0;
    hessImg = new Image(width, height);
    D1LImg  = new Image(width, height);
    //DxxImg  = new Image(width, height);
    //DxyImg  = new Image(width, height);
    //DyyImg  = new Image(width, height);
    for(x = 0; x < width; x++)
    {
        for(y = 0; y < height; y++)
        {
            dxx = _Dxx(src, x, y);
            dyy = _Dyy(src, x, y);
            dxy = _Dxy(src, x, y);
            s2  = sigma * sigma;
            trace = s2 * (dxx + dyy);
            det = s2 * (dxx * dyy - dxy * dxy);
            ///det = _sigma*_sigma*(dxx*dyy - dxy*dxy);
            d1L = det - D1L_ALPHA * trace * trace;
            hessImg->setPixel(x, y, det);
            D1LImg->setPixel(x,  y, d1L);
            //DxxImg->setPixel(x, y, dxx);
            //DxyImg->setPixel(x, y, dxy);
            //DyyImg->setPixel(x, y, dyy);
        }
    }

}

void HesAff::BuildHessianPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &HessPyramid,
                                 vector<vector<Image*> > &D1LPyramid)
{
    assert(gaussianPyramid.size() > 0);
    unsigned int i = 0, layer = 0;
    unsigned numLayers = gaussianPyramid.front().size();
    assert(numLayers > 0);
    ///cout<<"numb: "<<numLayers<<endl;

    float sigma;
    Image *crntImg, *hessImg = NULL, *D1LImg = NULL;
    Image *DxxImg = NULL, *DxyImg = NULL, *DyyImg = NULL;

    for (i = 0; i < gaussianPyramid.size(); i++)
    {
        vector<Image*> hessOctave, D1LOctave;
        vector<Image*> dxxOctave, dxyOctave, dyyOctave;
        for (layer = 0; layer < numLayers; layer++)
        {
            crntImg = gaussianPyramid[i][layer];
            sigma   = _SIGMA * pow(2.0, (float) layer / SCALES);
            BuildHessianImage(crntImg, sigma, hessImg, D1LImg);
            hessOctave.push_back(hessImg);
            D1LOctave.push_back(D1LImg);
            //dxxOctave.push_back(DxxImg);
            //dxyOctave.push_back(DxyImg);
            //dyyOctave.push_back(DyyImg);
        }
        HessPyramid.push_back(hessOctave);
        D1LPyramid.push_back(D1LOctave);
        //dxxPyramid.push_back(dxxOctave);
        //dxyPyramid.push_back(dxyOctave);
        //dyyPyramid.push_back(dyyOctave);
    }
}

void HesAff::buildDxPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &DxPyramid)
{
    assert(gaussianPyramid.size() > 0);
    unsigned int i = 0, layer = 0;
    unsigned numLayers = gaussianPyramid.front().size();
    assert(numLayers > 0);
    Image *crntImg = NULL, *DxImg = NULL;

    for (i = 0; i < gaussianPyramid.size(); i++)
    {
        vector<Image*> dxOctave;
        for (layer = 0; layer < numLayers; layer++)
        {
            crntImg = gaussianPyramid[i][layer];
            DxImg   = new Image(crntImg->width, crntImg->height);
            buildDxImage(crntImg, DxImg);
            dxOctave.push_back(DxImg);
        }
        DxPyramid.push_back(dxOctave);
    }
}

void HesAff::buildDyPyramid(vector<vector<Image*> > &gaussianPyramid, vector<vector<Image*> > &DyPyramid)
{
    assert(gaussianPyramid.size() > 0);
    unsigned int i = 0, layer = 0;
    unsigned numLayers = gaussianPyramid.front().size();
    assert(numLayers > 0);
    Image *crntImg = NULL, *DyImg = NULL;

    for (i = 0; i < gaussianPyramid.size(); i++)
    {
        vector<Image*> dyOctave;
        for (layer = 0; layer < numLayers; layer++)
        {
            crntImg = gaussianPyramid[i][layer];
            DyImg   = new Image(crntImg->width, crntImg->height);
            buildDyImage(crntImg, DyImg);
            dyOctave.push_back(DyImg);
        }
        DyPyramid.push_back(dyOctave);
    }
}

void HesAff::buildDxImage(Image *src, Image *dxImg)
{
    float dx = 0;
    int width = src->width, height = src->height, x = 0, y = 0;
    for(x = 0; x < width; x++)
    {
        for(y = 0; y < height; y++)
        {
            dx = _Dx(src, x, y);
            dxImg->setPixel(x, y, dx);
        }
    }
}

void HesAff::buildDyImage(Image *src, Image *dyImg)
{
    float dy = 0;
    int width = src->width, height = src->height, x = 0, y = 0;
    for(x = 0; x < width; x++)
    {
        for(y = 0; y < height; y++)
        {
            dy = _Dy(src, x, y);
            dyImg->setPixel(x, y, dy);
        }
    }
}

bool HesAff::isLocalExtrema(const Image *img, const int x, const int y)
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
            if (!maximal && !minimal)
                return false;
        }

    return true;
}

bool HesAff::isLocalExtrema(const Image *prevImg, const Image *img, const Image *nextImg, const int x, const int y)
{
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
                if (!maximal && !minimal)
                    return false;
            }
    return true;
}

void HesAff::FindKeypoints(const int octIndex, vector<Image*> &hessOctave, vector<Image*> &D1Loctave, vector<KeyPoint *> &peaks)
{
    Image *kpfound = new Image(hessOctave[0]->width, hessOctave[0]->height);
    float fx = 0, fy = 0, fs = 0, funcVal = 0;
    unsigned int s = 0;
    int x = 0, y = 0;

    for (s = 1; s < hessOctave.size() - 1; s++)
    {

        for (y = BORDER; y < (hessOctave[0]->height - BORDER); y++)
        {
            for (x = BORDER; x < (hessOctave[0]->width - BORDER); x++)
            {

                if (useD1L && D1Loctave[s]->getPixel(x, y) <= 0)
                    continue;

                funcVal = hessOctave[s]->getPixel(x, y);
                if (useThresh && funcVal < THRESH)
                    continue;

                if (kpfound->getPixel(x, y) == 1)
                    continue;

                /// 也可以判断Hessian(det(H))(x,y,s)必须为正定或者负定, 即特征值全部同号, 排除三维空间中的saddle
                if (!isLocalExtrema(hessOctave[s - 1], hessOctave[s], hessOctave[s + 1], x, y))
                    continue;

                fx = x;
                fy = y;
                fs = s;

                if (INTERP_KEYS)
                {
                    if (!InterpKey(x, y, s, hessOctave, &fx, &fy, &fs, &funcVal))
                        continue;
                }
                // 舍弃fs<=0的特征点, 因为这个点的最优scale很可能在上个octave已经计算并添加了
                if (fs <= 0) continue;

                KeyPoint *peak = new KeyPoint();

                peak->x = (int) round(fx * pow(2.0, octIndex));
                peak->y = (int) round(fy * pow(2.0, octIndex));
                peak->dscale = HesAff::_SIGMA * pow(2.0, octIndex + fs / (float) SCALES);
                peak->octSigma = HesAff::_SIGMA * pow(2.0, fs / (float) SCALES);
                peak->iscale  = peak->dscale * HesAff::mag;
                peak->funcVal = funcVal;
                peak->ori    = 0;
                peak->scale  = s;
                peak->fscale = fs;
                peak->gscale = octIndex + fs / HesAff::SCALES;
                peak->sx     = fx;
                peak->sy     = fy;
                peak->octIndex = octIndex;
                //if(peak->iscale < 3.0)
                //    peak->KP = false;
                //else
                //    peak->KP     = true;

                leveli_kps.push_back(peak);
                kpfound->setPixel((int) (fx + 0.5), (int) (fy + 0.5), 1);

            }
        }
        peaks.insert(peaks.begin(), leveli_kps.begin(), leveli_kps.end());
        leveli_kps.clear();
    }

    delete kpfound;
}

bool HesAff::KeypointBuild(const char *fn, const char *dstfn, const char *descfn, const char *dvfn)
{
    clock_t start = clock();
    assert(fn);
    AbstractDetector::releaseKpList(this->kps);
    unsigned int ioctave = 0;
    this->AFF_OUT = true;

    Image rawImage(fn);

    if (!rawImage.isActive())
        return false;

    /// ZOOM_OUT;
    Image *initImg = CreateInitialImage(&rawImage);
    this->crntimg = initImg;

    ///Gaussian Pyramid
    vector<vector<Image*> > hessPyramid, D1LPyramid;
    vector<vector<Image*> > dxPyramid, dyPyramid;
    vector<vector<Image*> > gaussianPyramid;

    BuildGaussianPyramid(initImg, gaussianPyramid, HesAff::_SIGMA, HesAff::SCALES);
    BuildHessianPyramid(gaussianPyramid, hessPyramid, D1LPyramid);

    ///detect initial points' locations
    vector<KeyPoint *> initKeypoints, peaks;
    vector<KeyPoint *>::iterator vit;

    for (ioctave = 0; ioctave < hessPyramid.size(); ioctave++)
    {
        FindKeypoints(ioctave, hessPyramid[ioctave], D1LPyramid[ioctave], peaks);
        initKeypoints.insert(initKeypoints.end(), peaks.begin(), peaks.end());
        peaks.clear();
    }

    buildDxPyramid(gaussianPyramid, dxPyramid);
    buildDyPyramid(gaussianPyramid, dyPyramid);

    kps.insert(kps.begin(), initKeypoints.begin(), initKeypoints.end());
    initKeypoints.clear();

    switch (this->sel_option)
    {
    case 0:
    {
        AbstractDetector::topkSelect(kps, this->fix_kp_numb);
        break;
    }
    case 1:
    {
        break;
    }
    default:
    {
        break;
    }
    }

    ///clock_t t2, t1 = clock();
    ///affineAdapt(kps, gaussianPyramid, dxPyramid, dyPyramid);
    ///affineAdapt2D(kps, gaussianPyramid);
    ///t2 = clock();
    ///float t3 = (t2 - t1+0.0)/CLOCKS_PER_SEC;
    ///cout<<"\ttime 1: "<<t3<<endl;
    ///FindOrientByGrad(kps, gaussianPyramid);
    ///float t4 = (clock() - t2+0.0)/CLOCKS_PER_SEC;
    ///cout<<"\ttime 2: "<<t4<<endl;

    Cleaner::releaseOctaves(hessPyramid);
    Cleaner::releaseOctaves(D1LPyramid);
    Cleaner::releaseOctaves(dxPyramid);
    Cleaner::releaseOctaves(dyPyramid);
    Cleaner::releaseOctaves(gaussianPyramid);
    //Cleaner::releaseOctaves(dxxPyramid);
    //Cleaner::releaseOctaves(dyyPyramid);
    //Cleaner::releaseOctaves(dxyPyramid);

    if(strcmp(dstfn, ""))
    {
        ///cout<<"hi: \t"<<this->AFF_OUT<<endl;
        (this->*saveKpts)(this->kps, this->kps.size(), dstfn, this->resize_rate, this->AFF_OUT);
        ///this->saveKeypoints(dstfn, kps);
    }


    if (strcmp(descfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildDescriptor(kps.size(), kps, descfn, this->resize_rate);
    }


    if(strcmp(dvfn, "") && this->myDescriptor != NULL)
    {
        this->myDescriptor->setupImage(this->crntimg);
        this->myDescriptor->buildPatchView(kps.size(), kps, dvfn, this->resize_rate);
    }

    if(intImg != NULL)
    {
       delete this->intImg;
       this->intImg = NULL;
    }

    if(this->crntimg != NULL)
    {
       delete this->crntimg;
       this->crntimg = NULL;
    }

    return true;
}

bool HesAff::getdRds(const float block[5][5], const float sigma, float &dxds, float &dyds)
{
    float dxx, dyy, dxy, xn, xp, xy, yn, yp, xpyp, xnyn, xnyp, xpyn;
    float det = 0, v = 0, xnn, xpp, ynn, ypp, dxxx, dxxy, dxyy, dyyy;

    xy = block[2][2];
    xn = block[2][1];    xp = block[2][3];
    yn = block[1][2];    yp = block[3][2];
    xpyp = block[3][3];  xnyn = block[1][1];
    xnyp = block[3][1];  xpyn = block[1][3];

    dxx = xn + xp - 2 * xy;
    dyy = yn + yp - 2 * xy;
    dxy = (xpyp + xnyn - xnyp - xpyn) / 4;
    det = dxx * dyy - dxy * dxy;
    v   = -sigma / det;

    xnn = block[2][0];    xpp = block[2][4];
    ynn = block[0][2];    ypp = block[4][2];

    dxxx = (xpp - 2 * xp + 2 * xn - xnn) / 2;
    dyyy = (ypp - 2 * yp + 2 * yn - ynn) / 2;
    dxyy = (2 * (xn - xp) + xpyn + xpyp - xnyn - xnyp) / 2;
    dxxy = (2 * (yn - yp) + xnyp + xpyp - xnyn - xpyn) / 2;

    dxds = v * (dyy * dxxx + dyy * dxyy - dxy * dxxy - dxy * dyyy);
    dyds = v * (dxx * dyyy + dxx * dxxy - dxy * dxxx - dxy * dxyy);

    return true;
}

float HesAff::getdHds(const float block[5][5], float sigma)
{
    float t = sigma;
    float dxx, dyy, dxy, dxxxx, dyyyy, dxxyy, dxxxy, dxyyy, dhs = 0;
    float xn = block[2][1]; //window->getPixel(cx - 1, cy); // x negative 1
    float xp = block[2][3];//window->getPixel(cx + 1, cy); // x positive 1
    float xy = block[2][2];//window->getPixel(cx, cy);
    float yn = block[1][2];//window->getPixel(cx, cy - 1);
    float yp = block[3][2];//window->getPixel(cx, cy + 1);
    float xpyp = block[3][3];//window->getPixel(cx + 1, cy + 1);
    float xnyn = block[1][1];//window->getPixel(cx - 1, cy - 1);
    float xnyp = block[3][1];//window->getPixel(cx - 1, cy + 1);
    float xpyn = block[1][3];//window->getPixel(cx + 1, cy - 1);

    dxx   = xn + xp - 2 * xy;
    dyy   = yn + yp - 2 * xy;
    dxy   = (xpyp + xnyn - xnyp - xpyn) / 4;
    dxxxx = block[2][0]  + block[2][4] - 4 * (xn + xp) + 6 * xy;
    dyyyy = block[0][2]  + block[4][2] - 4 * (yn + yp) + 6 * xy;
    dxxyy = xpyp + xpyn  + xnyp + xnyn - 2 * (xn + xp + yn + yp) + 4 * xy;
    dxxxy = (block[1][0] + block[3][4] - block[3][0] - block[1][4])/4 + (xnyp + xpyn - xnyn - xpyp)/2;
    dxyyy = (block[0][1] + block[4][3] - block[0][3] - block[4][1])/4 + (xnyp + xpyn - xnyn - xpyp)/2;
    dhs   = ((dxxxx + dxxyy) * dyy + (dxxyy + dyyyy) * dxx)/2 - dxy*(dxxxy + dxyyy);
    return 2 * t * (dxx * dyy - dxy * dxy) + t * t * dhs;
}

void HesAff::affineAdapt(vector<KeyPoint *> &peaks, vector<vector<Image*> > &gaussPyramid, vector<vector<Image*> > &dxPyramid, vector<vector<Image*> > &dyPyramid)
{
    int times = 0, num = 0, nPxs = 0, cx = 0, cy = 0;
    int oct = 0, si = 0, xi = 0, yi = 0, ri = 0, x, y, wd = 0, ci = 0;
    unsigned int ioct;

    vector<KeyPoint*>::iterator it;
    vector<Cords>::iterator cit;
    vector<Cords> visited;
    vector<Board *> octFlags;
    Board *bdFlag = NULL;

    Image *dxImg = NULL, *dyImg = NULL, *gsImg = NULL;
    KeyPoint *crntKp = NULL;

    float fx = 0, fy = 0, dx = 0, dy = 0, px = 0, sx, sy = 0, dr = 0, r = 4.81;
    float itSigma = 0, ds = 0, tmpdt = 0, dt = 0, dxds, dyds, dta, theta;
    float U[2][2] = {0}, ivU[2][2] = {0}, nwU[2][2]={0}, tU[2][2];
    float M[2][2] = {0}, lm1 = 0, lm2 = 0, dlta = 0, sc = 0, tr = 0;
    float tmpU[2][2] = {0}, es[2] = {0}, block[5][5];
    bool _CNVG_ = false;

    for(ioct = 0; ioct < gaussPyramid.size(); ioct++)
    {
        gsImg  = gaussPyramid[ioct][0];
        bdFlag = new Board(gsImg->width, gsImg->height);
        octFlags.push_back(bdFlag);
    }
    ///cout<<"\t";
    for(it = peaks.begin(); it != peaks.end(); it++, num++)
    {
        crntKp  = *it;
        if(crntKp->KP == false)
        continue;

        oct     = crntKp->octIndex;
        si      = crntKp->scale;
        itSigma = crntKp->octSigma;
        U[0][0] = U[1][1] = 1;
        U[0][1] = U[1][0] = 0;
        _CNVG_  = false;
        r       = 4.5;
        times   = 0;
        bdFlag  = octFlags[oct];
        es[0]   = es[1] = 1;
        ds      = 0.08;
        dt      = 1.0;
        sx      = crntKp->sx; sy = crntKp->sy;
        dxImg   = dxPyramid[oct][si];
        dyImg   = dyPyramid[oct][si];
        gsImg   = gaussPyramid[oct][si];
        while(!_CNVG_ && times < 10)
        {
            M[0][0] = M[0][1] = M[1][0] = M[1][1] = 0;
            sc = (itSigma*3.0)/r;
            wd = (int)floor((sc*es[1]*es[0]*itSigma*3.0)/r);
            wd = wd*2+1;
            nPxs  = 0;
            visited.clear();
            for(theta = 0; theta < 6.32; theta += 0.45)
            {
                for(dr = 0; dr < 4.81; dr += 0.8)
                {
                   xi = dr*cos(theta);
                   yi = dr*sin(theta);
                   fx = sc*xi*U[0][0] + sc*yi*U[0][1];
                   fy = sc*xi*U[1][0] + sc*yi*U[1][1];
                   x  = (int)floor(sx + fx);
                   y  = (int)floor(sy + fy);
                   for(ri = -wd; ri <= wd; ri++)
                   {
                       cy = y + ri;
                       for(ci = -wd; ci <= wd; ci++)
                       {
                          cx = x + ci;
                          if(bdFlag->getPixel(cx, cy) == 0)
                          {
                             dx = dxImg->getPixel(cx, cy);
                             dy = dyImg->getPixel(cx, cy);
                             M[0][0] += dx*dx;
                             M[0][1] += dx*dy;
                             M[1][1] += dy*dy;
                             nPxs++;
                             bdFlag->setPixel(cx, cy, 1);
                             Cords pt; pt.x = cx; pt.y = cy;
                             visited.push_back(pt);
                          }
                       }
                    }
                }
            }

            for(cit = visited.begin(); cit != visited.end(); cit++)
            {
                  Cords & pt = *cit;
                  bdFlag->setPixel(pt.x, pt.y, 0);
            }
            visited.clear();
            for(yi = 0; yi < 5; yi++)
            {
               for(xi = 0; xi < 5; xi++)
               {
                    x  = xi - 2; y = yi - 2;
                    fx = sc*x*U[0][0] + sc*y*U[0][1];
                    fy = sc*x*U[1][0] + sc*y*U[1][1];
                    x  = crntKp->sx + fx;
                    y  = crntKp->sy + fy;
                    px = gsImg->getPixelBI(x, y);
                    block[yi][xi] = px;
               }
            }

            tmpdt = getdHds(block, itSigma);
            if(tmpdt*dt < 0)
            {
                ds = ds/2;
            }
            if(tmpdt > 0)
            {
               itSigma = itSigma + ds;
            }else{
               itSigma = itSigma - ds;
            }
            dt = tmpdt;
            this->getdRds(block, itSigma, dxds, dyds);
            sx = sx + ds * dxds;
            sy = sy + ds * dyds;
            ///cout<<"tm: "<<times<<"\tds: "<<ds<<"\tsig: "<<itSigma<<"\txy: "<<dxds<<"\t"<<dyds<<endl;

            M[1][0] = M[0][1];
            M[0][0] = M[0][0]/nPxs; M[0][1] = M[0][1]/nPxs;
            M[1][0] = M[1][0]/nPxs; M[1][1] = M[1][1]/nPxs;

            tmpU[0][0] = M[0][0]; tmpU[0][1] = M[0][1];
            tmpU[1][0] = M[1][0]; tmpU[1][1] = M[1][1];
            ///cout<<"M: "<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<nPxs<<endl;

            VMath::sqrtSymMat(tmpU[0][0], tmpU[0][1], tmpU[1][1], nwU, es);
            Inverse2D(U, ivU);
            tU[0][0] = nwU[0][0]*ivU[0][0] + nwU[1][0]*ivU[1][0];
            tU[0][1] = nwU[0][0]*ivU[0][1] + nwU[1][0]*ivU[1][1];
            tU[1][0] = nwU[1][0]*ivU[0][0] + nwU[1][1]*ivU[1][0];
            tU[1][1] = nwU[1][0]*ivU[0][1] + nwU[1][1]*ivU[1][1];
            tr   = tU[0][0] + tU[1][1];
            dlta = tU[0][0]*tU[1][1] - tU[1][0]*tU[0][1];
            lm1  = tr/2 + sqrt(tr*tr/4-dlta);
            lm2  = tr/2 - sqrt(tr*tr/4-dlta);
            ///cout<<"lm1: "<<lm1<<"\t"<<lm2<<endl;
            if(lm2/lm1 > 0.96)
            {
                _CNVG_ = true;
            }
            U[0][0] = nwU[0][0]; U[0][1] = nwU[0][1];
            U[1][0] = nwU[1][0]; U[1][1] = nwU[1][1];
            ///exit(0);
            times++;
        }///while(!_CNVG_)
        ///exit(0);
        if(times >= 10)
        {
           crntKp->KP = false;
        }
        ///cout<<num<<"\t"<<times<<"\t"<<lm2/lm1<<"\tdt = "<<dt<<endl;
        crntKp->a = U[0][0]; crntKp->b = U[0][1]; crntKp->c = U[1][1];
        ///cout<<crntKp->a<<"\t"<<crntKp->b<<"\t"<<crntKp->c<<endl;
        crntKp->octSigma = itSigma;
        crntKp->iscale   =  pow(2.0, crntKp->octIndex)*itSigma*HesAff::mag;
        crntKp->dscale   =  pow(2.0, crntKp->octIndex)*itSigma;
        crntKp->x = (int) round(sx * pow(2.0, crntKp->octIndex));
        crntKp->y = (int) round(sy * pow(2.0, crntKp->octIndex));
        ///cout<<"\r\r\r\r\t"<<num<<"/"<<peaks.size();
    }
    Cleaner::releaseBoards(octFlags);
    ///cout<<"done\n";
    ///cout<<endl;
}

void HesAff::affineAdapt2D(vector<KeyPoint *> &peaks, vector<vector<Image*> > &gaussPyramid)
{
    int times = 0, num = 0, nPxs = 0, cx = 0, cy = 0;
    int oct = 0, si = 0, xi = 0, yi = 0, ri = 0, x, y, wd = 0, ci = 0;
    unsigned int ioct;

    vector<KeyPoint*>::iterator it;
    vector<Cords>::iterator cit;
    vector<Cords> visited;
    vector<Board *> octFlags;
    Board *bdFlag = NULL;

    Image *dxxImg = NULL, *dyyImg = NULL, *dxyImg = NULL, *gsImg = NULL;
    KeyPoint *crntKp = NULL;

    float fx = 0, fy = 0, dxx = 0, dyy = 0, dxy = 0, px = 0, sx, sy = 0, dr = 0, r = 4.81;
    float itSigma = 0, ds = 0, tmpdt = 0, dt = 0, dxds, dyds, dta, theta;
    float U[2][2] = {0}, ivU[2][2] = {0}, nwU[2][2]={0}, tU[2][2];
    float M[2][2] = {0}, lm1 = 0, lm2 = 0, dlta = 0, sc = 0, tr = 0;
    float tmpU[2][2] = {0}, es[2] = {0}, block[5][5];
    bool _CNVG_ = false;

    for(ioct = 0; ioct < gaussPyramid.size(); ioct++)
    {
        gsImg  = gaussPyramid[ioct][0];
        bdFlag = new Board(gsImg->width, gsImg->height);
        octFlags.push_back(bdFlag);
    }
    ///cout<<"\t";
    for(it = peaks.begin(); it != peaks.end(); it++, num++)
    {
        crntKp  = *it;
        if(crntKp->KP == false)
        continue;

        oct     = crntKp->octIndex;
        si      = crntKp->scale;
        itSigma = crntKp->octSigma;
        U[0][0] = U[1][1] = 1;
        U[0][1] = U[1][0] = 0;
        _CNVG_  = false;
        r       = 4.5;
        times   = 0;
        bdFlag  = octFlags[oct];
        es[0]   = es[1] = 1;
        ds      = 0.08;
        dt      = 1.0;
        sx      = crntKp->sx; sy = crntKp->sy;
        dxxImg  = dxxPyramid[oct][si];
        dyyImg  = dyyPyramid[oct][si];
        dxyImg  = dxyPyramid[oct][si];
        gsImg   = gaussPyramid[oct][si];
        while(!_CNVG_ && times < 10)
        {
            M[0][0] = M[0][1] = M[1][0] = M[1][1] = 0;
            sc = (itSigma*3.0)/r;
            wd = (int)floor((sc*es[1]*es[0]*itSigma*3.0)/r);
            wd = wd*2+1;
            nPxs  = 0;
            visited.clear();
            for(theta = 0; theta < 6.32; theta += 0.45)
            {
                for(dr = 0; dr < 4.81; dr += 0.8)
                {
                   xi = dr*cos(theta);
                   yi = dr*sin(theta);
                   fx = sc*xi*U[0][0] + sc*yi*U[0][1];
                   fy = sc*xi*U[1][0] + sc*yi*U[1][1];
                   x  = (int)floor(sx + fx);
                   y  = (int)floor(sy + fy);
                   for(ri = -wd; ri <= wd; ri++)
                   {
                       cy = y + ri;
                       for(ci = -wd; ci <= wd; ci++)
                       {
                          cx = x + ci;
                          if(bdFlag->getPixel(cx, cy) == 0)
                          {
                             dxx = dxxImg->getPixel(cx, cy);
                             dyy = dyyImg->getPixel(cx, cy);
                             dxy = dxyImg->getPixel(cx, cy);
                             M[0][0] += dxx*dxx;
                             M[0][1] += dxy*dxy;
                             M[1][1] += dyy*dyy;
                             nPxs++;
                             bdFlag->setPixel(cx, cy, 1);
                             Cords pt; pt.x = cx; pt.y = cy;
                             visited.push_back(pt);
                          }
                       }
                    }
                }
            }

            for(cit = visited.begin(); cit != visited.end(); cit++)
            {
                  Cords & pt = *cit;
                  bdFlag->setPixel(pt.x, pt.y, 0);
            }
            visited.clear();
            for(yi = 0; yi < 5; yi++)
            {
               for(xi = 0; xi < 5; xi++)
               {
                    x  = xi - 2; y = yi - 2;
                    fx = sc*x*U[0][0] + sc*y*U[0][1];
                    fy = sc*x*U[1][0] + sc*y*U[1][1];
                    x  = crntKp->sx + fx;
                    y  = crntKp->sy + fy;
                    px = gsImg->getPixelBI(x, y);
                    block[yi][xi] = px;
               }
            }

            tmpdt = getdHds(block, itSigma);
            if(tmpdt*dt < 0)
            {
                ds = ds/2;
            }
            if(tmpdt > 0)
            {
               itSigma = itSigma + ds;
            }else{
               itSigma = itSigma - ds;
            }
            dt = tmpdt;
            this->getdRds(block, itSigma, dxds, dyds);
            sx = sx + ds * dxds;
            sy = sy + ds * dyds;
            ///cout<<"tm: "<<times<<"\tds: "<<ds<<"\tsig: "<<itSigma<<"\txy: "<<dxds<<"\t"<<dyds<<endl;

            M[1][0] = M[0][1];
            M[0][0] = M[0][0]/nPxs; M[0][1] = M[0][1]/nPxs;
            M[1][0] = M[1][0]/nPxs; M[1][1] = M[1][1]/nPxs;

            tmpU[0][0] = M[0][0]; tmpU[0][1] = M[0][1];
            tmpU[1][0] = M[1][0]; tmpU[1][1] = M[1][1];
            ///cout<<"M: "<<M[0][0]<<"\t"<<M[0][1]<<"\t"<<M[1][0]<<"\t"<<M[1][1]<<"\t"<<nPxs<<endl;

            VMath::sqrtSymMat(tmpU[0][0], tmpU[0][1], tmpU[1][1], nwU, es);
            Inverse2D(U, ivU);
            tU[0][0] = nwU[0][0]*ivU[0][0] + nwU[1][0]*ivU[1][0];
            tU[0][1] = nwU[0][0]*ivU[0][1] + nwU[1][0]*ivU[1][1];
            tU[1][0] = nwU[1][0]*ivU[0][0] + nwU[1][1]*ivU[1][0];
            tU[1][1] = nwU[1][0]*ivU[0][1] + nwU[1][1]*ivU[1][1];
            tr   = tU[0][0] + tU[1][1];
            dlta = tU[0][0]*tU[1][1] - tU[1][0]*tU[0][1];
            lm1  = tr/2 + sqrt(tr*tr/4-dlta);
            lm2  = tr/2 - sqrt(tr*tr/4-dlta);
            ///cout<<"lm1: "<<lm1<<"\t"<<lm2<<endl;
            if(lm2/lm1 > 0.96)
            {
                _CNVG_ = true;
            }
            U[0][0] = nwU[0][0]; U[0][1] = nwU[0][1];
            U[1][0] = nwU[1][0]; U[1][1] = nwU[1][1];
            ///exit(0);
            times++;
        }///while(!_CNVG_)
        ///exit(0);
        if(times >= 10)
        {
           crntKp->KP = false;
        }
        ///cout<<num<<"\t"<<times<<"\t"<<lm2/lm1<<"\tdt = "<<dt<<endl;
        crntKp->a = U[0][0]; crntKp->b = U[0][1]; crntKp->c = U[1][1];
        ///cout<<crntKp->a<<"\t"<<crntKp->b<<"\t"<<crntKp->c<<endl;
        crntKp->octSigma = itSigma;
        crntKp->iscale   =  pow(2.0, crntKp->octIndex)*itSigma*HesAff::mag;
        crntKp->dscale   =  pow(2.0, crntKp->octIndex)*itSigma;
        crntKp->x = (int) round(sx * pow(2.0, crntKp->octIndex));
        crntKp->y = (int) round(sy * pow(2.0, crntKp->octIndex));
        ///cout<<"\r\r\r\r\t"<<num<<"/"<<peaks.size();
    }
    Cleaner::releaseBoards(octFlags);
    ///cout<<"done\n";
    ///cout<<endl;
}


vector<KeyPoint*> HesAff::FindOrientByGrad(vector<KeyPoint *> &kps, vector<vector<Image *> > & GaussianOctaves)
{
    vector<vector<float> > gmat;
    vector<KeyPoint*> newkps;
    KeyPoint *newkp;

    float sigma, m, theta, degree, fx, fy;
    float U[2][2];
    int c, index;
    int indexa, indexb, indexc, j;
    float thetaa, thetab, thetac;
    unsigned int i;
    float maxval,maxp;

    bool valid;
    float thetas[NOrient] = {0};
    float weight;

    for (i = 0; i < kps.size(); i++)
    {
        if(kps[i]->KP == false)
        continue;

        sigma = 1.5 * pow(2.0, (kps[i]->fscale)/(float) SCALES) *HesAff::_SIGMA;
        Filter::GaussianKernel2D(sigma, gmat);
        c = gmat.size()/2;

        for(j = 0; j < NOrient; j++)
        {
            thetas[j] = 0;
        }
        U[0][0] = kps[i]->a; U[0][1] = kps[i]->b;
        U[1][0] = kps[i]->b; U[1][1] = kps[i]->c;

        for (int y = -c; y <= c; y++)
        {
            for (int x = -c; x <= c; x++)
            {
                if (sqrt((float) x*x + y*y) > 3.0 * sigma)
                    continue;

                fx = x*U[0][0] + y*U[0][1];
                fy = x*U[1][0] - y*U[1][1];

                valid = Filter::GetPixOrientation((int) (kps[i]->sx + fx + 0.5), (int) (kps[i]->sy + fy + 0.5),
                                                  GaussianOctaves[kps[i]->octave][kps[i]->scale], m, theta);

                if(valid)
                {
                    degree = theta / PI * 180.0 + 180.0;
                    index  = ((int) (degree / DEGREE));
                    weight = m*gmat[y + c][x + c];
                    index  = index%NOrient;
                    thetas[index] +=  weight;

                }
            }
        }

        vector<vector<float> >::iterator it;
        vector<float> crntvect;
        for(it = gmat.begin(); it != gmat.end(); it++)
        {
            crntvect = *it;
            crntvect.clear();
        }
        gmat.erase(gmat.begin(),gmat.end());


        for (int j = 0; j < 6; j++)
            VMath::SmoothHistogram(thetas, NOrient);

        maxval = VMath::maxVec(thetas, NOrient);

        for (int j = 0; j < NOrient; j++)
        {
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

            if(thetas[j] == maxval)
            {
                kps[i]->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;
            }
            else
            {
                newkp = new KeyPoint();
                memcpy(newkp, kps[i], sizeof(KeyPoint));
                if(/**/this->mydesc == ERIFT || this->mydesc == NERIFT || /**/
                        this->mydesc == NSPIN || this->mydesc == SPIN)
                {
                    newkp->KP = false;
                }
                newkp->ori = ((float) j + maxp + 0.5) * 2.0 * PI / (float) NOrient - PI;

                newkps.push_back(newkp);

            }
        }
    }
    return  newkps ;
}

void HesAff::test()
{
    const char *srcFn = "/home/wlzhao/datasets/vgg/graf/graf1.jpg";
    const char *dstFn = "/home/wlzhao/datasets/vgg/graf1.keys";
    const char *drdir = "/home/wlzhao/datasets/vgg/";
    HesAff *myhes = new HesAff();
    myhes->KeypointBuild(srcFn, dstFn, "", "");

}

HesAff::~HesAff()
{
    //dtor
}
