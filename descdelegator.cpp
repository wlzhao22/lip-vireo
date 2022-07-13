#include "descdelegator.h"

#include "index_template.h"
#include "cleaner.h"
#include "filter.h"
#include "vmath.h"

#include <cmath>

///for SIFT descriptor
const int DescDelegator::GRID        = 4;
const int DescDelegator::NumOrient   = 8;
const int DescDelegator::DSize       = GRID*GRID*NumOrient;
const int DescDelegator::DEGREE      = (360/NumOrient);
const int DescDelegator::PSIFTLen    = 36;
const int DescDelegator::PatchMag    = 20;


///for SPIN descriptor
const int DescDelegator::Ints      = 10;
const int DescDelegator::Dist      = 10;
const float DescDelegator::alpha   = 1.8f;
const float DescDelegator::belta   = 3.6f;
const int DescDelegator::CLR_DEPTH = 255;

///for RIFT descriptor
const int     DescDelegator::Ori     = 8;
const float   DescDelegator::Theta_per_bin = PI2/DescDelegator::Ori;
const int     DescDelegator::DistF   = 4;

DescDelegator::DescDelegator(DESC desc0, const int PatchSize0)
{
    switch(desc0)
    {
    case SPIN:
    {
        cout<<"Descriptor .............................. un-normalized SPIN\n";
        this->extractFeat = &DescDelegator::getSPINDescriptor;
        this->featLen = DescDelegator::Ints*DescDelegator::Dist;
        break;
    }
    case NSPIN:
    {
        cout<<"Descriptor .............................. normalized SPIN\n";
        this->extractFeat = &DescDelegator::getSPINDescriptor;
        this->featLen = DescDelegator::Ints*DescDelegator::Dist;

        break;
    }
    case RIFT:
    {
        cout<<"Descriptor .............................. un-normalized RIFT\n";
        this->extractFeat = &DescDelegator::getRIFTDescriptor;
        this->featLen = DescDelegator::DistF*DescDelegator::Ori;
        break;
    }
    case NRIFT:
    {
        cout<<"Descriptor .............................. normalized RIFT\n";
        this->extractFeat = &DescDelegator::getRIFTDescriptor;
        this->featLen = DescDelegator::DistF*DescDelegator::Ori;
        break;
    }
    case CVSIFT:
    {
        cout<<"Descriptor .............................. un-normalized SIFT\n";
        this->extractFeat = &DescDelegator::getSIFTDescriptor;
        this->featLen = DSize;
        break;
    }
    case NCVSIFT:
    {
        cout<<"Descriptor .............................. normalized SIFT\n";
        this->extractFeat = &DescDelegator::getSIFTDescriptor;
        this->featLen = DSize;

        break;
    }
    case CTSIFT:
    {
        cout<<"Descriptor .............................. um-normalized SIFT\n";
        this->extractFeat = &DescDelegator::getSIFTDescriptor;
        this->featLen = DSize;
        break;
    }
    case NCTSIFT:
    {
        cout<<"Descriptor .............................. normalized SIFT\n";
        this->extractFeat = &DescDelegator::getSIFTDescriptor;
        this->featLen = DSize;

        break;
    }
    case GBLUR:
    {
        cout<<"Descriptor .............................. un-normalized GBLUR\n";
        this->extractFeat = &DescDelegator::getGBLURDescriptor;
        this->featLen = DSize;
        break;
    }
    case NGBLUR:
    {
        cout<<"Descriptor .............................. normalized GBLUR\n";
        this->extractFeat = &DescDelegator::getGBLURDescriptor;
        this->featLen = DSize;

        break;
    }
    default:
    {
        cout<<"No descriptor has been chosen!\n";
    }
    }

    gmat = DescDelegator::GaussianWeight2D(PatchSize0);

    this->descoption = desc0;
    this->PatchSize  = PatchSize0;
}

int DescDelegator::getFeatDim()
{
    return this->featLen;
}

vector<vector<float> > DescDelegator::GaussianWeight2D(const int winSize)
{
    int   c     = winSize/2, i;
    int   dim   = c*2 + 1;
    float sigma = c/2;
    vector<vector<float> > mat;
    for (i = 0; i < dim; i++)
    {
        vector<float> row(dim);
        mat.push_back(row);
    }

    float s2 = sigma * sigma;
    float v;
    int   x, y;

    for (y = -c; y <= c; y++)
    {
        for (x = -c; x <= c; x++)
        {
            v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
            mat[c+y][c+x] = v;
        }
    }

    return mat;
}

int DescDelegator::getDescriptor(float *feats, const float *myWin)
{
    (*this.*extractFeat)(feats, myWin);
    return 1;
}

int DescDelegator::getSIFTDescriptor(float *feats, const float *myWin)
{
    int i = 0, j = 0;
    assert(feats);
    memset(feats, 0, sizeof(float)*DSize);

    float angle = 0, mag = 0;
    float histWidth = PatchSize/(GRID + 0.0f);
    float Theta_per_bin = NumOrient/PI2;
    float rbin, cbin, obin;
    float rdr, rdc, rda, a, b;
    int rbin0, cbin0, obin0, ir, ic, ia;
    int rloc, cloc, oloc, row, bloc;
    int crnt, prev, next;
    float dx, dy, da;

    for(i = 1; i < PatchSize-1; i++)
    {
        rbin = i/histWidth;
        rbin = rbin - 0.5;
        rbin0 = (int)floor(rbin);
        dy = rbin - rbin0;

        for(j = 1; j < PatchSize-1; j++)
        {
            cbin = j/histWidth;
            cbin = cbin - 0.5;
            cbin0 = (int)floor(cbin);
            dx = cbin - cbin0;

            crnt = i*PatchSize;
            a = myWin[crnt + j + 1] - myWin[crnt + j - 1]; //dx
            next = crnt + PatchSize;
            prev = crnt - PatchSize;
            b = myWin[next + j] - myWin[prev +j];//dy
            mag = sqrt(a*a +  b*b);
            angle = atan2(b, a);
            angle = angle < 0? (angle + PI2):angle;

            obin  = angle*Theta_per_bin;
            obin0 = (int) floor(obin);
            da = obin - obin0;

            mag = mag*gmat[i][j];

            for(ir = 0; ir <= 1; ir++)
            {
                rloc = rbin0 + ir;

                if(rloc >= GRID || rloc < 0)
                    continue;

                rdr = mag*((ir == 0)?(1 - dy):dy);
                row = rloc*GRID;

                for(ic = 0; ic <= 1; ic++)
                {
                    cloc = cbin0 + ic;

                    if(cloc >= GRID || cloc < 0)
                        continue;

                    rdc  = rdr*((ic == 0)?(1 - dx):dx);
                    bloc = (row + cloc)*NumOrient;

                    for(ia = 0; ia <= 1; ia++)
                    {
                        oloc = (obin0 + ia)%NumOrient;
                        rda = rdc*((ia == 0)?(1 - da):da);
                        feats[bloc+oloc] = feats[bloc+oloc] + rda;
                    }
                }
            }//end-for, end of interpolation
        }//end inner for
    }//end outer for

    ///bool GOOD = VMath::sqrtSIFTNorm(this->featsBin, DSize);
    bool GOOD = VMath::SIFTNorm(feats, DSize);

    for(i = 0; i < DSize; i++)
    {
        feats[i] = round(feats[i]);
    }

    if(GOOD)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}


int DescDelegator::getSPINDescriptor(float *feats, const float *myWin)
{
    int i, j, k, m;
    for(i = 0; i < this->featLen; i++)
    {
        feats[i] = 0.0f;
    }

    float alpha2, belta2;
    int pos1, pos2, radius = PatchSize/2;
    float  intense, _intense, dist, _dist, step_dist, step_ints;
    unsigned int  _idx;

    step_dist = PatchSize/DescDelegator::Dist;
    step_ints = DescDelegator::CLR_DEPTH/DescDelegator::Ints;
    alpha2 = -2.0f*DescDelegator::alpha*DescDelegator::alpha;
    belta2 = -2.0f*DescDelegator::belta*DescDelegator::belta;

    for(i = 0; i < PatchSize; i++)
    {
        pos1 = i*PatchSize;
        for(j = 0; j< PatchSize; j++)
        {
            pos2 = pos1 + j;
            intense = myWin[pos2];
            dist = (i - radius)*(i - radius)  + (j - radius)*(j - radius);
            dist = sqrt(dist);

            if(dist > radius)
            continue;

            for(k = 0; k < DescDelegator::Ints; k++)
            {
                _intense = intense - k*step_ints;
                _intense = (_intense*_intense)/belta2;
                for(m = 0; m < DescDelegator::Dist; m++)
                {
                    _dist = dist - m*step_dist;
                    _dist = (_dist*_dist)/alpha2;
                    _idx = k*DescDelegator::Ints + m;
                    feats[_idx] +=  10*exp(_intense+_dist);
                }
            }
        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        feats[i] = round(feats[i]);
    }

    return 0;
}


int DescDelegator::getRIFTDescriptor(float *feats, const float *myWin)
{
    int i, j;
    float a, b, mag;
    float dx, dy, angle, angle1, angle2, a_r, dist;
    int pos, row, next_y, prev_y, radius = PatchMag;
    unsigned int idx, idx_d, idx_a;

    memset(feats, 0, sizeof(float)*this->featLen);

    for(i = 1; i < PatchSize-1; i++)
    {
        row = i*PatchSize;
        for(j = 1; j < PatchSize-1; j++)
        {
            pos    = row + j;
            next_y = row + PatchSize;
            prev_y = row - PatchSize;
            a      = myWin[ pos + 1]   - myWin[pos - 1]; ///dx
            b      = myWin[next_y + j] - myWin[prev_y + j];///dy

            mag    = sqrt(a*a +  b*b);
            angle1 = atan2(b, a);

            dx     = j - radius;
            dy     = i - radius;
            idx_d  = idx_table[pos];
            dist   = sqrt(dx*dx+dy*dy);
            if(dist <= 5)
            {
                idx_d = 0;
            }else if(dist <= 10)
            {
                idx_d = 1;
            }else if(dist < 15)
            {
                idx_d = 2;
            }else if(dist <= 20)
            {
                idx_d = 3;
            }else{
                continue;
            }

            angle2 = atan2(dy, dx);
            angle = angle1 - angle2;

            while(angle < 0.0)
                angle += PI2;

            while(angle >= PI2)
                angle -= PI2;

            angle = angle/DescDelegator::Theta_per_bin;
            idx_a = (int)floor(angle);
            a_r   = angle - idx_a;
            idx_a = idx_a < 0?(DescDelegator::Ori - 1):idx_a;

            idx   = idx_d*DescDelegator::Ori + idx_a;

            feats[idx] += 10*mag*(1 - a_r);
            idx = idx_d*DescDelegator::Ori + (idx_a + 1)%DescDelegator::Ori;
            feats[idx] += 10*mag*a_r;

        }
    }

    for(i = 0; i < this->featLen; i++)
    {
        feats[i] = round(feats[i]);
    }
    return 1;
}


int DescDelegator::getGBLURDescriptor(float *feats, const float *myWin)
{
    ///TO-DO by Ali
    return 1;
}

DescDelegator::~DescDelegator()
{
    Cleaner::clear2DArray(gmat);
}
