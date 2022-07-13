#include "kpdrawer.h"
#include "vmath.h"

#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

const unsigned int KPDrawer::_cwith    = 14;
const unsigned int KPDrawer::_cheight  = 14;
const unsigned int KPDrawer::_awidth   = 10; //arrow width
const unsigned int KPDrawer::_aheight  = 5; //arrow height
const float KPDrawer::PI0 = 3.1415926f;

/**************************************
const float KPDrawer::WHITE[0] = 255;
const float KPDrawer::WHITE[1] = 255;
const float KPDrawer::WHITE[2] = 255;
const float KPDrawer::RED[0] = 255;
const float KPDrawer::RED[1] = 0;
const float KPDrawer::RED[2] = 0;
**************************************/

KPDrawer::KPDrawer(const bool CIRCLE, const Detector det_option)
{
    WHITE[0] = 255;
    WHITE[1] = 255;
    WHITE[2] = 255;
    RED[0]   = 0;
    RED[1]   = 0;
    RED[2]   = 255;
    PURPL[0] = 128;
    PURPL[1] = 0;
    PURPL[2] = 128;
    YELLW[0] = 0;
    YELLW[1] = 255;
    YELLW[2] = 255;

    this->draw_circle = CIRCLE;
    this->myoption    = det_option;
}

void KPDrawer::draw_shapes(vector <KeyPoint *> kps, const char *srcimgfn, const float scale_rate0, const char *dstimgfn)
{
    CImage *img = new CImage(srcimgfn);
    if(img == NULL)
        return ;

    vector<KeyPoint*>::iterator vit;
    KeyPoint* crnt_pt;
    unsigned int tmpx, tmpy, tx, ty, i = 0;
    int x0, y0;
    float si, co;

    for(vit = kps.begin(); vit != kps.end(); vit++)
    {
        crnt_pt = *vit;
        if(!crnt_pt->KP)
            continue;

        tmpx = (unsigned int)(crnt_pt->x * scale_rate0);
        tmpy = (unsigned int)(crnt_pt->y * scale_rate0);
        draw_cross(img, tmpx, tmpy);
        //draw_ellipse(img, tmpx, tmpy, 1, 0, 1, 3);
        i++;
        /**/
        si   = sin(crnt_pt->ori);
        co   = cos(crnt_pt->ori);
        x0   = crnt_pt->iscale;
        y0 = 0;
        tx   = round(co*x0 - si*y0) + tmpx;
        ty   = round(si*x0 + co*y0) + tmpy;
        /**/

        if(this->draw_circle)
        {
            draw_ellipse(img, tmpx, tmpy, crnt_pt->a, crnt_pt->b, crnt_pt->c, crnt_pt->iscale+1);
            draw_ellipse(img, tmpx, tmpy, crnt_pt->a, crnt_pt->b, crnt_pt->c, crnt_pt->iscale);
            ///draw_ellipse(img, tmpx, tmpy, crnt_pt->a, crnt_pt->b, crnt_pt->c, crnt_pt->iscale-1);
            if(crnt_pt->a == 1 && crnt_pt->c == 1)
            {
                lineto(img, tmpx, tmpy, tx, ty);
                //lineto(img, tmpx-1, tmpy-1, tx-1, ty-1);
                draw_arrow(img, crnt_pt->ori, tx, ty);
            }
        }
    }
    img->save(dstimgfn);
    delete  img;
}

void KPDrawer::draw_ellipse(CImage *img, const unsigned int x0, const unsigned int y0, const float a,
                            const float b, const float c, const float sc)
{

    float theta = 0, x = 0, y = 0, xbar = 0, ybar = 0;
    float eigns[2], eignv[2][2];
    unsigned int xi, yi, i;

    VMath::eigvtSymMat(a, b, c, eignv, eigns);
    ///cout<<a<<"\t"<<b<<"\t"<<c<<endl;
    theta =  atan2(eignv[1][0], eignv[0][0]);
    ///cout<<theta<<"\t"<<eigns[0]<<"\t"<<eigns[1]<<endl;
    //theta = 0;

    unsigned int div  = 400*2*KPDrawer::PI0*sqrt(0.5*(eigns[0]*eigns[0]+eigns[1]*eigns[1]));
    float delta_alpha = KPDrawer::PI0/div, alpha = 0;

    div = 2*div;
    for(i = 0; i < div; i++)
    {
        x = sc*eigns[0]*cos(alpha)*0.33;
        y = sc*eigns[1]*sin(alpha)*0.33;

        xbar = x*cos(theta) - y*sin(theta);
        ybar = y*cos(theta) + x*sin(theta);

        xi = x0 + (unsigned int)round(xbar);
        yi = y0 + (unsigned int)round(ybar);
        img->setPixel(xi-1, yi, YELLW);
        img->setPixel(xi, yi, YELLW);
        img->setPixel(xi, yi-1, YELLW);
        alpha += delta_alpha;
    }
    return ;
}

void KPDrawer::draw_rects(vector <KeyPoint *> kps, const char *srcimgfn,
                          const float scale_rate0, const char *dstimgfn)
{
    int i, j, rd;
    int x_cord[5], y_cord[5];
    float rectx[5], recty[5], si, co;

    CImage *img = new CImage(srcimgfn);
    if(img == NULL)
        return ;

    vector<KeyPoint*>::iterator vit;
    KeyPoint* crnt_pt;
    int tmpx, tmpy;

    for(vit = kps.begin(); (vit != kps.end()); vit++)
    {
        crnt_pt = *vit;
        if(!crnt_pt->KP)
            continue;

        tmpx = (int)(crnt_pt->x * scale_rate0);
        tmpy = (int)(crnt_pt->y * scale_rate0);
        draw_cross(img, tmpx, tmpy);

        si       = sin(crnt_pt->ori);
        co       = cos(crnt_pt->ori);
        rd       = round(crnt_pt->iscale/2);
        rectx[0] = recty[0] = -1*rd;
        rectx[1] =  rd;
        recty[1] = -1*rd;
        rectx[2] = recty[2] = rd;
        rectx[3] = -1*rd;
        recty[3] =  rd;
        rectx[4] =  rd;
        recty[4] =  0;
        for(i = 0; i < 5; i++)
        {
            x_cord[i] = round(co*rectx[i] - si*recty[i]) + tmpx;
            y_cord[i] = round(si*rectx[i] + co*recty[i]) + tmpy;
        }

        if(!this->draw_circle)
        {
            continue;
        }
        for(i = 0; i < 4; i++)
        {
            j = (i+1)%4;
            this->lineto(img, x_cord[i], y_cord[i], x_cord[j], y_cord[j]);
        }
        this->lineto(img, tmpx, tmpy, x_cord[4], y_cord[4]);
    }
    img->save(dstimgfn);
    delete img;
    return ;
}

void KPDrawer::lineto(CImage *Img, const int x0, const int y0, const int x1, const int y1)
{
    float sy, sx, ey, ex;
    int i, height, width;
    float dx, dy, absval, ratio, px, py;
    int crnt_y, crnt_x;

    height = Img->height;
    width = Img->width;
    sx = x0;
    sy = y0;
    ex = x1;
    ey = y1;

    dx = ex - sx;
    dy = ey - sy;

    if(dx == 0)
    {
        absval = VMath::absx(dy);
        crnt_x = ex;
        if(crnt_x < 0)
            return ;
        for(i = 1; i < absval; i++)
        {
            crnt_y = ey - VMath::Sign(dy)*i;
            if(crnt_y < 0)
                continue;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x,   crnt_y, YELLW);
                Img->setPixel(crnt_x+1, crnt_y, YELLW);
            }
        }
        return ;
    }

    if(dy == 0)
    {
        absval = VMath::absx(dx);
        crnt_y = ey;
        if(crnt_y < 0)
            return ;
        for(i = 1; i < absval; i++)
        {
            crnt_x = ex - VMath::Sign(dx)*i;
            if(crnt_x < 0)
                continue ;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, YELLW);
                Img->setPixel(crnt_x+1, crnt_y, YELLW);
            }
        }
        return;
    }

    ratio = dy/dx;
    py = 0;
    px = 0;
    ratio = VMath::absx(ratio);

    if(ratio > 1)
    {
        absval = dy >= 0?dy:(-dy);
        for(i = 1; i < absval; i++)
        {
            px = i/ratio;
            crnt_x = sx + VMath::Sign(dx)*(int)round(px);
            crnt_y = sy + VMath::Sign(dy)*i;
            if(crnt_x < 0 || crnt_y < 0)
                continue ;
            if((crnt_x < width) && (crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, YELLW);
                Img->setPixel(crnt_x+1, crnt_y, YELLW);
            }
        }
    }
    else  //step by x
    {
        absval = dx>=0?dx:(-dx);
        for(i = 1; i < absval; i++)
        {
            py = i*ratio;
            crnt_y = sy  + VMath::Sign(dy)*(int)round(py);
            crnt_x = sx  + VMath::Sign(dx)*i;
            if(crnt_x < 0 || crnt_y < 0)
                continue ;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y,   YELLW);
                Img->setPixel(crnt_x+1, crnt_y, YELLW);
            }
        }
    }
    return ;
}


void KPDrawer::lineto(CImage *Img, const int x0, const int y0,
                      const int x1, const int y1, const float *color)
{
    float sy, sx, ey, ex;
    int i, height, width;
    float dx, dy, absval, ratio, px, py;
    int crnt_y, crnt_x;

    height = Img->height;
    width = Img->width;
    sx = x0;
    sy = y0;
    ex = x1;
    ey = y1;

    dx = ex - sx;
    dy = ey - sy;

    if(dx == 0)
    {
        absval = VMath::absx(dy);
        crnt_x = ex;
        if(crnt_x < 0)
            return ;
        for(i = 1; i < absval; i++)
        {
            crnt_y = ey - VMath::Sign(dy)*i;
            if(crnt_y < 0)
                continue;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x,   crnt_y, color);
                Img->setPixel(crnt_x+1, crnt_y, color);
            }
        }
        return ;
    }

    if(dy == 0)
    {
        absval = VMath::absx(dx);
        crnt_y = ey;
        if(crnt_y < 0)
            return ;
        for(i = 1; i < absval; i++)
        {
            crnt_x = ex - VMath::Sign(dx)*i;
            if(crnt_x < 0)
                continue ;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, color);
                Img->setPixel(crnt_x+1, crnt_y, color);
            }
        }
        return;
    }

    ratio = dy/dx;
    py = 0;
    px = 0;
    ratio = VMath::absx(ratio);

    if(ratio > 1)
    {
        absval = dy >= 0?dy:(-dy);
        for(i = 1; i < absval; i++)
        {
            px = i/ratio;
            crnt_x = sx + VMath::Sign(dx)*(int)round(px);
            crnt_y = sy + VMath::Sign(dy)*i;
            if(crnt_x < 0 || crnt_y < 0)
                continue ;
            if((crnt_x < width) && (crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, color);
                Img->setPixel(crnt_x+1, crnt_y, color);
            }
        }
    }
    else  //step by x
    {
        absval = dx>=0?dx:(-dx);
        for(i = 1; i < absval; i++)
        {
            py = i*ratio;
            crnt_y = sy  + VMath::Sign(dy)*(int)round(py);
            crnt_x = sx  + VMath::Sign(dx)*i;
            if(crnt_x < 0 || crnt_y < 0)
                continue ;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y,   color);
                Img->setPixel(crnt_x+1, crnt_y, color);
            }
        }
    }
    return ;
}

void KPDrawer::draw_arrow(CImage *img, const float theta0,
                          const unsigned int x0, const unsigned int y0)
{
    float si  = sin(theta0);
    float co  = cos(theta0);
    float R[2][2], iR[2][2];
    R[0][0]  = co;
    R[0][1]  = si;
    R[1][0]  = -1*si;
    R[1][1]  = co;
    iR[0][0] = co;
    iR[0][1] = -1*si;
    iR[1][0] = si;
    iR[1][1] = co;
    float upcrd[2],  lwcrd[2];
    int   iupcrd[2], ilwcrd[2];
    int width  = img->width;
    int height = img->height;

    //cout<<"bug-1:\t"<<R[0][0]<<"\t"<<R[0][1]<<"\t"<<R[1][0]<<"\t"<<R[1][1]<<endl;
    upcrd[0] = x0*R[0][0] + y0*R[0][1] - KPDrawer::_awidth;
    upcrd[1] = x0*R[1][0] + y0*R[1][1] + KPDrawer::_aheight;

    lwcrd[0] = x0*R[0][0] + y0*R[0][1] - KPDrawer::_awidth;
    lwcrd[1] = x0*R[1][0] + y0*R[1][1] - KPDrawer::_aheight;

    ///cout<<"bug0:\t"<<x0<<"\t"<<y0<<"\t"<<iupcrd[0]<<"\t"<<iupcrd[1]<<endl;

    iupcrd[0] = (int)round(upcrd[0]*iR[0][0] + upcrd[1]*iR[0][1]);
    iupcrd[1] = (int)round(upcrd[0]*iR[1][0] + upcrd[1]*iR[1][1]);

    ///cout<<"bug1:\t"<<iupcrd[0]<<"\t"<<iupcrd[1]<<endl;

    iupcrd[0] = iupcrd[0]<0?0:iupcrd[0];
    iupcrd[0] = iupcrd[0]>=width?(width-1):iupcrd[0];
    iupcrd[1] = iupcrd[1]<0?0:iupcrd[1];
    iupcrd[1] = iupcrd[1]>=height?(height-1):iupcrd[1];

    ///cout<<"bug2:\t"<<iupcrd[0]<<"\t"<<iupcrd[1]<<endl;

    ilwcrd[0] = (int)round(lwcrd[0]*iR[0][0] + lwcrd[1]*iR[0][1]);
    ilwcrd[1] = (int)round(lwcrd[0]*iR[1][0] + lwcrd[1]*iR[1][1]);

    ilwcrd[0] = ilwcrd[0]<0?0:ilwcrd[0];
    ilwcrd[0] = ilwcrd[0]>width?(width-1):ilwcrd[0];
    ilwcrd[1] = ilwcrd[1]<0?0:ilwcrd[1];
    ilwcrd[1] = ilwcrd[1]>=height?(height-1):ilwcrd[1];
    ///cout<<"bug3:\t"<<iupcrd[0]<<"\t"<<iupcrd[1]<<endl;
    KPDrawer::lineto(img, iupcrd[0], iupcrd[1], x0, y0, RED);
    KPDrawer::lineto(img, ilwcrd[0], ilwcrd[1], x0, y0, RED);

    return ;
}

void KPDrawer::draw_cross(CImage *img, const unsigned int x0, const unsigned int y0)
{
    int dx, dy, x, y;

    for(dx = -5; dx <= 5; dx++)
    {
        x = x0 + dx;
        img->setPixel(x, y0, RED);
    }
    for(dy = -5; dy <= 5; dy++)
    {
        y = y0 + dy;
        img->setPixel(x0, y, RED);
    }
    return ;
}

void KPDrawer::load_siftgeo(vector<KeyPoint*> &kps, const char *srcfn)
{
    ifstream inStrm;
    inStrm.open(srcfn, ios::in|ios::binary);
    unsigned int  d0, nline = 0;
    unsigned char *feat  = NULL;
    KeyPoint    *crnt_pt = NULL;

    if(!inStrm.is_open())
    {
        cout<<"File '"<<srcfn<<"' cannot open for read!\n";
        return ;
    }
    unsigned long bg = inStrm.tellg();
    inStrm.seekg(0, ios::end);
    unsigned long sz = inStrm.tellg();
    sz = sz - bg;
    unsigned int dim = 128, row = 0;

    feat   = new unsigned char[dim];
    float *properties = new float[9];
    row    = sz/(sizeof(float)*9+sizeof(int)+sizeof(unsigned char)*dim);
    inStrm.close();

    inStrm.open(srcfn, ios::in|ios::binary);
    while(!inStrm.eof() && nline < row)
    {
        inStrm.read((char*)properties, 9*sizeof(float));
        inStrm.read((char*)&d0, sizeof(int));
        assert(dim == d0);

        inStrm.read((char*)feat, dim*sizeof(unsigned char));

        crnt_pt    = new KeyPoint();
        crnt_pt->x = properties[0];
        crnt_pt->y = properties[1];

        crnt_pt->iscale   = 2*properties[2];
        crnt_pt->ori      = properties[3];
        crnt_pt->funcVal  = properties[8];
        kps.push_back(crnt_pt);
        nline++;
    }
    inStrm.close();
    assert(nline == row);

    delete [] feat;
    delete [] properties;
    return ;
}

void KPDrawer::load_keys(vector<KeyPoint*> &kps, const char *srcfn)
{
    KeyPoint *crnt_pt;
    ifstream *inStrm = new ifstream(srcfn, ios::in);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<srcfn<<"'cannot open for read!\n"<<endl;
        exit(0);
    }

    unsigned int num = 0, dim = 0;
    float fx = 0, fy = 0;
    (*inStrm)>>dim;
    (*inStrm)>>num;
    cout<<num<<"\t"<<dim<<endl;
    int counter = 0;
    while(!inStrm->eof() && counter < 60)
    {
        crnt_pt = new KeyPoint();
        (*inStrm)>>fx;
        (*inStrm)>>fy;
        //cout<<fx<<"\t"<<fy<<endl;
        crnt_pt->x = (int)round(fx);
        crnt_pt->y = (int)round(fy);
        (*inStrm)>>crnt_pt->a;
        (*inStrm)>>crnt_pt->b;
        (*inStrm)>>crnt_pt->c;
        /**
        (*inStrm)>>crnt_pt->iscale;
        (*inStrm)>>crnt_pt->ori;
        (*inStrm)>>crnt_pt->funcVal;
        **/
        crnt_pt->KP = true;
        crnt_pt->iscale = 900.0f;
        counter++;

        kps.push_back(crnt_pt);
    }

    inStrm->close();
}

void KPDrawer::clear_kps(vector<KeyPoint*> &kps)
{
    vector<KeyPoint*>::iterator vit;
    KeyPoint* crnt_pt;

    for(vit = kps.begin(); vit != kps.end(); vit++)
    {
        crnt_pt = *vit;
        delete crnt_pt;
    }
    kps.clear();
}

void KPDrawer::test()
{
    const char *srcimgfn = "/home/wlzhao/datasets/vgg/graf/graf3.jpg";
    const char *srcfn    = "/home/wlzhao/datasets/vgg/graf3.keys";
    const char *dstimgfn = "/home/wlzhao/datasets/vgg/graf3_plot.jpg";

    KPDrawer *drawer = new KPDrawer(true, hessian);
    vector<KeyPoint*> kps;
    KPDrawer::load_keys(kps, srcfn);
    ///KPDrawer::load_siftgeo(kps, srcfn);
    cout<<kps.size()<<endl;
    drawer->draw_shapes(kps, srcimgfn, 1.0, dstimgfn);
    KPDrawer::clear_kps(kps);

}
