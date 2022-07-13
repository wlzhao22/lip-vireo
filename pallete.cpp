#include "pallete.h"
#include "vmath.h"

#include <iostream>
#include <cmath>
#include <map>


const unsigned int Pallete::WHITE = 0x00FFFFFF; //255, 255, 255
const unsigned int Pallete::RED   = 0x00FF0000; //0,0, 255
const unsigned int Pallete::YELLW = 0x00FFFF00; //0, 255, 255
const unsigned int Pallete::PURPL = 0x00800080; //128, 0, 128
const unsigned int Pallete::BLACK = 0x00000000; //0, 0, 0

const float Pallete::PI0  = 3.14159265358979323846;
const float Pallete::PI3  = 6.28318531;

void Pallete::InitView(CImage *img, VIREO_COLOR color)
{
    int i = 0, j = 0;
    float bgcolor[3];

    Pallete::fetchColor(color, bgcolor);
    for(i = 0; i < img->height; i++)
    {
        for(j = 0; j < img->width; j++)
        {
            img->setPixel(j, i, bgcolor);
        }
    }
}

bool Pallete::fetchColor(VIREO_COLOR color, float mycolor[3])
{
    unsigned char *pt;
    switch(color)
    {
    case _vireo_white:
    {
        pt = (unsigned char*)&Pallete::WHITE;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
        break;
    }
    case _vireo_red:
    {
        pt = (unsigned char*)&Pallete::RED;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
        break;
    }
    case _vireo_purple:
    {
        pt = (unsigned char*)&Pallete::PURPL;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
        break;
    }
    case _vireo_black:
    {
        pt = (unsigned char*)&Pallete::BLACK;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
        break;
    }
    case _vireo_yellow:
    {
        pt = (unsigned char*)&Pallete::YELLW;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
        break;
    }
    default:
    {
        pt = (unsigned char*)&Pallete::WHITE;
        mycolor[0] = (float)(*pt);
        pt++;
        mycolor[1] = (float)(*pt);
        pt++;
        mycolor[2] = (float)(*pt);
    }
    }
    return true;
}

void Pallete::display(CImage *img, const unsigned int x0, const unsigned int y0, CornType mytype)
{
    int dx, dy, x, y;
    float cred[3];
    float cblk[3];
    unsigned char *pt;
    pt = (unsigned char*)&Pallete::RED;
    cred[0] = (float)(*pt);
    pt++;
    cred[1] = (float)(*pt);
    pt++;
    cred[2] = (float)(*pt);
    pt = (unsigned char*)&Pallete::BLACK;
    cblk[0] = (float)(*pt);
    pt++;
    cblk[1] = (float)(*pt);
    pt++;
    cblk[2] = (float)(*pt);

    if(mytype == _JUNC)
    {
        for(dx = -5; dx <= 5; dx++)
        {
            x = x0 + dx;
            img->setPixel(x, y0, cred);
        }
        for(dy = -5; dy <= 5; dy++)
        {
            y = y0 + dy;
            img->setPixel(x0, y, cred);
        }
    }
    else if (mytype == _CORNER)
    {
        /**/
        for(dx = -5; dx <= 5; dx++)
        {
            x = x0 + dx;
            img->setPixel(x, y0, cred);
        }
        for(dy = -5; dy <= 5; dy++)
        {
            y = y0 + dy;
            img->setPixel(x0, y, cred);
        }
        /**/
        img->setPixel(x0, y0, cred);
    }
    else
    {
        //cout<<"hello\n";
        img->setPixel(x0, y0, cblk);
    }
}

void Pallete::drawEllipse(CImage *img, KeyPoint *crnt_pt, VIREO_COLOR fcolor)
{
    float angle = 0, x, y, xbar, ybar;
    float eigns[2];
    unsigned int xi, yi, i;
    float dx, dy;
    float cred[3];
    float cblk[3];
    Pallete::fetchColor(_vireo_red, cred);
    Pallete::fetchColor(fcolor, cblk);

    eigns[0] = crnt_pt->e1;
    eigns[1] = crnt_pt->e2;
    unsigned int x0 = crnt_pt->x;
    unsigned int y0 = crnt_pt->y;

    unsigned int div = 100*2*Pallete::PI0*sqrt(0.5*(eigns[0]*eigns[0]+eigns[1]*eigns[1]));
    float delta_alpha = Pallete::PI0/div, alpha = 0;
    angle = crnt_pt->sori;

    for(dx = -5; dx <= 5; dx++)
    {
        x = x0 + dx;
        img->setPixel(x, y0, cred);
    }
    for(dy = -5; dy <= 5; dy++)
    {
        y = y0 + dy;
        img->setPixel(x0, y, cred);
    }
    div = 2*div;
    for(i = 0; i < div; i++)
    {
        y = eigns[1]*sin(alpha);
        x = eigns[0]*cos(alpha);

        xbar = x*cos(angle) + y*sin(angle);
        ybar = y*cos(angle) - x*sin(angle);

        xi = x0 + (unsigned int)round(xbar);
        yi = y0 + (unsigned int)round(ybar);
        img->setPixel(xi-1, yi, cblk);
        img->setPixel(xi, yi, cblk);
        img->setPixel(xi, yi-1, cblk);
        alpha += delta_alpha;
    }
}


void Pallete::lineto(CImage *Img, const unsigned int x0, const unsigned int y0,
                     const unsigned int x1, const unsigned int y1, VIREO_COLOR fcolor)
{
    float sy, sx, ey, ex;
    int i, height, width;
    float dx, dy, absval, ratio, px, py;
    int crnt_y, crnt_x;
    float mycolor[3];
    Pallete::fetchColor(fcolor, mycolor);

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
        for(i = 1; i < absval; i++)
        {

            crnt_y = ey - VMath::Sign(dy)*i;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, mycolor);
                Img->setPixel(crnt_x+1, crnt_y, mycolor);
            }
        }
        return ;
    }

    if(dy == 0)
    {
        absval = VMath::absx(dx);
        crnt_y = ey;
        for(i = 1; i < absval; i++)
        {

            crnt_x = ex - VMath::Sign(dx)*i;
            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, mycolor);
                Img->setPixel(crnt_x+1, crnt_y, mycolor);
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
        //step by y
        absval = dy >= 0?dy:(-dy);
        for(i = 1; i < absval; i++)
        {
            px = i/ratio;
            crnt_x = sx + VMath::Sign(dx)*(int)round(px);
            crnt_y = sy + VMath::Sign(dy)*i;
            if((crnt_x < width) && (crnt_y < height))
            {

                Img->setPixel(crnt_x, crnt_y, mycolor);
                Img->setPixel(crnt_x+1, crnt_y, mycolor);
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
            crnt_y = crnt_y <= 0?1:crnt_y;

            if((crnt_x < width)&&(crnt_y < height))
            {
                Img->setPixel(crnt_x, crnt_y, mycolor);
                Img->setPixel(crnt_x+1, crnt_y, mycolor);
            }
        }
    }
    return ;
}


void Pallete::buildView(const char *srcfn, const char *dstfn, vector<KeyPoint*> &kps,
                        VIREO_COLOR fcolor, KeyPtType flt)
{
    vector<KeyPoint*>::iterator itkp;
    KeyPoint* crnt_kpt;
    CImage *view = new CImage(srcfn);
    for(itkp = kps.begin(); itkp != kps.end(); itkp++)
    {
        crnt_kpt = *itkp;

        if(crnt_kpt->_type == flt)
        {
            Pallete::drawEllipse(view, crnt_kpt, fcolor);
        }
    }

    view->save(dstfn);
    delete  view;
}


void Pallete::buildEdgeView(const unsigned char* edgeimg, const int width, const int height, const char *dstfn)
{
    Image *myimg = new Image(width, height);
    int count = 0;
    for(int i = 0; i < height; i++)
    {
        for(int j = 0; j < width; j++)
        {
            myimg->pix[count] = edgeimg[count];
            count++;
        }
    }
    myimg->save(dstfn);
    delete myimg;
    myimg  = NULL;
}


void Pallete::test()
{
    float mycolor[3];
    Pallete::fetchColor(_vireo_purple, mycolor);
    cout<<"Color: "<<mycolor[0]<<"\t"<<mycolor[1]<<"\t"<<mycolor[2]<<endl;
}
