#ifndef KEYPOINT_H
#define KEYPOINT_H

#include <cassert>
#include <cstring>
#include <vector>

using namespace std;

enum KeyPtType {FLM, CNR, CCV};

/*************************************
flm: functionally local maximum
cnr: corner
ccv: closed curve
***************************************/

enum CornType {_EDGE = 0, _JUNC, _CORNER, _STOP};

struct CPoint{

public:
    float x, y;
    /**sidx: index on the slave **/
    /** idx: index on the list position **/
    /**cidx: indx on the cordinates (x, y), which
     is the unique identification on the point **/
    unsigned int sidx, idx, cidx;
    unsigned int slaves[8];
    /**slaves: point to slaves**/
    unsigned int neibs[8];
    /** point  to all neighbors ***/

    unsigned char nn;
    /**number of neighbors**/

    int cid;   /**id of corner, inrementally starts from 1 for corner, **/
               /**decrementally starts -1 for closed curve             **/
    unsigned int eid; /**id of edge in which the corn has been found, corner
                         identified along thesame edge will share the same id**/
    float curv;
    CornType cntype;
    CPoint *next;

    CPoint()
    {
        y      = x = 0;
        curv   = 0;
        next   = NULL;
        cntype = _EDGE;
        idx    = sidx = cidx = 0;
        nn     = 0;
        for(unsigned int i = 0; i < 8; i++)
        slaves[i] = 0;
    }

    CPoint(const unsigned int x0, const unsigned int y0)
    {
        x    = x0;   y = y0;
        curv = 0;
        idx    = sidx = cidx = 0;
        nn   = 0;
        cntype = _EDGE;
        for(unsigned int i = 0; i < 8; i++)
        slaves[i] = 0;
    }

    static int llcomp(const CPoint *a,const CPoint *b)
    {
         if(a->y < b->y)
         {
             return 1;
         }else if (a->y == b->y)
         {
             return (a->x < b->x);
         }else{
             return 0;
         }
    }
};

struct LPoint {
    float x, y;
    unsigned int cidx;
};

struct XPoint{

    float x, y;
    bool valid;
    XPoint(float x0, float y0)
    {
        x     = x0;
        y     = y0;
        valid = true;
    }
};

struct HCPoint{

    unsigned int pcount;
    CPoint       *next;
    bool visited;
    unsigned int p_cidx;
    unsigned int eid; //keep the edge id, actually it later becomes sub-graph id
    HCPoint()
    {
        pcount  = 0;
        visited = false;
        p_cidx  = 0;
        next    = NULL;
    }
};

struct PPoint{
    unsigned int pidx;
    CPoint*      parent;
};

struct Cords{
public:
    unsigned short x, y;
};


class KeyPoint
{
public:
    int x, y;
    float fx, fy;
    float sx; //for dog
    float sy; //for dog
    float a;
    float b;
    float c;
    float e1, e2; //eigen values for the ellipse
    float u1, u2, u3, u4;
    float iscale;
    float dscale; //scale before amplify
    float gscale; //for dog
    int   scale, img_width; //for dog
    float fscale; //for dog
    float funcVal;
    bool  KP;
    float octSigma;
    int octIndex;
    float ori;
    float sori; //structural orientation
    int octave;
    KeyPtType _type;
    short div, flip;

    KeyPoint()
    {
        x = y = 0;
        fx = fy = 0.0f;
        a = 1;
        b = 0;
        c = 1;
        e1 = e2 =1;
        iscale = dscale = 0;
        KP = true;
        ori  = sori = 0.0f;
        octave = 0;
        funcVal = 0.0f;
        div = flip = 0;
        img_width = 0;
    }
    ~KeyPoint()
    {
    }

    static int keypCompF (const KeyPoint *a,const KeyPoint *b)
    {
        if(a->funcVal < b->funcVal)
            return 0;
        else
            return 1;
    }

    static KeyPoint *light_clone(const KeyPoint* kp)
    {
        assert(kp);
        KeyPoint *clone_pt = new KeyPoint();
        clone_pt->x  = kp->x;
        clone_pt->y  = kp->y;
        clone_pt->fx = kp->fx;
        clone_pt->fy = kp->fy;
        clone_pt->sx = kp->sx;
        clone_pt->sy = kp->sy;
        clone_pt->a  = kp->a;
        clone_pt->b  = kp->b;
        clone_pt->c  = kp->c;
        clone_pt->e1 = kp->e1;
        clone_pt->e2 = kp->e2;
        clone_pt->KP   = kp->KP;
        clone_pt->ori  = kp->ori;
        clone_pt->sori = kp->sori;
        clone_pt->div  = kp->div;
        clone_pt->flip = kp->flip;
        clone_pt->iscale = kp->iscale;
        clone_pt->dscale = kp->dscale;
        clone_pt->gscale = kp->gscale;
        clone_pt->scale  = kp->scale;
        clone_pt->fscale = kp->fscale;
        clone_pt->octave = kp->octave;
        //clone_pt->_type  = kp->_type;
        clone_pt->img_width = kp->img_width;
        clone_pt->funcVal   = kp->funcVal;
        return clone_pt;
    }

    static KeyPoint *deep_clone(const KeyPoint* kp)
    {
        assert(kp);
        KeyPoint *clone_pt = new KeyPoint();
        clone_pt->x  = kp->x;
        clone_pt->y  = kp->y;
        clone_pt->fx = kp->fx;
        clone_pt->fy = kp->fy;
        clone_pt->sx = kp->sx;
        clone_pt->sy = kp->sy;
        clone_pt->a  = kp->a;
        clone_pt->b  = kp->b;
        clone_pt->c  = kp->c;
        clone_pt->e1 = kp->e1;
        clone_pt->e2 = kp->e2;
        clone_pt->KP   = kp->KP;
        clone_pt->ori  = kp->ori;
        clone_pt->sori = kp->sori;
        clone_pt->div  = kp->div;
        clone_pt->flip = kp->flip;
        clone_pt->iscale = kp->iscale;
        clone_pt->dscale = kp->dscale;
        clone_pt->gscale = kp->gscale;
        clone_pt->scale  = kp->scale;
        clone_pt->fscale = kp->fscale;
        clone_pt->octave = kp->octave;
        ///clone_pt->_type  = kp->_type;
        clone_pt->img_width = kp->img_width;
        clone_pt->funcVal   = kp->funcVal;
        return clone_pt;
    }
};

#endif
