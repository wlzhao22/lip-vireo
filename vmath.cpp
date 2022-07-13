#include <iostream>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "vmath.h"

using namespace std;

/** this code for SIFT norm (from Lowe's paper) does not really work well
       sum = sqrt(sum);
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;

            if(vect[i] > 0.20)
                vect[i] = 0.20;

        }
        sum = 0;
        for(i = 0; i < dim; i++)
        {
            sum = sum + vect[i]*vect[i];
        }
        sum = sqrt(sum)/512.0f;
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;
        }
**/

void VMath::l2norm(float *vect, const int dim)
{
    assert(vect);
    assert(dim > 0);
    float sum = 0;
    int i = 0;

    for(i = 0; i < dim; i++)
    {
        sum = sum + vect[i]*vect[i];
    }
    sum = sqrt(sum);
    if(sum > 0.05)
    {
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;
        }
    }
}

int VMath::l2norm(float *vectr, const int dim, const int nseg)
{
    assert(vectr);
    int i = 0, si = 0, loci = 0, sd = dim/nseg;
    float sum = 0;

    for(si = 0; si < nseg; si++)
    {
        sum = 0;
        for(i = 0; i < sd; i++)
        {
            sum = sum + vectr[loci+i]*vectr[loci+i];
        }
        sum = sqrt(sum);
        if(sum > 0.0f)
        {
            for(i = 0; i < sd; i++)
            {
                vectr[loci+i] = vectr[loci+i]/sum;
                ///cout<<vectr[loci+i]<<" ";
            }
        }
        loci += sd;
    }
    return 1;
}

void VMath::l1norm(float *vect, const int dim)
{
    assert(vect);
    assert(dim > 0);
    float sum = 0;
    int i;
    for(i = 0; i < dim; i++)
    {
        sum = sum + fabs(vect[i]);
    }

    if(sum != 0)
    {
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;
        }
    }
}

void VMath::l1norm(float *vect, const int dim, const float factor0)
{
    float total = 0;
    int i = 0;
    assert(dim > 0);
    for(i = 0; i < dim; i++)
    {
        total += fabs(vect[i]);
    }

    total /= dim;


    for(i = 0; i < dim; i++)
    {
        vect[i] = (vect[i] / total) * factor0;
    }
}

bool VMath::SIFTNorm(float *vect, const int dim)
{
    assert(vect);

    float sum = 0;
    int i = 0;
    for(i = 0; i < dim; i++)
    {
        sum = sum + vect[i]*vect[i];
    }
    if(sum > 0.00)
    {
        /**standard way**/
        sum = sqrt(sum)/512.0;
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;
        }
        /**/
        return true;
    }
    else
    {
        return false;
    }
}

bool VMath::sqrtSIFTNorm(float *vect, const int dim)
{
    assert(vect);
    float sum = 0;
    int i;
    for(i = 0; i < dim; i++)
    {
        vect[i] = sqrt(vect[i]);
        sum += vect[i]*vect[i];
    }
    sum = sqrt(sum)/4092.0f;
    if(sum > 0.00)
    {
        for(i = 0; i < dim; i++)
        {
            vect[i] = vect[i]/sum;
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool VMath::powerLaw(float *vectr, const unsigned int d0, const float p0)
{
    assert(vectr);
    float tmp = 0;
    for(unsigned int i = 0; i < d0; i++)
    {
             tmp = fabs(vectr[i]);
        vectr[i] = VMath::Sign(vectr[i])*pow(tmp, p0);
    }
    return true;
}

float VMath::dst_cos(const float vect1[], const float vect2[], const unsigned int d0)
{
    float w1  = 0.0f, w2 = 0.0f;
    float val = 0.0f;
    unsigned int i = 0;
    for(i = 0; i < d0; i++)
    {
        w1  = w1  + vect1[i]*vect1[i];
        w2  = w2  + vect2[i]*vect2[i];
        val = val + vect1[i]*vect2[i];
    }
    w1 = sqrt(w1);
    w2 = sqrt(w2);

    if(w1 == 0 )
    {
        return 0;
    }
    else if(w2 == 0)
    {
        return 0;
    }
    else
    {
        val = val/(w1*w2);
        return val;
    }
}

double VMath::besseli(const unsigned int nu, const float x)
{
   double val   = x*0.5;
   double term  = 0.0, bessel = 0;
   double gamma = 0;
   int m = 0, p = 1, i;

   do{
    gamma = factorial(m)*factorial(nu+m);
    gamma = 1.0/gamma;
    p     = 2*m + nu;
    term  = val;
    for(i = 1; i < p; i++)
    {
        term = term*val;
    }
    term  = term*gamma;
    bessel += term;
    m++;
   }while(term > 1e-7);

    return bessel;
}

double VMath::factorial(const unsigned int n)
{
   unsigned int i = 0;
   double result = 1;
   if(n == 0)
   return 1;

   for(i = 1; i <= n; i++)
   {
      result = result*i;
   }
   return result;
}


float VMath::lgx(const float a, const float base)
{

    return log(a)/log(base);

}

/**
 * Calculates the inverse of a 3x3 matrix.
 *
 * @param a Input matrix
 * @param b Output inverse matrix.
 */
void VMath::mInv33(float a[3][3], float b[3][3])
{
    assert(a);
    assert(b);

    float det = a[0][0]*(a[2][2]*a[1][1] - a[2][1]*a[1][2])
                - a[1][0]*(a[2][2]*a[0][1] - a[2][1]*a[0][2])
                + a[2][0]*(a[1][2]*a[0][1] - a[1][1] * a[0][2]);

    b[0][0] = (a[2][2]*a[1][1] - a[2][1]*a[1][2]) / det;
    b[0][1] = (-a[2][2]*a[0][1] + a[2][1]*a[0][2]) / det;
    b[0][2] = (a[1][2]*a[0][1] - a[1][1]*a[0][2]) / det;

    b[1][0] = (-a[2][2]*a[1][0] + a[2][0]*a[1][2]) / det;
    b[1][1] = (a[2][2]*a[0][0] - a[2][0]*a[0][2]) / det;
    b[1][2] = (-a[1][2]*a[0][0] + a[1][0]*a[0][2]) / det;

    b[2][0] = (a[2][1]*a[1][0] - a[2][0]*a[1][1]) / det;
    b[2][1] = (-a[2][1]*a[0][0] + a[2][0]*a[0][1]) / det;
    b[2][2] = (a[1][1]*a[0][0] - a[1][0]*a[0][1]) / det;
}

/**
 * Finds the x coordinate of the maximum (or minimum) of the parabola described by 3 points.
 *
 * @return x coordinate of the maximum (or minimum) of the parabola.
 */
float VMath::getMaxParabola(float x1, float y1, float x2, float y2,float x3, float y3)
{

    float mb[3] = {y1, y2, y3};

    float mA[3][3] = {{x1*x1, x1, 1},
        {x2*x2, x2, 1},
        {x3*x3, x3, 1}
    };
    float mAInv[3][3];

    mInv33(mA, mAInv);

    float a = 0;

    // mx = maInv * mb
    for (int i = 0; i < 3; i++)
        a += mAInv[0][i] * mb[i];

    float b = 0;
    for (int i = 0; i < 3; i++)
        b += mAInv[1][i] * mb[i];

    return -b / (2*a);
}

/**8
Input 2 by 2 matrix
a: m11
b: m12 and m21
c: m22

output two eigenvalues
eigen1 >= eigen2
**/

bool VMath::eigvlSymMat(const float a,const float b, const float c, float &eign1, float &eign2)
{
    float tmp = a + c;
    float delta = sqrt(a*a+c*c+4*b*b-2*a*c);
    eign2 = (tmp + delta)/2;
    eign1 = (tmp - delta)/2;
    return true;
}

bool VMath::eigvtSymMat(const float a, const float b, const float c, float eignv[2][2], float eigs[2])
{
    ///in column order for eigenvector
    float tmp = a + c;
    float delta = sqrt(a*a+c*c+4*b*b-2*a*c);
    float eign1 = (tmp + delta)/2;
    float eign2 = (tmp - delta)/2;

    eigs[0] = eign1;
    eigs[1] = eign2;
    ///cout<<eign1<<"\t"<<eign2<<endl;

    if(fabs(b) < 0.0001)
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
    }
    else
    {
        float len1, len2, ratio;
        //get eigen vectors for the normal case
        eignv[0][0] = eign1 - c;
        eignv[0][1] = eign2 - c;
        eignv[1][0] = eignv[1][1] = b;
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize

        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    return true;
}



bool VMath::sqrtSymMat(const float a, const float b, const float c, float U[2][2], float eigs[2])
{
    ///in column order for eigenvector
    float tmp = a + c;
    float delta = sqrt(a*a+c*c+4*b*b-2*a*c);
    float eign1 = (tmp + delta)/2;
    float eign2 = (tmp - delta)/2;
    float eignv[2][2] = {0};
    float m11 = 0, m12 = 0, m21 = 0, m22 = 0;

    eigs[0] = eign1;
    eigs[1] = eign2;

    if(fabs(b) < 0.0001)
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
    }
    else
    {
        float len1, len2, ratio;

        ///get eigen vectors for the normal case
        eignv[0][0] = eign1 - c;
        eignv[0][1] = eign2 - c;
        eignv[1][0] = eignv[1][1] = b;

        ///orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);

        ///normalize
        len1 = sqrt(len1);
        m11  = eignv[0][0]/len1;
        m21  = eignv[1][0]/len1;

        m12  = eignv[0][1]/len2;
        m22  = eignv[1][1]/len2;

        eigs[0] = sqrt(fabs(eigs[0]));
        eigs[1] = sqrt(fabs(eigs[1]));
        eigs[1] = eigs[1]/eigs[0];
        eigs[0] = 1;
        eigs[1] = 1.0/eigs[1];
        ///eigs[1] = 1.0;

        U[0][0] = m11*m11*eigs[0] + m12*m12*eigs[1];
        U[0][1] = m11*m21*eigs[0] + m12*m22*eigs[1];
        U[1][0] = m11*m21*eigs[0] + m12*m22*eigs[1];
        U[1][1] = m21*m21*eigs[0] + m22*m22*eigs[1];
    }
    return true;
}



float VMath::Det(const float *mat, const unsigned int row)
{
    assert(row == 3);
    float det = 0;

    det += mat[0]*mat[4]*mat[8];
    det += mat[1]*mat[5]*mat[6];
    det += mat[2]*mat[3]*mat[7];

    det -= mat[2]*mat[4]*mat[6];
    det -= mat[0]*mat[5]*mat[7];
    det -= mat[1]*mat[3]*mat[8];

    return det;
}

/****************************************************
@return inversed square root of input the symmetric matrix M
****************************************************/

bool VMath::inv_sqrtSymMat(const float a,const float b, const float c, float mat[2][2], const float thr0)
{
    float eignv[2][2];
    float tmp, delta, eign1, eign2;
    tmp = a + c;
    delta = sqrt(a*a+c*c+4*b*b-2*a*c);

    if(a < c)
    {
        eign2 = (tmp + delta)/2;
        eign1 = (tmp - delta)/2;
    }
    else
    {
        eign1 = (tmp + delta)/2;
        eign2 = (tmp - delta)/2;
    }

    if(eign1 <= 0 || eign2 <= 0)
    {
        mat[0][0] = mat[1][1] = 1;
        mat[0][1] = mat[1][0] = 0;
        return true;
    }

    if(b == 0)
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
    }
    else
    {
        //get eigen vectors for the normal case
        float len1, len2, ratio;

        eignv[0][0] = eign1 - c;
        eignv[0][1] = eign2 - c;
        eignv[1][0] = eignv[1][1] = b;
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }

    //get square root inverse of eigenvalues
    eign1 = 1.0/sqrt(eign1);
    eign2 = 1.0/sqrt(eign2);
    tmp = eign1> eign2?(eign1/eign2):(eign2/eign1);

    //get inversed square root of the matrix by combining three parts v*d*v'
    mat[0][0] = eign1*eignv[0][0]*eignv[0][0] + eign2*eignv[0][1]*eignv[0][1];
    mat[0][1] = eign1*eignv[0][0]*eignv[1][0] + eign2*eignv[0][1]*eignv[1][1];
    mat[1][0] = eign1*eignv[1][0]*eignv[0][0] + eign2*eignv[1][1]*eignv[0][1];
    mat[1][1] = eign1*eignv[1][0]*eignv[1][0] + eign2*eignv[1][1]*eignv[1][1];

    if(tmp >= thr0)
    {
        return true;
    }
    else
        return false;
}


void VMath::inv_sqrtSymMat(const float a, const float b, const float c, float mat[2][2], float ei[2], bool _norm_)
{
    float eignv[2][2] = {0, 0, 0, 0};
    float len1 = 0, len2 = 0, ratio = 0;
    float tmp = 0, delta = 0, eign1 = 0, eign2 = 0;
    float fa = a, fb = b, fc =c ;

    if(a < 0 && c < 0)
    {
        fa = abs(a);
        fc = abs(c);
        ///fb = -1*b;
    }
    else if(a > 0 && c > 0)
    {
        fb = -1*b;
    }

    tmp = fa + fc;
    delta = sqrt(fa*fa+fc*fc+4*fb*fb-2*fa*fc);

    eign1 = (tmp + delta)/2;
    eign2 = (tmp - delta)/2;

    if((eign1*eign2) < 0)
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
        return ;
    }

    //get eigen vectors for the normal case
    eignv[0][0] = eign1 - fc;
    eignv[0][1] = eign2 - fc;
    eignv[1][0] = eignv[1][1] = fb;

    //orthogonalize
    len1  = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
    ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
    ratio = ratio/len1;
    eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
    eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
    len2  = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);

    //normalize
    len1 = sqrt(len1);
    eignv[0][0] = eignv[0][0]/len1;
    eignv[1][0] = eignv[1][0]/len1;
    eignv[0][1] = eignv[0][1]/len2;
    eignv[1][1] = eignv[1][1]/len2;

    //get square root inverse of eigenvalues
    ei[0] = eign1;
    ei[1] = eign2;
    eign1 = 1.0/eign1;
    eign2 = 1.0/eign2;

    if(_norm_)
    {
        if(eign1 < eign2)
        {
            eign1 = eign1/eign2;
            eign2 = 1;
        }
        else
        {
            eign2 = eign2/eign1;
            eign1 = 1;
        }
    }

    //get inversed square root of the matrix by combining three parts v*d*v'
    mat[0][0] = eign1*eignv[0][0]*eignv[0][0] + eign2*eignv[0][1]*eignv[0][1];
    mat[0][1] = eign1*eignv[0][0]*eignv[1][0] + eign2*eignv[0][1]*eignv[1][1];
    mat[1][0] = eign1*eignv[1][0]*eignv[0][0] + eign2*eignv[1][1]*eignv[0][1];
    mat[1][1] = eign1*eignv[1][0]*eignv[1][0] + eign2*eignv[1][1]*eignv[1][1];

    return ;
}


bool VMath::normMat(const float m[2][2], float nmat[2][2])
{
    float eignv[2][2];
    float tmp, delta, eign1, eign2;
    float len1, len2, ratio;

    tmp = m[0][0] + m[1][1];
    delta = sqrt(m[0][0]*m[0][0]+m[1][1]*m[1][1]+4*m[1][0]*m[0][1]-2*m[0][0]*m[1][1]);
    if(m[0][0] < m[1][1])
    {
        eign2 = (tmp + delta)/2;
        eign1 = (tmp - delta)/2;
    }
    else
    {
        eign1 = (tmp + delta)/2;
        eign2 = (tmp - delta)/2;
    }

    if (m[1][0] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[0][0] = eign1 - m[1][1];
        eignv[0][1] = eign2 - m[1][1];
        eignv[1][0] = eignv[1][1] = m[1][0];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[1][0] = eign1 - m[0][0];
        eignv[1][1] = eign2 - m[0][0];
        eignv[0][0] = eignv[0][1] = m[0][1];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] == 0 && m[1][0] == 0 )
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
    }

    /**
    VMath::printMat(eignv);
    cout<<eign1<<"\t"<<eign2<<endl;
    **/

    //normalize eigenvalues by dividing it with the larger one
    if(eign1 < eign2)
    {
        eign1 = eign1/eign2;
        eign2 = 1;
    }
    else
    {
        eign2 = eign2/eign1;
        eign1 = 1;
    }

    //normalize eigenvalues by dividing it with the larger one
    //ratio = eign1*eign2;
    //eign1 = eign1/ratio;
    //eign2 = eign2/ratio;

    //get inversed square root of the matrix by combining three parts v*d*v'
    nmat[0][0] = eign1*eignv[0][0]*eignv[0][0] + eign2*eignv[0][1]*eignv[0][1];
    nmat[0][1] = eign1*eignv[0][0]*eignv[1][0] + eign2*eignv[0][1]*eignv[1][1];
    nmat[1][0] = eign1*eignv[1][0]*eignv[0][0] + eign2*eignv[1][1]*eignv[0][1];
    nmat[1][1] = eign1*eignv[1][0]*eignv[1][0] + eign2*eignv[1][1]*eignv[1][1];
    return true;
}


bool VMath::normMat(const float m[2][2], float nmat[2][2], const float  thr0)
{
    float eignv[2][2];
    float tmp, delta, eign1, eign2;
    float len1, len2, ratio;

    tmp = m[0][0] + m[1][1];
    delta = sqrt(m[0][0]*m[0][0]+m[1][1]*m[1][1]+4*m[1][0]*m[0][1]-2*m[0][0]*m[1][1]);
    if(m[0][0] < m[1][1])
    {
        eign2 = (tmp + delta)/2;
        eign1 = (tmp - delta)/2;
    }
    else
    {
        eign1 = (tmp + delta)/2;
        eign2 = (tmp - delta)/2;
    }
    if (m[1][0] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[0][0] = eign1 - m[1][1];
        eignv[0][1] = eign2 - m[1][1];
        eignv[1][0] = eignv[1][1] = m[1][0];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[1][0] = eign1 - m[0][0];
        eignv[1][1] = eign2 - m[0][0];
        eignv[0][0] = eignv[0][1] = m[0][1];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] == 0 && m[1][0] == 0 )
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
    }

    //get square root inverse of eigenvalues
    //eign1 = eign1/eign2;
    //eign2 = 1;

    ratio = eign1*eign2;
    eign1 = eign1/ratio;
    eign2 = eign2/ratio;

    //get inversed square root of the matrix by combining three parts v*d*v'
    nmat[0][0] = eign1*eignv[0][0]*eignv[0][0] + eign2*eignv[0][1]*eignv[0][1];
    nmat[0][1] = eign1*eignv[0][0]*eignv[1][0] + eign2*eignv[0][1]*eignv[1][1];
    nmat[1][0] = eign1*eignv[1][0]*eignv[0][0] + eign2*eignv[1][1]*eignv[0][1];
    nmat[1][1] = eign1*eignv[1][0]*eignv[1][0] + eign2*eignv[1][1]*eignv[1][1];

    ratio = (eign2 > eign1)?(eign2/eign1):(eign1/eign2);

    if(ratio > thr0)
    {
        return false;
    }
    else
    {
        return true;
    }

}


bool VMath::normMat(const float m[2][2], float nmat[2][2],  float &e1, float &e2, const float  thr0)
{
    float eignv[2][2];
    float tmp, delta, eign1, eign2;
    float len1, len2, ratio;

    tmp = m[0][0] + m[1][1];
    delta = sqrt(m[0][0]*m[0][0]+m[1][1]*m[1][1]+4*m[1][0]*m[0][1]-2*m[0][0]*m[1][1]);
    if(m[0][0] < m[1][1])
    {
        eign2 = (tmp + delta)/2;
        eign1 = (tmp - delta)/2;
    }
    else
    {
        eign1 = (tmp + delta)/2;
        eign2 = (tmp - delta)/2;
    }
    if (m[1][0] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[0][0] = eign1 - m[1][1];
        eignv[0][1] = eign2 - m[1][1];
        eignv[1][0] = eignv[1][1] = m[1][0];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] != 0 )
    {
        //get eigen vectors for the normal case
        eignv[1][0] = eign1 - m[0][0];
        eignv[1][1] = eign2 - m[0][0];
        eignv[0][0] = eignv[0][1] = m[0][1];
        //orthogonalize
        len1 = eignv[0][0]*eignv[0][0]+eignv[1][0]*eignv[1][0];
        ratio = eignv[0][0]*eignv[0][1]+eignv[1][0]*eignv[1][1];
        ratio = ratio/len1;
        eignv[0][1] = eignv[0][1] - ratio*eignv[0][0];
        eignv[1][1] = eignv[1][1] - ratio*eignv[1][0];
        len2 = sqrt(eignv[0][1]*eignv[0][1] + eignv[1][1]*eignv[1][1]);
        //normalize
        len1 = sqrt(len1);
        eignv[0][0] = eignv[0][0]/len1;
        eignv[1][0] = eignv[1][0]/len1;

        eignv[0][1] = eignv[0][1]/len2;
        eignv[1][1] = eignv[1][1]/len2;
    }
    else if(m[0][1] == 0 && m[1][0] == 0 )
    {
        eignv[0][0] = eignv[1][1] = 1.0f;
        eignv[0][1] = eignv[1][0] = 0.0f;
        e1 = e2  = 1;
    }

    //get square root inverse of eigenvalues

    ratio = eign1*eign2;
    e1 = eign1 = eign1/ratio;
    e2 = eign2 = eign2/ratio;

    //get inversed square root of the matrix by combining three parts v*d*v'
    nmat[0][0] = eign1*eignv[0][0]*eignv[0][0] + eign2*eignv[0][1]*eignv[0][1];
    nmat[0][1] = eign1*eignv[0][0]*eignv[1][0] + eign2*eignv[0][1]*eignv[1][1];
    nmat[1][0] = eign1*eignv[1][0]*eignv[0][0] + eign2*eignv[1][1]*eignv[0][1];
    nmat[1][1] = eign1*eignv[1][0]*eignv[1][0] + eign2*eignv[1][1]*eignv[1][1];

    ratio = (eign2 > eign1)?(eign2/eign1):(eign1/eign2);
    if(ratio > thr0)
    {
        return false;
    }
    else
    {
        return true;
    }
}


void VMath::normMat(vector<vector<float> > & mat)
{

    float sum = 0;

    for (unsigned int j = 0; j < mat.size(); j++)
    {
        for (unsigned int i = 0; i < mat[j].size(); i++)
        {
            assert(mat[j].size() == mat[0].size());
            sum += mat[j][i];
        }
    }

    for (unsigned int j = 0; j < mat.size(); j++)
    {
        for (unsigned int i = 0; i < mat[j].size(); i++)
        {
            mat[j][i] /= sum;
        }
    }

}

bool VMath::invSymMat(const float a,const float b, const float c, float mat[2][2])
{
    return true;
}

float *VMath::project(float *vect, const int dim, const float *prjMat,
                      const int row, const int col, const int dstdim)
{
    assert(vect);
    assert(dim == col);
    assert(dstdim > 0);

    int i = 0, k = 0, loc = 0;
    float *prjVct = new float[dstdim];
    memset(prjVct, 0, sizeof(float)*dstdim);

    for(i = 0; i < dim; i++)
    {
        loc = i*dim;
        for(k = 0; k < dim; k++)
        {
            prjVct[i] += prjMat[loc + k]*vect[k];
        }
    }

    return prjVct;
}

void VMath::printMat(const float mat[2][2])
{
    cout<<"--------------------------------------\n";
    cout<<mat[0][0]<<"\t"<<mat[0][1]<<endl;
    cout<<mat[1][0]<<"\t"<<mat[1][1]<<endl;
    cout<<"--------------------------------------\n";
}

int VMath::Sign(const float val)
{
    if(val >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

float VMath::maxVec(const float * v, const int len)
{
    assert(v);

    float maxval = v[0];
    for (int j = 1; j < len; j++)
    {
        if (v[j] > maxval)
        {
            maxval = v[j];
        }
    }
    return maxval;
}

/**
 * Mod function that returns all positive values.
 *
 * @param x value to take mod
 * @param b modulo
 * @return positive value of x%b.
 */
int VMath::mod(const int x, const int b)
{
    int val = x % b;
    if (val < 0)
        val += b;
    return val;
}

void VMath::printVect(const float *vect, const int dim)
{
    cout<<"Vector: ";
    for(int i  = 0; i < dim; i++)
    {
        cout<<vect[i]<<" ";
    }
    cout<<endl;
}

void VMath::SmoothHistogram(float *x,const int n)
{
    assert(x);
    if(n == 1)
        return ;
    float temp = 0, result;
    int prev, next;

    for ( int i = 0; i < n; i ++)
    {
        prev = (i-1) < 0? n-1: i-1;
        next =  (i+1)%n;
        result = (x[prev]+ 5*x[i]  + x[next])/7.0;
        x[prev] = temp;
        temp = result;
    }

    x[n-1] = temp;
}

void VMath::smoothHist(float *vector,const int n)
{
    assert(vector);
    if(n == 1)
        return ;
    float temp = 0, result;
    int prev, next;
    temp = vector[n-1];
    for ( int i = 0; i < n; i ++)
    {
        prev = (i-1) < 0? n-1: i-1;
        next =  (i+1)%n;
        result = (vector[prev]+ vector[i]  + vector[next])/3.0;
        vector[prev] = temp;
        temp = result;
    }

    vector[n-1] = temp;
}


void VMath::smoothonePass(float *vector,const int n)
{
    assert(vector);
    if(n == 1)
        return ;

    for ( int i = 0; i < n-1; i ++)
    {
        vector[i] = vector[i+1]+vector[i]/2;
    }

    vector[n-1] = (vector[n-1]+vector[n-2])/2;
}

void VMath::smooth(float vect[], const unsigned int idxS, const unsigned int idxE)
{
    float rslt, tmp;
    rslt  = vect[idxS];
    for(unsigned int i = idxS; i < idxE; i++)
    {
        tmp       = (vect[i-1] + 2*vect[i] + vect[i+1])/4.0f;
        vect[i-1] = rslt;
        rslt      = tmp;
    }
}

float VMath::absx(const float val)
{
    float tmpval = val<0?(-1*val):val;
    return tmpval;
}

/** convert int to hex number **/
void VMath::dec2hex(int number, int* hex_array)
{
    for(int i=0; i<4; i++)
    {
        hex_array[i] = number%256;
        number = number/256;
    }
}

int  VMath::min(const int x, const int y)
{
    return (x>y?y:x);
}

int  VMath::max(const int x, const int y)
{
    return (x>y?x:y);
}

void VMath::test()
{
    float eigv[2][2] = {0};
    float eigs[2] = {0};
    VMath::eigvtSymMat(0.476, 0.1587, 0.952, eigv, eigs);
    cout<<eigs[0]<<"\t"<<eigs[1]<<endl;
    cout<<eigv[0][0]<<"\t"<<eigv[0][1]<<endl;
    cout<<eigv[1][0]<<"\t"<<eigv[1][1]<<endl;
    return ;
}
