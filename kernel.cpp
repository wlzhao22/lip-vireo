#include "kernel.h"
#include "vmath.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>

/***
constants definition
*/

const int Kernel::mu = 8;
const int Kernel::N  = 8;

using namespace std;

void Kernel::besselAngle(const float theta, const int k0, const int N, float *mapFeat)
{
    assert(mapFeat);
    memset(mapFeat, 0, (2*N+1)*sizeof(float));

    int i = 0;
    float snh0   = sinh(k0);
    float r      = (float)(VMath::besseli(0, k0) - exp(-k0))/(2*snh0);
    mapFeat[0]   = sqrt(r);
    for(i = 1; i <= N; i++)
    {
        r =  VMath::besseli(i, k0)/snh0;
        mapFeat[i]     = cos(i*theta)*sqrt(r);
        mapFeat[N+i] = sin(i*theta)*sqrt(r);
    }
}

void Kernel::besselScale(const float scale, const int k0, const int N, float *mapFeat)
{
    assert(mapFeat);
    memset(mapFeat, 0, (2*N+1)*sizeof(float));

    int i = 0;
    float snh0 = sinh(k0);
    int sg     = VMath::Sign(scale);
    float r    = (float)(VMath::besseli(0, k0) - exp(-k0))/(2*snh0);
    mapFeat[0]   = sg*sqrt(r);
    for(i = 1; i <= N; i++)
    {
        r =  VMath::besseli(i, k0)/snh0;
        mapFeat[i]     = sg*cos(i*scale)*sqrt(r);
        mapFeat[N+i] = sg*sin(i*scale)*sqrt(r);
    }
}

void Kernel::test()
{
   float a[17] = {0};
   int k0 = 8, N = 8, i = 0;
   Kernel::besselScale(-0.3, k0, N, a);
   for(i = 0; i < 17; i++)
   {
       cout<<i<<", "<<a[i]<<"\n";
   }
   cout<<" "<<i<<endl;
   cout<<endl;
}
