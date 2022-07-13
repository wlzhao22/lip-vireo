#ifndef KERNEL_H
#define KERNEL_H


/**
Implementation of kernel mapping, in particular with modified Bessel function

Follow the idea shown in paper:
Giogos Tolias, Adnraei Bursuc, Teddy Furon, Herve Jegou: Rotation and translation covariant match kernels for image retrieval

**/

class Kernel
{
    private:
        const static int N;
        const static int mu;
    public:
        Kernel(){}
        virtual ~Kernel(){}
        static void besselAngle(const float angle, const int mu0, const int N, float *mapFeat);
        static void besselScale(const float scale, const int mu0, const int N, float *mapFeat);
        static void test();
};

#endif // KERNEL_H
