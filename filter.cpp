#include "filter.h"

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>

/***********************************************
 * Normalizes a matrix to unit length.
 *
 * @param mat Matrix to normalize
 **********************************************/

void Filter::SmoothHistogram(float *hist, const int n)
{
	assert(hist);
	float temp = 0, result;
	int prev, next;

	for ( int i = 0; i < n; i ++)
	{
		prev = (i-1) < 0? n-1: i-1;
		next =  (i+1)%n;
		result = (hist[prev] + hist[i]  + hist[next])/3.0;
		hist[prev] = temp;
		temp = result;
	}

	hist[n-1] = temp;
}

void Filter::SmoothHistogram(float *hist, const int n,const int start)
{
	assert(hist);
	float temp = 0, result;
	int prev, next, i;

	for (i = 0; i < n; i ++)
	{
		prev = (i-1) < 0? (n-1): i-1;
		prev = prev + start;
		next =  start + ((i+1)%n);
		result = (hist[prev] + hist[i]  + hist[next])/3.0;
		hist[prev] = temp;
		temp = result;
	}

	hist[start+n-1] = temp;
}

float Filter::Average(float *hist, int n)
{
    assert(hist);
    float ave = 0.0;
    for(int i = 0; i < n; i++)
    {
        ave += hist[i];
    }

    ave = ave/n;
    return ave;

}

void Filter::normalizeMat(vector<vector<float> > & mat)
{
	float sum = 0;
	unsigned int i, j;
	for (j = 0; j < mat.size(); j++)
	{
		for (i = 0; i < mat[j].size(); i++)
		{
			assert(mat[j].size() == mat[0].size());
			sum += mat[j][i];
		}
	}

	for (j = 0; j < mat.size(); j++)
	{
		for (i = 0; i < mat[j].size(); i++)
		{
			mat[j][i] /= sum;
		}
	}
}

/**
 * Normalizes a vector to unit length.
 *
 * @param vec vector to normalize
 */
void Filter::normalizeVec(vector<float> & vec)
{

	float sum = 0;
	for (unsigned int i = 0; i < vec.size(); i++)
		sum += vec[i];

	for (unsigned int i = 0; i < vec.size(); i++)
		vec[i] /= sum;
}

/**
 * Creates a 1D Gaussian kernel.
 * Kernel is always odd.
 *
 * @param sigma Sigma of Gaussian weighting function.
 * @return 1D Gaussian kernel
 */

vector<float> Filter::GaussianKernel1D(const float sigma)
{
    assert (sigma > 0);

	//ranges from [-3.0*sigma,+3.0*sigma];

	int dim = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));
	dim = dim*2+1;

	vector<float> kern(dim);

	float s2 = sigma * sigma;

	int c = dim / 2;

	float v;

	for (int i = 0; i < (dim + 1) / 2; i++)
	{
		v = (1.0 / (sqrt(2 * PI) * sigma) )* exp(-(i*i) / (2 * s2));
		kern[c+i] = v;
		kern[c-i] = v;

	}

	normalizeVec(kern);
	return kern;
}

void Filter::GaussianKernel1D(const float sigma, const float ratio, vector<float> &kern)
{
    assert (sigma > 0);
    kern.clear();
	int dim = (int) max(1.0f, (float)floor(ratio*sigma + 0.5f));
	dim = dim*2+1;
	for(int i = 0; i < dim; i++)
	kern.push_back(0);

	float s2 = sigma * sigma;
	int c = dim / 2;
	float v;

	for (int i = 0; i < (dim + 1) / 2; i++)
	{
		v = (1.0 / (sqrt(2 * PI) * sigma) )* exp(-(i*i) / (2 * s2));
		kern[c+i] = v;
		kern[c-i] = v;

	}
	normalizeVec(kern);
	///cout<<"Kernel size2: "<<kern.size()<<endl;
	return ;
}

/**
 * Creates a 1D Gaussian derivative kernel.
 * Kernel is always odd.
 *
 * @param sigma Sigma of Gaussian weighting function.
 * @return 1D Gaussian kernel
 */
vector<float> Filter::GaussianDxKernel1D(const float sigma)
{
    assert (sigma > 0);
	int dim = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));

	/**make odd*/
	dim = dim*2+1;

	vector<float> kern(dim);

	float v, s2 = sigma * sigma;
	int   i, c = dim / 2;
	for (i = -c; i <= c; i++)
	{
	    v = -2.0*(i*exp(-i*i/(2*s2))) /(sqrt(2*PI)*sigma*s2);
		kern[c+i] = v;
	}

	return kern;
}


/**
*creat a kernel of dxx with specified sigma
*
*@param sigma derivative parameter for gaussian convulation
*
**/

vector<float> Filter::GaussianDxxKernel1D(const float sigma)
{
    assert (sigma > 0);

	//ranges from [-3.0*sigma,+3.0*sigma];
	int dim = (int) max(1.0f, (float)floor(3*sigma + 0.5));
	dim = dim*2+1;
	vector<float> kern(dim);
	float v, s2 = sigma * sigma;
	int   c = dim / 2;

	for (int i = -c; i <= c; i++)
	{
	    v = (i*i/s2-1)*exp(-i*i/(2*s2))/s2;
		kern[c+i] = v;
	}

	return kern;

}

vector<vector<float> > Filter::GaussianDxyKernel2D(const float sigma)
{
    assert (sigma > 0);
	//ranges from [-3.0*sigma,+3.0*sigma];

	int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));
    int dim = c*2 + 1;
	vector<vector<float> > mat;
	for (int i = 0; i < dim; i++)
	{
		vector<float> row(dim);
		mat.push_back(row);
	}

	float s2 = sigma * sigma;
	float v;
	int y,x;

	for (y = -c; y <= c; y++)
	{
		for (x = -c; x <= c; x++)
		{
			v =  1.0/(s2*s2) * exp(-(x*x + y*y) / (2 * s2));
			v = v*x*y;
			mat[c+y][c+x] = v;
		}
	}
	return mat;
}

/**
 * Creates a 2D Gaussian kernel.
 * Kernel is always odd.
 *
 * @param sigma Sigma of Gaussian weighting function.
 * @return 2D Gaussian kernel
 */
vector<vector<float> > Filter::GaussianKernel2D(const float sigma)
{
	assert (sigma > 0);
	//ranges from [-3.0*sigma,+3.0*sigma];
	int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));
    int dim = c*2 + 1;
	vector<vector<float> > mat;
	for (int i = 0; i < dim; i++)
	{
		vector<float> row(dim);
		mat.push_back(row);
	}

	float v, s2 = sigma * sigma;
	int x, y;

	for (y = -c; y <= c; y++)
	{
		for (x = -c; x <= c; x++)
		{
			v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
			mat[c+y][c+x] = v;
		}
	}

	normalizeMat(mat);
	return mat;
}

/**
 * Creates a 2D Gaussian kernel.
 * Kernel is always odd.
 *
 * @param sigma Sigma of Gaussian weighting function.
 * @param radius Radius of Gassian window
 * @return 2D Gaussian kernel
 */

vector<vector<float> > Filter::GaussianKernel2D(const float sigma, const int radius)
{
	assert (sigma > 0);
	int c = radius;
    int dim = c*2 + 1;
	vector<vector<float> > mat;
	for (int i = 0; i < dim; i++)
	{
		vector<float> row(dim);
		mat.push_back(row);
	}

	float v = 0, s2 = sigma * sigma;
	int x = 0, y = 0;

	for (y = -c; y <= c; y++)
	{
		for (x = -c; x <= c; x++)
		{
			v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
			mat[c+y][c+x] = v;
		}
	}

	normalizeMat(mat);
	return mat;
}

void Filter::GaussianKernel2D(const float sigma, vector<vector<float> > &G2)
{
	assert (sigma > 0);

	int c = (int) max(1.0f, (float)floor(3.0*sigma + 0.5f));
	int dim = c*2 + 1;

	for (int i = 0; i < dim; i++)
	{
		vector<float> row(dim);
		G2.push_back(row);
	}

	float s2 = sigma * sigma;
	float v;
	int x,y;

	for (y = -c; y <= c; y++)
	{
		for (x = -c; x <= c; x++)
		{
			v = 1 / (2 * PI * s2) * exp(-(x*x + y*y) / (2 * s2));
			G2[c+y][c+x] = v;
		}
	}
	normalizeMat(G2);
	return ;
}

/**
 * Creates a 2D Gaussian kernel.
 * Kernel is always odd.
 *
 * @param sigma Sigma of Gaussian weighting function.
 * @return 2D Gaussian kernel
 */
float *Filter::GaussianKernel2D(const float sigma, int &size)
{
	int dim = (int) max(3.0f, GAUSSKERN * sigma);
	float *mat;
	// make dim odd
	if (dim % 2 == 0)
		dim++;


      mat = new float[dim*dim];

	float s2 = sigma * sigma;

	int c = dim / 2;
	int locx1, locx2, locy1, locy2;
	float v, sum = 0;

	//cout<<"In Gaussian Kernel begin"<<endl;
	for (int i = 0; i < (dim + 1) / 2; i++)
	{
	      locy1 = dim*(c+i);
	      locy2 = dim*(c-i);
		for (int j = 0; j < (dim + 1) / 2; j++)
		{
		      locx1 = c+j;
		      locx2 = c-j;
			v = 1 / (2 * PI * s2) * exp(-(i*i + j*j) / (2 * s2));
			mat[locy1+locx1] = v;
			mat[locy2+locx1] = v;
			mat[locy1+locx2] = v;
			mat[locy2+locx2] = v;

		}
	}
	//cout<<"2"<<endl;
	for (int j = 0; j < dim; j++)
	{
	      locy1 = j*dim;
		for (int i = 0; i < dim; i++)
		 {
			sum += mat[locy1+i];
		}
	}
	//cout<<"3"<<endl;
	for (int j = 0; j < dim; j++)
	{
	      locy1 = j*dim;
		for (int i = 0; i < dim; i++)
		{
			mat[locy1+i] /= sum;
		}
	}
	size = dim;
	return mat;
}


/** Convolves an entire image along the Y direction.
 *
 * @param kern Kernel to convolve with.
 * @param src Source image
 * @param dst Destination image
 **********************************************/

float Filter::fast_expn (const float x)
{
  float a, b, r, xi;
  int i;
  if (x > 25.0f) return 0.0 ;

  xi = x*256.0 / 25.0f;
  i  = (int)floor(xi);
  r  = xi - i ;
  //a  = expn_tab [i    ] ;
  //b  = expn_tab [i + 1] ;
  a = b = 1.0f;
  return a + r * (b - a) ;
}

void Filter::Convolve1DHeight(vector<float> & kernel, Image * src, Image * dst)
{
    float pixel = 0.0f, tmp = 0;
    unsigned int cen = kernel.size() / 2;
    unsigned int sz = kernel.size();
    int x = 0, y = 0, row = 0;
    unsigned int dim = 0;

	for (y = 0; y < src->height; y++)
	{
		for (x = 0; x < src->width; x++)
		{
		    pixel = 0.0;
		    for (dim = 0; dim < sz; dim++)
		    {
		        row = y + (dim - cen);
		        tmp = src->getPixel(x, row);
		        pixel += kernel[dim] * tmp;
            }

			dst->setPixel(x, y, pixel);
		}
	}
	return ;
}

/** Convolves an entire image along the X direction.
 *
 * @param kern Kernel to convolve with.
 * @param src Source image
 * @param dst Destination image
 **/

void Filter::Convolve1DWidth(vector<float> & kernel, Image * src, Image * dst)
 {
    float pixel = 0.0f, tmp = 0.0f;
    unsigned int cen = kernel.size() / 2;
    unsigned int sz = kernel.size();
    int col = 0, x = 0, y = 0;
    unsigned int dim = 0;

    for(y = 0; y < src->height; y++)
    {
        for(x = 0; x < src->width; x++)
		{
		    pixel = 0.0;
		    for (dim = 0; dim < sz; dim++)
		    {
		        col = x + (dim - cen);
		        tmp = src->getPixel(col, y);
		        pixel += kernel[dim] * tmp;
            }
			dst->setPixel(x, y, pixel);
		}
	}
	return ;
}


/**
 * Does Gaussian smoothing of an image.
 *
 * @param src Source image
 * @param dst Destination image
 * @param sigma Sigma of Gaussian weighting function.
 *
 void BlurImage(Image * src, Image * dst, float sigma);
 */
void Filter::BlurImage(Image * src, Image * dst, const float sigma)
 {
	assert (src && dst);

	assert(src->width == dst->width);
	assert(src->height == dst->height);

	Image * tmpImage = new Image(src->width, src->height);

	vector<float> convkernel = GaussianKernel1D(sigma);

	Convolve1DWidth(convkernel, src, tmpImage);
	Convolve1DHeight(convkernel,tmpImage,dst);
	convkernel.clear();

	delete tmpImage;
}

void Filter::BlurImage(Image * src, const float sigma)
 {
	assert (src);
	Image * tmpImage = new Image(src->width, src->height);
	vector<float> convkernel = GaussianKernel1D(sigma);

	Convolve1DWidth(convkernel, src, tmpImage);
	Convolve1DHeight(convkernel,tmpImage,src);
	convkernel.clear();

	delete tmpImage;
}


/**
*Scaling pixel value of one image to [0,1.0]
*
*@param Image src, source image to be scaled
* Note that it only make sense in some cases,
*not all images under this scaling operation are valid or meaningful
*/

void Filter::ScaleImage(Image *src)
{
    float max = 0;
    float min = 1;
    float scale;
    float pixel;

    int x,y;
    for(y = 0; y < src->height; y++)
    {
        for(x=0; x < src->width; x++)
        {
            pixel = src->getPixel(x,y);
            min = pixel<min?pixel:min;
            max = pixel<max?max:pixel;
        }
    }

    scale = max - min;
    if(scale == 0)
    return;

    for(y = 0; y < src->height; y++)
    {
        for(x=0; x < src->width; x++)
        {
            pixel = (src->getPixel(x,y)-min)/scale;
            pixel = 255- pixel*255;
            src->setPixel(x,y,pixel);
        }
    }

}

/**
*simply a test on Image::getConValWidth
*and on Image::getConValHeight
*/

Image* Filter::getDx(const Image *srcimg)
{
    Image *dstImg = new Image(srcimg->width, srcimg->height);
    int x = 0,y = 0;
    float val = 0;
    for(y = 0; y < srcimg->height; y++)
    {
        for(x = 0; x < srcimg->width; x++)
        {
            val = srcimg->getPixel(x+1,y) - srcimg->getPixel(x,y);
            dstImg->setPixel(x,y,val);
        }
    }

    return dstImg;
}

Image* Filter::getDy(const Image *srcimg)
{
    Image *dstImg = new Image(srcimg->width,srcimg->height);
    int x,y;
    float val;
    for(x = 0; x <srcimg->width;x++)
    {
        for(y = 0; y < srcimg->height;y++)
        {
            val = srcimg->getPixel(x,y+1) - srcimg->getPixel(x,y);
            dstImg->setPixel(x,y,val);
        }
    }

    return dstImg;
}

Image* Filter::getDxx(const Image *srcimg)
{
    Image *dstImg = new Image(srcimg->width, srcimg->height);
    int x = 0, y = 0;
    float val = 0;
    for(x = 0; x < srcimg->width; x++)
    {
        for(y = 0; y < srcimg->height; y++)
        {
            val = srcimg->getPixel(x+1, y) + srcimg->getPixel(x-1, y) - 2*srcimg->getPixel(x, y);
            dstImg->setPixel(x, y, val);
        }
    }

    return dstImg;
}

Image* Filter::getDxy(const Image *srcimg)
{
    Image *dstImg = new Image(srcimg->width, srcimg->height);
    int x = 0, y = 0;
    float val = 0;
    for(x = 0; x < srcimg->width; x++)
    {
        for(y = 0; y < srcimg->height; y++)
        {
            val = srcimg->getPixel(x+1, y+1) + srcimg->getPixel(x-1, y-1);
            val = val - srcimg->getPixel(x - 1, y + 1) - srcimg->getPixel(x + 1, y - 1);
            dstImg->setPixel(x, y, val);
        }
    }

    return dstImg;
}


Image* Filter::getDyy(const Image *srcimg)
{
    Image *dstImg = new Image(srcimg->width,srcimg->height);
    int x,y;
    float val;
    for(x = 0; x <srcimg->width;x++)
    {
        for(y = 0; y < srcimg->height;y++)
        {
            val = srcimg->getPixel(x,y+1) + srcimg->getPixel(x,y-1) - 2*srcimg->getPixel(x,y);
            dstImg->setPixel(x,y,val);
        }
    }

    return dstImg;
}

bool Filter::getDxx(const float sigma,Image *srcimg,Image *dstimg)
{
    assert(dstimg);
    assert(srcimg);
    //assert(srcimg->width <= dstimg->width);
    //assert(srcimg->height <= dstimg->height);

    vector<float> dxxkern = Filter::GaussianDxxKernel1D(sigma);
    vector<float> ikern = Filter::GaussianKernel1D(2*sigma);
    vector<float> dxkern = Filter::GaussianDxKernel1D(sigma);
    float umatrix[4] = {1,0,0,1};
    for(int i = 0; i< 4;i++)
    {
        umatrix[i] =umatrix[i];
    }

    Image *tmpimg = new Image(srcimg->width,srcimg->height);

    Filter::Convolve1DHeight(dxkern,srcimg,tmpimg);
    Filter::Convolve1DWidth(ikern,tmpimg,dstimg);

    Filter::Convolve1DWidth(ikern,dstimg,tmpimg);
    Filter::Convolve1DHeight(dxkern,tmpimg,dstimg);

    //Filter::BlurImage(srcimg,tmpimg,sigma);
    //cout<<"hellow\n";

    /**
    for(y = 0; y <srcimg->height; y++)
    {
        for(x = 0; x < srcimg->width; x++)
        {
            //pixel = tmpimg->getConValWidth(x,y,dkern);
            pix1 = tmpimg->getPixel(x-1,y);
            pixel = srcimg->getPixel(x,y);
            pix2 = tmpimg->getPixel(x+1,y);
            pixel = pix1+pix2-2*pixel;
            //pixel += srcimg->getConValHeight(x,y,dkern,umatrix);
            //pixel = pixel < 0?(-1*pixel):pixel;
            //pixel = pixel*sigma;
            dstimg->setPixel(x,y,pixel);
        }
    }
    //Filter::Convolve1DHeight(ikern,tmpimg,dstimg);
    cout<<"end\n";
    **/
    delete tmpimg;
    return true;
}

bool Filter::GetPixOrientation(const float x, const float y, const Image * image, float &m,  float &theta)
{
	assert(image);

	if (x < 1 || y < 1 || x > image->width - 2 || y > image->height - 2)
		return false;

	float a = image->getPixelBI(x + 1, y) - image->getPixelBI(x - 1, y);
	float b = image->getPixelBI(x, y - 1) - image->getPixelBI(x, y + 1);

	m = sqrt(a*a +  b*b);

	theta = atan2(b, a);

	return true;
}

bool Filter::GetPixCurl(const float x, const float y, const Image * image, float &curl)
{
	assert(image);

	if (x < 1 || y < 1 || x > image->width - 2 || y > image->height - 2)
		return false;

	curl = image->getPixelBI(x - 1, y - 1) + image->getPixelBI(x + 1, y + 1)
	       - image->getPixelBI(x, y-1) - image->getPixelBI(x-1, y);
	curl = curl*2;

	return true;
}

bool Filter::GetPixOrientation(const int x, const int y, const Image * image, float &m, float &m1, float &theta)
{
	assert(image);

	if (x < 1 || y < 1 || x > image->width - 2 || y > image->height - 2)
		return false;

	float a = image->getPixel(x + 1, y) - image->getPixel(x - 1, y);
	float b = image->getPixel(x, y - 1) - image->getPixel(x, y + 1);

	m = sqrt(a*a +  b*b);
    float val = image->getPixel(x, y);

    if(val !=0)
    m1 = m/val;

	theta = atan2(b, a);

	return true;
}

bool Filter::multiply(Image *img1,Image *img2,Image *dstimg)
{
    assert(img1);
    assert(img2);
    assert(dstimg);
    assert(img1->width == img2->width);
    assert(img1->height == img2->height);
    assert(img1->width == dstimg->width);
    assert(img1->height == dstimg->height);

    float pixel = 0.0f;
    for(int j = 0; j < img1->height; j++)
    {
        for(int i = 0; i < img1->width; i++)
        {
            pixel = img1->getPixel(i,j)*img2->getPixel(i,j);
            dstimg->setPixel(i,j,pixel);
        }
    }
    return true;
}

/**
*simply a test on Image::getConVal2
*/
bool Filter::getDxy(const float sigma,Image *srcimg,Image *dstimg)
{
    int x = 0, y = 0;
    float pixel = 0.0f;
    vector<float> dkern = Filter::GaussianDxKernel1D(sigma);
    vector<float> row;

    for(y = 0; y < dstimg->height; y++)
    {
        for(x = 0; x < dstimg->width; x++)
        {
            dstimg->setPixel(x, y, pixel);
        }
        printf("%3d\t",y);
    }

      dkern.clear();
      return true;
}

void Filter::printMat(vector<vector<float> > &mat)
{
    vector<float> row;
    unsigned int x = 0, y = 0;
    for(y = 0; y < mat.size();y++)
    {
        row = mat[y];
        for(x = 0; x <row.size(); x++)
        {
            printf("%8.6f ",row[x]);
        }
        cout<<endl;
    }
}

Image *Filter::ScaleImage(Image *srcImg,const int zoomout,const float sigma,const float sigma0)
{
	assert(srcImg);
	Image * dstImg,*tmpImg;
	float isigma;

	if (zoomout)
	{
		tmpImg = Image::doubleSizeImage(srcImg);
		dstImg = new Image(tmpImg->width, tmpImg->height);
		isigma = sqrt(sigma *sigma - sigma0* sigma0 * 4);
		Filter::BlurImage(tmpImg, dstImg, isigma);
		delete tmpImg;
	} else {
		dstImg = new Image(srcImg->width, srcImg->height);
		isigma = sqrt(sigma *sigma - sigma0* sigma0);
		Filter::BlurImage(srcImg, dstImg, isigma);
	}

	return dstImg;
}

void Filter::nonMaxSuppr(Image* srcImg)
{
    int nbs[8];
    nbs[0] = -1*srcImg->width -1;
    nbs[1] = -1*srcImg->width;
    nbs[2] = -1*srcImg->width +1;
    nbs[3] = -1;
    nbs[4] = 1;
    nbs[5] = srcImg->width -1;
    nbs[6] = srcImg->width;
    nbs[7] = srcImg->width +1;
    int x, y, loc, row, nb, k;
    bool Maximum = true;
    Image *dstImg = new Image(srcImg->width, srcImg->height);

    for(y = 1; y < (srcImg->height-1); y++)
    {
        row = y*srcImg->width;
        for(x = 1; x < (srcImg->width-1); x++)
        {
            loc = row + x;
            Maximum = true;
            for(k = 0; k < 8; k++)
            {
                nb = nbs[k] + loc;
                if(srcImg->pix[loc] <= srcImg->pix[nb])
                {
                    Maximum = false;
                }
            }
            if(Maximum)
            {
                dstImg->pix[loc] = 255;//srcImg->pix[loc];
            }
        }
    }

    memcpy(srcImg->pix, dstImg->pix, sizeof(float)*srcImg->height*srcImg->width);

    delete dstImg;
}


void Filter::getLoG()
{
    const char *fn = "/home/wlzhao/datasets/vgg/lenna.png";
    const char *dstfn = "/home/wlzhao/datasets/vgg/lenna_dx.jpg";
    float sigma = 2;
    Image *myimg  = new Image(fn);
    Image *imgdxx = new Image(myimg->width,myimg->height);
    Image *imgdyy = new Image(myimg->width,myimg->height);
    Image *imglog = new Image(myimg->width,myimg->height);

    vector<float> gdx  = Filter::GaussianDxKernel1D(sigma);
    vector<float> g    = Filter::GaussianKernel1D(sigma);
    vector<float> gdxx = Filter::GaussianDxxKernel1D(sigma);

    Filter::Convolve1DWidth(gdx, myimg, imgdxx);

    Filter::ScaleImage(imgdxx);
    imgdxx->save(dstfn);
}

void Filter::clear2Dvector(vector<vector<float> > &G2)
{
	vector<vector<float> >::iterator it;
	vector<float> crntvect;
	for(it = G2.begin(); it != G2.end(); it++)
	{
		crntvect = *it;
		crntvect.clear();
	}
	G2.erase(G2.begin(),G2.end());
}

void Filter::test()
{
    const char *fn = "/home/wlzhao/datasets/vgg/lenna.png";
    const char *dstfn = "/home/wlzhao/datasets/vgg/lenna_hess.jpg";
    float sigma = 2;

    Image *myimg  = new Image(fn);
    Image *imgdxx = new Image(myimg->width,myimg->height);
    Image *imgdyy = new Image(myimg->width,myimg->height);
    Image *imgdxy = new Image(myimg->width,myimg->height);
    Image *hess = new Image(myimg->width, myimg->height);

    vector<float> gdx  = Filter::GaussianDxKernel1D(sigma);
    vector<float> g    = Filter::GaussianKernel1D(sigma);
    vector<float> gdxx = Filter::GaussianDxxKernel1D(sigma);
    float dxx, dyy, dxy;


    Filter::Convolve1DWidth(gdx, myimg, imgdxx);
    Filter::Convolve1DHeight(gdx, imgdxx, imgdxy);
    Filter::Convolve1DWidth(gdxx, myimg, imgdxx);
    Filter::Convolve1DHeight(gdxx, myimg, imgdyy);

    for(int y = 0; y < (imgdxx->height); y++)
    {
        for(int x = 0; x < (imgdxx->width); x++)
        {
            dxx = imgdxx->getPixel(x, y);
            dyy = imgdyy->getPixel(x, y);
            dxy = imgdxy->getPixel(x, y);
            dxy = dxx*dyy - dxy*dxy;
            hess->setPixel(x, y, dxy);
        }
    }


    Filter::ScaleImage(hess);
    hess->save(dstfn);
    delete imgdxx;
    delete imgdyy;
    delete imgdxy;
}
