#include "haar.h"
#include "filter.h"

#include <iostream>
#include <cmath>


using namespace std;

int Haar::upbound = 10;

/** The 1D Haar Transform **/
void Haar::haar1d(float *vec, int n)
{
	int i=0;
	int w=n;
	float *vecp = new float[n];
	for(i=0;i<n;i++)
		vecp[i] = 0;

	while(w > Haar::upbound)
	{
		w/=2;
		for(i=0;i<w;i++)
		{
			vecp[i] = (vec[2*i] + vec[2*i+1])/sqrt(2.0);
			vecp[i+w] = (vec[2*i] - vec[2*i+1])/sqrt(2.0);
		}

		for(i=0;i<(w*2);i++)
			vec[i] = vecp[i];
	}

	delete [] vecp;
}


/** A Modified version of 1D Haar Transform, used by the 2D Haar Transform function **/
void Haar::haar1(float *vec,const int n, const int w)
{
	int i=0;
	//cout<<"before haar1\n";
	float *vecp = new float[n];
	for(i=0;i<n;i++)
    {
		    vecp[i] = 0;
    }
    int len = w/2;
    //cout<<len<<"\t"<<n<<endl;

	for(i=0;i<len;i++)
	{
		vecp[i] = (vec[2*i] + vec[2*i+1])/sqrt(2.0);
		vecp[i+len] = (vec[2*i] - vec[2*i+1])/sqrt(2.0);
	}
	//cout<<"end of assigment\n";

	for(i=0;i<(len*2);i++)
	{
	    vec[i] = vecp[i];
	    //cout<<vec[i]<<"\t";
	}
	//cout<<endl;

    //cout<<"end haar1\n";
		delete [] vecp;
}


/** The 2D Haar Transform **/
void Haar::haar2(float *matrix, const int rows, const int cols)
{
	float *temp_row = new float[cols];
	float *temp_col = new float[rows];

	int i=0,j=0;
	int w = cols, h=rows;
	int loc;
	while(w > Haar::upbound || h > Haar::upbound)
	{
		if(w > Haar::upbound)
		{
			for(i=0;i<h;i++)
			{
			    loc = i*cols;
				for(j=0;j<cols;j++)
				{
					    temp_row[j] = matrix[loc+j];
				}

				haar1(temp_row,cols,w);

				for(j=0;j<cols;j++)
                {
					    matrix[loc+j] = temp_row[j];
				}
			}
		}


		if(h > Haar::upbound)
		{
		    //cout<<h<<"\t"<<w<<endl;
			for(i=0;i<w;i++)
			{
				for(j=0;j<rows;j++)
				{
				    loc =  j*cols;
					temp_col[j] = matrix[loc+i];
				}
				haar1(temp_col, rows, h);
				for(j=0;j<rows;j++)
				{
				    loc = j*cols;
				    matrix[loc+i] = temp_col[j];
				}
			}
		}

		if(w > Haar::upbound)
			w/=2;
		if(h > Haar::upbound)
			h/=2;
	}

	delete [] temp_row;
	delete [] temp_col;
}


void Haar::haar2(Image *srcimg)
{
    int rows = 2*(srcimg->height/2) + 1;
    int cols = 2*(srcimg->width/2) + 1;

	float *temp_row = new float[cols];
	float *temp_col = new float[rows];

	int i=0,j=0;
	int w = cols, h=rows;
	//cout<<"debug 1\n";
	//cout<<"w: "<<cols;
	//cout<<"h: "<<rows<<endl;

	while(w > Haar::upbound || h > Haar::upbound)
	{
		if(w > Haar::upbound)
		{
			for(i=0;i<h;i++)
			{
				for(j=0;j<cols;j++)
				{
					temp_row[j] = srcimg->getPixel(j,i);
				}

				haar1(temp_row,cols,w);

				for(j=0;j<cols;j++)
				{
					srcimg->setPixel(j,i,temp_row[j]);
				}
			}
		}

		//cout<<"========================\n";

		if(h > Haar::upbound)
		{
		    //cout<<w<<"\t"<<h<<endl;
			for(i=0;i<w;i++)
			{
			    //cout<<w<<"\n";
				for(j=0;j<rows;j++)
				{
				    temp_col[j] = srcimg->getPixel(i,j);
				    //cout<<j<<endl;
				}

				haar1(temp_col, rows, h);
				//cout<<rows<<"i\t"<<h<<endl;
				for(j=0;j<rows;j++)
				{
				    //cout<<j<<endl;
					srcimg->setPixel(i,j,temp_col[j]);

				}
			}

		}

		//cout<<w<<"\t"<<h<<endl;

		if(w > Haar::upbound)
			w/=2;
		if(h > Haar::upbound)
			h/=2;

			break;
	}

    //cout<<"debug 2\n";

	delete [] temp_row;
	delete [] temp_col;
}

/** Here's an example on how to use these functions **/
int Haar::test1()
{
	int i=0;
	float vec3[4] = {4,2,5,5};

	haar1d(vec3,4);

	cout << "The 1D Haar Transform: " << endl;
	for(i=0;i<4;i++)
		cout << vec3[i] << " ";
	cout << endl;

	cout << "\n\nThe 2D Haar Transform: " << endl;


	cout << endl;
	return 0;
}

int Haar::test()
{
    const char *fn = "e:/siftlab/src/graf_img1.pgm";
    const char *dstfn = "e:/siftlab/src/graf_img1_haar.pgm";

    Image *srcimg = new Image(fn);

    Haar::haar2(srcimg);
    srcimg->multiply(1);
    //Filter::ScaleImage(srcimg);

    srcimg->save(dstfn);
	cout << endl;
	return 0;
}
