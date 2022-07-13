#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#ifndef __APPLE__
#include <malloc.h>
#endif

#include <setjmp.h>
#include <cassert>
#include <cmath>

#include "vstring.h"
#include "ioimage.h"
#include "vmath.h"
#include "bmp.h"


using namespace std;

struct my_error_mgr
{
    struct jpeg_error_mgr pub;
    /** "public" fields */
    jmp_buf setjmp_buffer;
    /** for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void)my_error_exit (j_common_ptr cinfo)
{
    my_error_ptr myerr = (my_error_ptr) cinfo->err;
    (*cinfo->err->output_message) (cinfo);
    longjmp(myerr->setjmp_buffer, 1);
}

char  *IOImage::read_pgm(const char *srcfn ,int &width, int &height,int *maxVal)
{
    int  i, int_tmp, num = 0, count = 0;
    vector<int> digits;
    width = height = 0;
    char header[512];
    bool Binary;
    char *data;
    ifstream inStrm;
    inStrm.open(srcfn, fstream::in|ifstream::binary);

    if(!inStrm.is_open())
    {
        cerr << "File '" << srcfn << "' cannot be open!\n";
        return NULL;
    }
    inStrm.getline(header, 512);
    if(!strncmp(header, "P5", 2))
    {
        Binary = true;
    }
    else if(!strncmp(header, "P2", 2))
    {
        Binary = false;
    }
    else
    {
        cerr<<"'"<<header<<"'"<<" is not valid magic number for PGM format!\n";
        exit(0);
    }

    do
    {
        if(header[0] != '#')
        {
            VString::trimAfter(header, '#');
            count = VString::parse_words(header, digits);
            num  += count;
        }
        if(num < 3)
        {
            inStrm.getline(header, 512);
        }
    }
    while(num < 3 && !inStrm.eof());

    width   = digits[0];
    height  = digits[1];
    *maxVal = digits[2];
    digits.clear();
    data = new char [width*height];

    if(Binary)
    {
        inStrm.read(data, width*height*sizeof(char));
        inStrm.close();
    }
    else
    {
        inStrm.close();
        inStrm.open(srcfn, ifstream::in);
        inStrm.getline(header, 512);
        do
        {
            if(header[0] != '#')
            {
                VString::trimAfter(header, '#');
                count = VString::parse_words(header, digits);
                num  += count;
            }
            if(num < 3)
            {
                inStrm.getline(header, 512);
            }
        }
        while(num < 3 && !inStrm.eof());

        digits.clear();
        for (i = 0; i < height*width; i++)
        {
            inStrm>>int_tmp;
            data[i]  = int_tmp;
        }
        inStrm.close();
    }
    return data;
}

char *IOImage::read_bmp(const char *fn, int &width, int &height,
                        int *maxVal, const int channel)
{
    ifstream file;
    file.open(fn, std::ios::in | std::ios::binary );

    if(file.fail())
    {
        printf("File %s can not open",fn);
        return NULL;
    }

    file.seekg(0, std::ios::end);
    long flen = file.tellg();
    file.seekg(0, std::ios::beg);
    width = height = 0;

    char* buffer = new char[flen + 1];
    file.read(buffer, flen);
    file.close();
    char *pBuff = buffer;

    BITMAPFILEHEADER* bmpheader = new BITMAPFILEHEADER;
    BITMAPINFOHEADER* bmih      = new BITMAPINFOHEADER;
    memcpy(bmpheader, pBuff, sizeof(BITMAPFILEHEADER));
    BITMAPFILEHEADER* bmfh = bmpheader;
    pBuff += sizeof(BITMAPFILEHEADER);

    /// verify that it's a BITMAP file
    if(bmfh->bfType != BITMAP_ID)
    {
        delete [] buffer;
        buffer = NULL;
        cout<<hex<<bmfh->bfType<<endl;
        cout<<"File is not standard BMP type or the program is unable to recorgnize it!\n";
        return NULL;
    }

    memcpy(bmih, pBuff, sizeof(BITMAPINFOHEADER));

    pBuff += sizeof(BITMAPINFOHEADER);

    int bitCount = bmih->biBitCount;
    int  colors = bitCount/8;
    ///int compression	= bmih->biCompression;
    unsigned char r = 0, g = 0, b = 0;
    unsigned int  pt = 0, pt3 = 0;
    *maxVal = 0;
    int i = 0, j = 0;
    width  = bmih->biWidth;
    height = bmih->biHeight;

    int jump;
    int sz = bmih->biSizeImage;
    if( sz > (colors*width*height))
    {
        jump = sz - colors*width*height;
        jump = jump/height;
    }
    else
    {
        jump = 0;
    }

    char *data = NULL;
    /// move the pixel data Pointer to the begining of bitmap data
    pBuff = buffer + (bmfh->bfOffBits);
    if(bmih->biCompression == 1)
    {
        char *rawData = IOImage::decmp_bmp(pBuff, width, height, maxVal);
        delete bmpheader;
        delete [] buffer;
        delete bmih;

        if(channel == 3)
        {
            data = new char[width*height*3];
            for(j = 0 ; j < height; j++)
            {
                pt  = width*j;
                pt3 = width*j*3;
                for(i = 0 ; i < width ; i++)
                {
                    data[pt3+i*3]   = rawData[pt+i];
                    data[pt3+i*3+1] = rawData[pt+i];
                    data[pt3+i*3+2] = rawData[pt+i];
                }
            }
            delete [] rawData;
            rawData = NULL;
            return data;
        }
        else
        {
            return rawData;
        }

    }
    else if(bmih->biCompression > 1)
    {
        cout<<"color BMP with RLE compression is not supported!\n";
        width = height = 0;
        return NULL;
    }

    if(channel == 1 || channel == 3)
    {
        data = new char[channel*width*height];
    }
    else
    {
        cout<<"Channel should be eigther 1 or 3\n";
        exit(1);
    }

    for(j = height ; j > 0  ; j--)
    {
        pt = width*(j-1);
        for(i = 0 ; i < width ; i++)
        {
            if(channel == 3 && colors == 3)
            {
                data[3 * (pt+ i)]    = *(pBuff++);
                data[3 * (pt + i)+1] = *(pBuff++);
                data[3 * (pt + i)+2] = *(pBuff++);
            }
            else if(channel == 1 && colors == 3)
            {
                r = *(pBuff++);
                g = *(pBuff++);
                b = *(pBuff++);
                data[pt] = (r+g+b)/3;
                *maxVal = (*maxVal)>data[pt]?(*maxVal):data[pt];
                pt++;
            }
            else if(channel == 3 && colors == 1)
            {
                data[3 * (pt+ i)]    = *(pBuff);
                data[3 * (pt + i)+1] = *(pBuff);
                data[3 * (pt + i)+2] = *(pBuff++);
            }
            else
            {
                data[pt+ i] = *(pBuff++);
            }
        }

        pBuff = pBuff+jump;
    }
    delete bmpheader;
    delete [] buffer;
    buffer = NULL;
    delete bmih;
    return data;
}

char *IOImage::decmp_bmp(char *pBuff, const int w, const int h, int *maxVal)
{
    unsigned char firstByte = 0, secondByte = 0;
    unsigned int x = 0, y = 0, i = 0, state = 0;
    char *data = new char[w*h];
    memset(data, 0, w*h);
    *maxVal = 0;
    while(state == 0)
    {
        firstByte = *pBuff;
        pBuff++;
        if( firstByte != 0)
        {
            secondByte = *pBuff;
            for(i = 0; i < firstByte; i++)
            {
                data[y*w+x] = secondByte;
                *maxVal = (*maxVal < secondByte)?secondByte:(*maxVal);
                x++;
            }
            pBuff++;
        }
        else
        {
            firstByte = *pBuff;            // store next byte
            pBuff++;
            switch (firstByte)
            {
            case 0:
            {
                x = 0;
                y++;                       // next row
                break;
            }
            case 1:
            {
                state = 1;
                break;                // end of bitmap
            }
            case 2:
            {
                x += (int)(*pBuff);    // read byte and add value to x value
                pBuff++;
                y += (int)(*pBuff);    // read byte and add value to y value
                pBuff++;
                break;
            }
            default:
            {
                for(i = 0; i < firstByte; i++)
                {
                    secondByte = *pBuff;
                    data[y*w+x] = secondByte;
                    *maxVal = (*maxVal < secondByte)?secondByte:(*maxVal);
                    pBuff++;
                    x++;
                }
                if((firstByte & 0x01) != 0)   // if the run doesn't end on a word boundary,
                    pBuff++;                 // advance the pointer
            }
            } //end-switch
        } //end-if-else
    } //end-while
    return data;
}

char *IOImage::read_jpg(const char *srcfn, int &w, int &h,
                        int *MaxVal, const int channel)
{
    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;
    JSAMPARRAY buffer;
    int row_stride;
    FILE *infile;
    /** physical row width in output buffer **/

    if ((infile = fopen(srcfn, "rb")) == NULL)
    {
        fprintf(stderr, "can't open %s\n", srcfn);
        return 0;
    }

    /** Step 1: allocate and initialize JPEG decompression object **/
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    if (setjmp(jerr.setjmp_buffer))
    {
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return 0;
    }

    jpeg_create_decompress(&cinfo);
    /** Step 2: specify data source (eg, a file) **/
    jpeg_stdio_src(&cinfo, infile);
    /** Step 3: read file parameters with jpeg_read_header() **/
    jpeg_read_header(&cinfo, TRUE);

    /** Step 4: set parameters for decompression **/

    /** Step 5: Start decompressor **/
    (void) jpeg_start_decompress(&cinfo);
    /** JSAMPLEs per row in output buffer **/
    row_stride = cinfo.output_width * cinfo.output_components;
    w = cinfo.output_width;
    h = cinfo.output_height;

    int ch = cinfo.output_components;
    unsigned char *body = new unsigned char[ w*h*ch];
    unsigned char *data = new unsigned char[ w*h*channel];

    /** Make a one-row-high sample array that will go away when done with image **/
    buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    /** Step 6: while (scan lines remain to be read) **/
    int k, row_number;
    *MaxVal = 0;
    float pix_val;
    int   ich, _ich;
    unsigned char *row;

    while (cinfo.output_scanline < cinfo.output_height)
    {
        jpeg_read_scanlines(&cinfo, buffer, 1);
        row_number = cinfo.output_scanline-1;
        row = body+row_number*row_stride;
        for(k = 0; k < row_stride; k+=ch )
        {
            for(ich = 0; ich < ch; ich++)
            {
                _ich = ch - ich - 1;
                row[k+ich] = buffer[0][k+_ich];
                if(row[k+ich] > *MaxVal)
                {
                    *MaxVal = row[k+ich];
                }
            }
        }
    }
    int rowloc, rowloc_dt, i, j;

    if(ch == channel)
    {
        memcpy(data, body, sizeof(unsigned char)*w*h*ch);
    }
    else if (ch == 1 && channel == 3)
    {
        for(i = 0; i < h; i++)
        {
            for(j = 0; j < w; j+=ch)
            {
                rowloc_dt = (i*w+j)*channel;
                rowloc    = i*w + j;
                data[rowloc_dt] = body[rowloc];
                data[rowloc_dt+1] = body[rowloc];
                data[rowloc_dt+2] = body[rowloc];
            }
        }
    }
    else if(ch == 3 && channel == 1)
    {
        unsigned char tmpval = 0;
        for(i = 0; i < h; i++)
        {
            for(j = 0; j < w; j++)
            {
                rowloc    = (i*w+j)*ch;
                rowloc_dt = i*w+j;
                pix_val   = (body[rowloc]+body[rowloc+1]+body[rowloc+2])/3.0f;
                tmpval    = (unsigned char)round(pix_val);
                if(tmpval > *MaxVal)
                {
                    *MaxVal = tmpval;
                }
                data[rowloc_dt] = tmpval;
            }
        }
    }

    /** Step 7: Finish decompression **/
    (void) jpeg_finish_decompress(&cinfo);

    /** Step 8: Release JPEG decompression object **/
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    delete [] body;
    return (char*)data;
}

char *IOImage::read_ppm(const char *srcfn, int &w, int &h,
                        int *maxVal, const int channel)
{
    assert(srcfn);
    char red, green, blue, val;
    int  count = 0, num = 0, color[3], xi, x, y;
    unsigned int rowloc;
    unsigned char *data = NULL;
    vector<int> digits;
    char header[256];
    ifstream inStrm;
    bool  Binary = false;
    float gray;

    inStrm.open(srcfn, fstream::in|ifstream::binary);

    if(!inStrm.is_open())
    {
        cerr << "File '" << srcfn << "' cannot be opened.\n";
        return NULL;
    }
    inStrm.getline(header, 256);
    if(!strncmp(header, "P6", 2))
    {
        Binary = true;
    }
    else if (!strncmp(header, "P3", 2))
    {
        Binary = false;
    }
    else
    {
        cerr<<"'"<<header<<"'"<<" is not valid magic number for PGM format!\n";
        exit(0);
    }

    do
    {
        if(header[0] != '#')
        {
            VString::trimAfter(header, '#');
            count = VString::parse_words(header, digits);
            num  += count;
        }
        if(num < 3)
        {
            inStrm.getline(header, 256);
        }
    }
    while(num < 3 && !inStrm.eof());

    w = digits[0];
    h = digits[1];
    if(channel == 3)
        *maxVal = digits[2];
    else
        *maxVal = 0;
    digits.clear();

    data = new unsigned char [channel*w*h];

    if(Binary)
    {
        if(channel == 1)
        {
            for(y = 0; y < h; y++)
            {
                rowloc = y*w;
                for (x = 0; x < w; x++)
                {
                    inStrm.get(red);
                    inStrm.get(green);
                    inStrm.get(blue);
                    color[0] = (red&127) + (red&128);
                    color[1] = (green&127) + (green&128);
                    color[2] = (blue&127) + (blue&128);
                    gray     = (color[0]+color[1]+color[2]+0.0f)/3.0f;
                    val      = (unsigned char)floor(gray);
                    *maxVal  = *maxVal< val?val:*maxVal;
                    data[rowloc+x] = val;
                }
            }
        }
        else if(channel == 3)
        {
            xi = 0;
            for(y = 0; y < h; y++)
            {
                rowloc = y*w*channel;
                for (x = 0, xi = 0; x < w; x++, xi+=3)
                {
                    inStrm.get(red);
                    inStrm.get(green);
                    inStrm.get(blue);
                    color[0] = (red&127) + (red&128);
                    color[1] = (green&127) + (green&128);
                    color[2] = (blue&127) + (blue&128);
                    data[rowloc+xi+2] = color[0];
                    data[rowloc+xi+1] = color[1];
                    data[rowloc+xi]   = color[2];
                }
            }
        }
        inStrm.close();
    }///end if (Binary)
    else
    {
        inStrm.close();
        inStrm.open(srcfn, fstream::in);
        inStrm.getline(header, 256);
        do
        {
            if(header[0] != '#')
            {
                VString::trimAfter(header, '#');
                count = VString::parse_words(header, digits);
                num  += count;
            }
            if(num < 3)
            {
                inStrm.getline(header, 256);
            }
        }
        while(num < 3 && !inStrm.eof());

        digits.clear();
        if(channel == 1)
        {
            for (xi = 0; xi < w*h; xi++)
            {
                inStrm>>color[0];
                inStrm>>color[1];
                inStrm>>color[2];
                gray     = (color[0] + color[1] + color[2])/3.0f;
                val      = (unsigned char)round(gray);
                *maxVal  = *maxVal< val?val:*maxVal;
                data[xi]  = val;
            }
        }
        else
        {
            xi = 0;
            for(y = 0; y < h; y++)
            {
                rowloc = y*w*3;
                for (x = 0, xi = 0; x < w; x++, xi+=3)
                {
                    inStrm>>color[0];
                    inStrm>>color[1];
                    inStrm>>color[2];
                    data[rowloc+xi+2] = color[0];
                    data[rowloc+xi+1] = color[1];
                    data[rowloc+xi]   = color[2];
                }
            }
        }
        inStrm.close();
    }
    return (char*)data;
}


int IOImage::write_pgm(const char *fn, const int w, const int h, const unsigned char *data,
                       const unsigned char maxval, const char* comment_string, const int channel)
{
    FILE *file = fopen(fn, "wb");
    unsigned char *tmpdata = new unsigned char[w*h];
    int pt = 0;
    float val;
    unsigned int uval, mval = 0;
    int i, j;

    if(channel == 1)
    {
        for(i = 0; i < h; i++)
        {
            for(j = 0; j < w; j++)
            {
                tmpdata[pt] = data[pt];
                mval = mval < tmpdata[pt]?tmpdata[pt]:mval;
                pt++;
            }
        }
    }
    else if(channel == 3)
    {
        for(int i = 0; i < h; i++)
        {
            for(int j = 0; j < w; j++)
            {
                val = data[3*(i*w+j)];
                val = data[3*(i*w+j)+1] + val;
                val = data[3*(i*w+j)+2] + val;
                uval = (unsigned int)floor(val/3.0);
                mval = uval > mval?uval:mval;
                tmpdata[pt] = (unsigned char)uval;
                pt++;
            }
        }
    }

    if (file== NULL)
    {
        printf("File %s open failed in write pgm\n",fn);
        return(-1);
    }
    fprintf(file,"P5\n");

    if (comment_string != NULL)
        fprintf(file,"# %s\n", comment_string);
    fprintf(file,"%d %d\n", w, h);

    fprintf(file, "%d\n", mval);
    fwrite((void*)tmpdata, sizeof(unsigned char), (w)*(h), file);
    fclose(file);
    delete [] tmpdata;
    return 0;
}

int IOImage::write_pgm(const char *fn, int w, int h, const float *data, const unsigned char maxval,
                       const char *comment_string, const int channel)
{
    FILE *file= fopen(fn, "wb");
    unsigned char *tmpdata = new unsigned char[w*h];
    int pt = 0;
    float val;
    unsigned int uval, mval = 0;
    int i, j;

    if(channel == 1)
    {
        for(i = 0; i < h; i++)
        {
            for(j = 0; j < w; j++)
            {
                tmpdata[pt] = (unsigned char)round(data[pt]);
                mval = mval < tmpdata[pt]?tmpdata[pt]:mval;
                pt++;
            }
        }
    }
    else if(channel == 3)
    {
        for(int i = 0; i < h; i++)
        {
            for(int j = 0; j < w; j++)
            {
                val = data[3*(i*w+j)];
                val = data[3*(i*w+j)+1] + val;
                val = data[3*(i*w+j)+2] + val;
                uval = (unsigned int)floor(val/3.0);
                mval = uval > mval?uval:mval;
                tmpdata[pt] = (unsigned char)uval;
                pt++;
            }
        }
    }

    if (file== NULL)
    {
        printf("File %s open failed in write pgm\n",fn);
        return(-1);
    }
    fprintf(file, "P5\n");

    if (comment_string != NULL)
        fprintf(file,"# %s\n", comment_string);
    fprintf(file,"%d %d\n", w, h);

    fprintf(file, "%d\n", mval);
    fwrite((void*)tmpdata, sizeof(unsigned char), (w)*(h), file);
    fclose(file);
    delete [] tmpdata;
    return 0;
}

void IOImage::write_bmp(const char *dstfn,const int w,const int h,const unsigned char *body,const int channel)
{
    int image_data_size = w * h;
    image_data_size *= 4;

    int hexIdent[4];
    int hexWidth[4];
    int hexHeight[4];
    int hexFileSize[4];

    VMath::dec2hex(w, hexWidth);
    VMath::dec2hex(h, hexHeight);
    VMath::dec2hex(image_data_size+54, hexFileSize);
    VMath::dec2hex(image_data_size, hexIdent);

    FILE * maskFile  = fopen(dstfn , "w+b");

    char headerArray[54];
    headerArray[0] =(char)0x42 ;
    headerArray[1] =(char)0x4D ;
    headerArray[2] =(char)hexFileSize[0] ;
    headerArray[3] =(char)hexFileSize[1] ;
    headerArray[4] =(char)hexFileSize[2] ;
    headerArray[5] =(char)hexFileSize[3] ;
    headerArray[6] = (char)0x0;
    headerArray[7] = (char)0x0;
    headerArray[8] = (char)0x0;
    headerArray[9] = (char)0x0;
    headerArray[10] = (char)0x36;
    headerArray[11] = (char)0x0;
    headerArray[12] = (char)0x0;
    headerArray[13] = (char)0x0;
    headerArray[14] = (char)0x28;
    headerArray[15] = (char)0x0;
    headerArray[16] = (char)0x0;
    headerArray[17] = (char)0x0;
    headerArray[18] = (char)hexWidth[0];
    headerArray[19] = (char)hexWidth[1];
    headerArray[20] = (char)hexWidth[2];
    headerArray[21] = (char)hexWidth[3];
    headerArray[22] = (char)hexHeight[0];
    headerArray[23] = (char)hexHeight[1];
    headerArray[24] = (char)hexHeight[2];
    headerArray[25] = (char)hexHeight[3];
    headerArray[26] = (char)0x01;
    headerArray[27] = (char)0x0;
    headerArray[28] = (char)0x20;
    headerArray[29] = (char)0x0;
    headerArray[30] = (char)0x0;
    headerArray[31] = (char)0x0;
    headerArray[32] = (char)0x0;
    headerArray[33] = (char)0x0;
    headerArray[34] = (char)hexIdent[0];
    headerArray[35] = (char)hexIdent[1];
    headerArray[36] = (char)hexIdent[2];
    headerArray[37] = (char)hexIdent[3];
    headerArray[38] = (char)0xC4;
    headerArray[39] = (char)0x0E;
    headerArray[40] = (char)0x0;
    headerArray[41] = (char)0x0;
    headerArray[42] = (char)0xC4;
    headerArray[43] = (char)0x0E;
    headerArray[44] = (char)0x0;
    headerArray[45] = (char)0x0;
    headerArray[46] = (char)0x0;
    headerArray[47] = (char)0x0;
    headerArray[48] = (char)0x0;
    headerArray[49] = (char)0x0;
    headerArray[50] = (char)0x0;
    headerArray[51] = (char)0x0;
    headerArray[52] = (char)0x0;
    headerArray[53] = (char)0x0;

    fwrite(headerArray, sizeof(char), 54, maskFile);
    fclose(maskFile);
    maskFile  = fopen( dstfn, "a+b");

    unsigned char* data = new unsigned char[image_data_size];

    int index=0;
    int indexM;

    for(int m=0; m<h; m++)
    {
        for(int n=0; n<w; n++)
        {
            index   = channel*(m*w+n);
            indexM  = 4*((h-m-1)*w+n);

            if( channel == 3 )
            {
                data[indexM  ] = (int)round(body[index]);
                data[indexM+1] = (int)round(body[index+1]); //(int)round(body[index+1]);
                data[indexM+2] = (int)round(body[index+2]); //(int)round(body[index+2]);
                data[indexM+3] = 0;
            }
            else if( channel == 1 )
            {
                data[indexM  ] = (int)round(body[index]);
                data[indexM+1] = (int)round(body[index]);
                data[indexM+2] = (int)round(body[index]);
                data[indexM+3] = 0;
            }
        }
    }
    fwrite(data, sizeof(char), image_data_size, maskFile);
    fclose(maskFile);
    delete []data;
    data = NULL;
}

void IOImage::write_bmp(const char *dstfn,const int w, const int h, const float *body,const int channel)
{
    int image_data_size = 4*w * h;
    int hexWidth[4];
    int hexHeight[4];
    int hexFileSize[4];
    int hexIdent[4];

    VMath::dec2hex(w, hexWidth);
    VMath::dec2hex(h, hexHeight);
    VMath::dec2hex(image_data_size+54, hexFileSize);
    VMath::dec2hex(image_data_size, hexIdent);

    FILE * maskFile  = fopen(dstfn , "w+b");

    char headerArray[54];
    headerArray[0] =(char)0x42 ;
    headerArray[1] =(char)0x4D ;
    headerArray[2] =(char)hexFileSize[0] ;
    headerArray[3] =(char)hexFileSize[1] ;
    headerArray[4] =(char)hexFileSize[2] ;
    headerArray[5] =(char)hexFileSize[3] ;
    headerArray[6] = (char)0x0;
    headerArray[7] = (char)0x0;
    headerArray[8] = (char)0x0;
    headerArray[9] = (char)0x0;
    headerArray[10] = (char)0x36;
    headerArray[11] = (char)0x0;
    headerArray[12] = (char)0x0;
    headerArray[13] = (char)0x0;
    headerArray[14] = (char)0x28;
    headerArray[15] = (char)0x0;
    headerArray[16] = (char)0x0;
    headerArray[17] = (char)0x0;
    headerArray[18] = (char)hexWidth[0];
    headerArray[19] = (char)hexWidth[1];
    headerArray[20] = (char)hexWidth[2];
    headerArray[21] = (char)hexWidth[3];
    headerArray[22] = (char)hexHeight[0];
    headerArray[23] = (char)hexHeight[1];
    headerArray[24] = (char)hexHeight[2];
    headerArray[25] = (char)hexHeight[3];
    headerArray[26] = (char)0x01;
    headerArray[27] = (char)0x0;
    headerArray[28] = (char)0x20;
    headerArray[29] = (char)0x0;
    headerArray[30] = (char)0x0;
    headerArray[31] = (char)0x0;
    headerArray[32] = (char)0x0;
    headerArray[33] = (char)0x0;
    headerArray[34] = (char)hexIdent[0];
    headerArray[35] = (char)hexIdent[1];
    headerArray[36] = (char)hexIdent[2];
    headerArray[37] = (char)hexIdent[3];
    headerArray[38] = (char)0xC4;
    headerArray[39] = (char)0x0E;
    headerArray[40] = (char)0x0;
    headerArray[41] = (char)0x0;
    headerArray[42] = (char)0xC4;
    headerArray[43] = (char)0x0E;
    headerArray[44] = (char)0x0;
    headerArray[45] = (char)0x0;
    headerArray[46] = (char)0x0;
    headerArray[47] = (char)0x0;
    headerArray[48] = (char)0x0;
    headerArray[49] = (char)0x0;
    headerArray[50] = (char)0x0;
    headerArray[51] = (char)0x0;
    headerArray[52] = (char)0x0;
    headerArray[53] = (char)0x0;

    fwrite(headerArray, sizeof(char), 54, maskFile);
    fclose(maskFile);
    maskFile  = fopen( dstfn, "a+b");

    unsigned char* data = new unsigned char[image_data_size];

    int index=0;
    int indexM;
    for(int m = 0; m < h; m++)
    {
        for(int n=0; n<w; n++)
        {
            index   = channel*(m*w+n);
            indexM  = 4*((h-m-1)*w+n);

            if( channel == 3 )
            {
                data[indexM  ] = (int)round(body[index]);
                data[indexM+1] = (int)round(body[index+1]); //(int)round(body[index+1]);
                data[indexM+2] = (int)round(body[index+2]); //(int)round(body[index+2]);
                data[indexM+3] = 0;
            }
            else if( channel == 1 )
            {
                data[indexM  ] = (int)round(body[index]);
                data[indexM+1] = (int)round(body[index]);
                data[indexM+2] = (int)round(body[index]);
                data[indexM+3] = 0;
            }
        }
    }
    fwrite(data, sizeof(char), image_data_size, maskFile);
    fclose(maskFile);
    delete []data;
    data = NULL;
}

void IOImage::write_jpg(const char *srcfn, const unsigned char *data, const int w, const int h, const int ch, const int quality)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    assert(data);

    FILE * outfile;
    JSAMPROW row_pointer[1];
    /** pointer to JSAMPLE row[s] **/

    int row_stride;
    /** physical row width in image buffer **/
    /** Step 1: allocate and initialize JPEG compression object **/
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    unsigned char* body = new unsigned char[h*w*ch];
    int i,  j,  rowloc;
    row_stride = w*ch;
    if(ch == 3)
    {
        for(i = 0; i < h; i++)
        {
            rowloc = i*w*ch;
            for(j = 0; j < row_stride; j+=ch)
            {
                body[rowloc+j] = (unsigned char) data[rowloc+j+2];
                body[rowloc+j+1] = (unsigned char) data[rowloc+j+1];
                body[rowloc+j+2] = (unsigned char) data[rowloc+j];
            }
        }
    }
    else
    {
        for(i = 0; i < h; i++)
        {
            rowloc = i*w*ch;
            for(j = 0; j < row_stride; j+=ch)
            {
                body[rowloc+j] = data[rowloc+j];
            }
        }
    }

    /** Step 2: specify data destination (eg, a file) **/
    /** Note: steps 2 and 3 can be done in either order. **/

    if ((outfile = fopen(srcfn, "wb")) == NULL)
    {
        fprintf(stderr, "can't open %s\n", srcfn);
        exit(1);
    }

    jpeg_stdio_dest(&cinfo, outfile);
    /** Step 3: set parameters for compression **/

    cinfo.image_width = w;
    cinfo.image_height = h;
    cinfo.input_components = ch;

    if( ch == 3 )
        cinfo.in_color_space = JCS_RGB;
    else
        cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);

    //cout<<"bug 1.1\n";

    /** Now you can set any non-default parameters you wish to.
    * Here we just illustrate the use of quality (quantization table) scaling:
    **/

    /** limit to baseline-JPEG values **/
    jpeg_set_quality(&cinfo, quality, TRUE);

    /** Step 4: Start compressor */
    jpeg_start_compress(&cinfo, TRUE);

    /** Step 5: while (scan lines remain to be written) */
    /**           jpeg_write_scanlines(...); */

    /** Here we use the library's state variable cinfo.next_scanline as the
    * loop counter, so that we don't have to keep track ourselves.
    * To keep things simple, we pass one scanline per call; you can pass
    * more if you wish, though.*/

    row_stride = w * ch;
    /** JSAMPLEs per row in image_buffer */
    while (cinfo.next_scanline < cinfo.image_height)
    {
        row_pointer[0] = & body[cinfo.next_scanline * row_stride];
        (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    /** Step 6: Finish compression **/
    jpeg_finish_compress(&cinfo);
    fclose(outfile);

    /** Step 7: release JPEG compression object **/
    jpeg_destroy_compress(&cinfo);
    delete [] body;
}

void IOImage::write_jpg(const char *srcfn, const float *data, const int w, const int h, const int ch, const int quality)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    assert(data);

    FILE * outfile;
    JSAMPROW row_pointer[1];
    /** pointer to JSAMPLE row[s] **/
    int row_stride;
    /** physical row width in image buffer **/
    /** Step 1: allocate and initialize JPEG compression object **/
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    unsigned char* body = new unsigned char[h*w*ch];

    int i, j, rowloc;
    row_stride = w*ch;
    if(ch == 3)
    {
        for(i = 0; i < h; i++)
        {
            rowloc = i*w*ch;
            for(j = 0; j < row_stride; j+=ch)
            {
                body[rowloc+j] = (unsigned char) round(data[rowloc+j+2]);
                body[rowloc+j+1] = (unsigned char) round(data[rowloc+j+1]);
                body[rowloc+j+2] = (unsigned char)round( data[rowloc+j]);
            }
        }
    }
    else
    {
        for(i = 0; i < h; i++)
        {
            rowloc    = i*w;
            for(j = 0; j < row_stride; j++)
            {
                body[rowloc+j] = data[rowloc+j];
            }
        }
    }

    /** Step 2: specify data destination (eg, a file) **/
    /** Note: steps 2 and 3 can be done in either order. **/

    if ((outfile = fopen(srcfn, "wb")) == NULL)
    {
        fprintf(stderr, "can't open %s\n", srcfn);
        exit(1);
    }

    jpeg_stdio_dest(&cinfo, outfile);
    /** Step 3: set parameters for compression **/

    cinfo.image_width = w;
    cinfo.image_height = h;
    cinfo.input_components = ch;

    if( ch == 3)
        cinfo.in_color_space = JCS_RGB;
    else
        cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    /** Now you can set any non-default parameters you wish to.
    * Here we just illustrate the use of quality (quantization table) scaling: **/

    /** limit to baseline-JPEG values **/
    jpeg_set_quality(&cinfo, quality, TRUE);

    /** Step 4: Start compressor */
    jpeg_start_compress(&cinfo, TRUE);

    /** Step 5: while (scan lines remain to be written) */
    /**           jpeg_write_scanlines(...); */

    row_stride = w * ch;
    /** JSAMPLEs per row in image_buffer */

    while (cinfo.next_scanline < cinfo.image_height)
    {
        row_pointer[0] = & body[cinfo.next_scanline * row_stride];
        (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    /** Step 6: Finish compression **/
    jpeg_finish_compress(&cinfo);
    fclose(outfile);

    /** Step 7: release JPEG compression object **/
    jpeg_destroy_compress(&cinfo);

    //cout<<"i saved jpeg\n";

    delete [] body;
}


void IOImage::write_ppm(const char *srcfn, const unsigned char *data, const int w, const int h, const int ch)
{
    FILE *pFile;
    int  x, y, xi;
    unsigned int rowloc;
    assert(srcfn);
    assert(data);
    pFile = fopen(srcfn, "wb");
    if(pFile == NULL)
    {
        cout<<"File '"<<srcfn<<"' cannot open for write!\n";
        return;
    }
    fprintf(pFile, "P6\n%d %d\n255\n", w, h);
    unsigned char *tmpdata = NULL;
    tmpdata = new unsigned char[w*h*3];
    /**source is in gray level**/
    if(ch == 1)
    {
        for(y = 0; y < h; y++)
        {
            rowloc = y*w*3;
            for(x = 0, xi = 0; x < w; x++, xi+=3)
            {
                tmpdata[rowloc+xi] = data[rowloc+x];
                tmpdata[rowloc+xi+1] = data[rowloc+x];
                tmpdata[rowloc+xi+2] = data[rowloc+x];
            }
        }
        for(y = 0; y < h; y++)
        {
            fwrite(tmpdata+y*w*3, 1, w*3, pFile);
        }
        delete [] tmpdata;
    }
    else
    {
        for(y = 0; y < h; y++)
        {
            rowloc = y*w*3;
            for(x = 0, xi = 0; x < w; x++, xi+=3)
            {
                tmpdata[rowloc+xi]   = data[rowloc+xi+2];
                tmpdata[rowloc+xi+1] = data[rowloc+xi+1];
                tmpdata[rowloc+xi+2] = data[rowloc+xi];
            }
        }
        for(y = 0; y < h; y++)
        {
            fwrite(tmpdata+y*w*3, 1, w*3, pFile);
        }
        delete [] tmpdata;
    }

    fclose(pFile);
    return ;
}

void IOImage::write_ppm(const char *srcfn, const float *data,const int w, const int h, const int ch)
{
    FILE *pFile;
    int  x, y, xi;
    unsigned int rowloc;
    assert(srcfn);
    assert(data);

    pFile = fopen(srcfn, "wb");

    if(pFile==NULL)
    {
        cout<<"File '"<<srcfn<<"' cannot open for write!\n";
        return;
    }

    fprintf(pFile, "P6\n%d %d\n255\n", w, h);

    unsigned char *tmpdata = NULL;
    tmpdata = new unsigned char[w*h*3];
    /**source is in gray level**/
    if(ch == 1)
    {
        for(y = 0; y < h; y++)
        {
            rowloc = y*w*3;
            for(x = 0, xi = 0; x < w; x++, xi+=3)
            {
                tmpdata[rowloc+xi]   = (unsigned char)round(data[rowloc+x]);
                tmpdata[rowloc+xi+1] = (unsigned char)round(data[rowloc+x]);
                tmpdata[rowloc+xi+2] = (unsigned char)round(data[rowloc+x]);
            }
        }
        for(y = 0; y < h; y++)
        {
            fwrite(tmpdata+y*w*3, 1, w*3, pFile);
        }
        delete [] tmpdata;
    }
    else
    {
        for(y = 0; y < h; y++)
        {
            rowloc = y*w*3;
            for(x = 0, xi = 0; x < w; x++, xi+=3)
            {
                tmpdata[rowloc+xi] = (unsigned char)round(data[rowloc+xi+2]);
                tmpdata[rowloc+xi+1] = (unsigned char)round(data[rowloc+xi+1]);
                tmpdata[rowloc+xi+2] = (unsigned char)round(data[rowloc+xi]);
            }
        }
        for(y = 0; y < h; y++)
        {
            fwrite(tmpdata+y*w*3, 1, w*3, pFile);
        }
        delete [] tmpdata;
    }

    fclose(pFile);
    return ;
}

char *IOImage::read_png(const char *srcfn, int &w, int &h, int *maxVal, const int channel)
{
    char *data = NULL;
    int x, y, width, height, locrow, val = 0, sz = 0;
    int bit_depth = 0;
    png_byte color_type;

    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep *row_pointers = NULL;
    png_byte header[PNG_BYTES_TO_CHECK];

    /** open file and test for it being a png */
    FILE *fp = fopen(srcfn, "rb");
    if(!fp)
    {
        cerr<<"Error happens when reading '"<<srcfn<<"'\n";
        return NULL;
    }
    int si = fread(header, 1, PNG_BYTES_TO_CHECK, fp);
    if(si != PNG_BYTES_TO_CHECK)
    {
        cerr<<"Error happens when reading PNG image header'"<<srcfn<<"'\n";
        return NULL;
    }

    png_sig_cmp(header, 0, PNG_BYTES_TO_CHECK);
    /**
    if(!png_sig_cmp(header, (png_size_t)0, PNG_BYTES_TO_CHECK))
    {
          cerr<<"Error happens when reading '"<<srcfn<<"'\n";
          return NULL;
    }
    **/

    /** initialize stuff */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
    {
        cerr<<"[read_png] png_create_read_struct failed"<<endl;
        return NULL;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        cerr<<"[read_png] png_create_info_struct failed"<<endl;
        return NULL;
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[read_png] Error during init_io"<<endl;
        png_destroy_info_struct(png_ptr, &info_ptr);
        ///png_destroy_read_struct(&png_ptr, &info_ptr);
        return NULL;
    }

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    width      = png_get_image_width(png_ptr,  info_ptr);
    height     = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr,   info_ptr);
    bit_depth  = png_get_bit_depth(png_ptr,    info_ptr);

    if(width == 0 || height == 0)
    {
        cerr<<"[read_png] Error size of the image"<<"\tbit_depth"<<bit_depth<<endl;
        return NULL;
    }
    /**
    cout<<width<<"\t"<<height<<endl;
    cout<<"type: "<<(int)color_type<<endl;
    cout<<"depth: "<<(int)bit_depth<<endl;
    **/

    int num_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);


    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[read_png_file] Error during png_jmpbuf(...)"<<"\tpasses"<<num_of_passes<<endl;
        return NULL;
    }

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (y = 0; y < height; y++)
    {
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr, info_ptr));
    }

    png_read_image(png_ptr, row_pointers);
    w = width;
    h = height;

    if(channel == 1)
    {
        data = new char[width*height];
        sz = width*height;
    }
    else if(channel == 3)
    {
        data = new char[width*height*3];
        sz = width*height*3;
    }
    else
    {
        data = new char[width*height];
        sz = width*height;
    }
    memset(data, 0, sizeof(char)*sz);
    int step = 0;

    if(color_type == PNG_COLOR_TYPE_GRAY)
    {
        step = 1;
    }
    else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
    {
        step = 4;
    }
    else if(color_type == PNG_COLOR_TYPE_RGB)
    {
        step = 3;
    }
    else if (color_type == PNG_COLOR_TYPE_RGBA)
    {
        step = 4;
    }
    else
    {
        w = h = 0;
        cout<<"Unknown png format!\n";
        return NULL;
    }
    if(channel == 1)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = channel*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*step]);
                if(step == 4 && ptr[3] == 0)
                {
                    data[locrow+x] = 255;
                }
                else
                {
                    if(step >= 3)
                    {
                        val = floor((ptr[0]+ptr[1]+ptr[2])/3.0);
                    }
                    else
                    {
                        val = ptr[0];
                    }
                    data[locrow+x] = val;
                }
            }
        }
    }
    else if(channel == 3)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = 3*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*step]);
                if(step == 4 && ptr[3] == 0)
                {
                    data[locrow+3*x]   = 255;
                    data[locrow+3*x+1] = 255;
                    data[locrow+3*x+2] = 255;
                }
                else if(step == 3)
                {
                    data[locrow+3*x]   = ptr[2];  //Blue
                    data[locrow+3*x+1] = ptr[1];  //Green
                    data[locrow+3*x+2] = ptr[0];  //Red
                }
                else if (step == 2)
                {
                    if(ptr[1] == 0)
                    {
                        data[locrow+3*x]   = 255;
                        data[locrow+3*x+1] = 255;
                        data[locrow+3*x+2] = 255;
                    }
                    else
                    {
                        data[locrow+3*x]   = ptr[0];  //Blue
                        data[locrow+3*x+1] = ptr[0];  //Green
                        data[locrow+3*x+2] = ptr[0];  //Red
                    }
                }
                else if(step == 1)
                {
                    data[locrow+3*x]   = ptr[0];  //Blue
                    data[locrow+3*x+1] = ptr[0];  //Green
                    data[locrow+3*x+2] = ptr[0];  //Red
                }
            }
        }
    }
    ///if-else
    /** cleanup heap allocation */
    for (y = 0; y < height; y++)
    {
        free(row_pointers[y]);
    }
    free(row_pointers);
    png_destroy_info_struct(png_ptr, &info_ptr);
    fclose(fp);
    return data;
}


void IOImage::write_png(const char* srcfn, const int width, const int height,
                        const char *data, const int channel)
{
    int x = 0, y = 0, locrow = 0;
    int bit_depth = 8;
    int color_type = PNG_COLOR_TYPE_RGB;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep *row_pointers = NULL;

    FILE *fp = fopen(srcfn, "wb");
    if(!fp)
    {
        cerr<<"[write_png] Error during png_jmpbuf(...)"<<endl;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
    {
        cerr<<"[write_png_file] png_create_write_struct failed"<<endl;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        cerr<<"[write_png_file] png_create_info_struct failed"<<endl;
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during init_io"<<endl;
    }

    png_init_io(png_ptr, fp);

    if(setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during writing header"<<endl;
    }

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (y = 0; y < height; y++)
    {
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr, info_ptr));
    }

    if(channel == 1)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = channel*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*3]);
                ptr[0] = ptr[1] = ptr[2] = data[locrow+x];
            }
        }
    }
    else if(channel == 3)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = 3*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*3]);
                ptr[2] = data[locrow+3*x];    //Blue
                ptr[1] = data[locrow+3*x+1];  //Green
                ptr[0] = data[locrow+3*x+2];  //Red
            }
        }
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    /** write bytes **/
    if(setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during writing bytes"<<endl;
    }

    png_write_image(png_ptr, row_pointers);

    /** end write */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during end of write"<<endl;
    }

    png_write_end(png_ptr, NULL);

    /** cleanup heap allocation */
    for (y = 0; y < height; y++)
    {
        free(row_pointers[y]);
    }
    free(row_pointers);
    png_destroy_info_struct(png_ptr, &info_ptr);
    fclose(fp);
}


void IOImage::write_png(const char* srcfn, const int width, const int height,
                        const float *data, const int channel)
{
    int x = 0, y = 0, locrow = 0;
    int bit_depth = 8;
    int color_type = PNG_COLOR_TYPE_RGB;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep *row_pointers = NULL;

    FILE *fp = fopen(srcfn, "wb");
    if(!fp)
    {
        cerr<<"[write_png] Error during png_jmpbuf(...)"<<endl;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
    {
        cerr<<"[write_png_file] png_create_write_struct failed"<<endl;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        cerr<<"[write_png_file] png_create_info_struct failed"<<endl;
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during init_io"<<endl;
    }

    png_init_io(png_ptr, fp);

    if(setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during writing header"<<endl;
    }

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (y = 0; y < height; y++)
    {
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr, info_ptr));
    }

    if(channel == 1)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = channel*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*3]);
                ptr[0] = ptr[1] = ptr[2] = (int)floor(data[locrow+x]);
            }
        }
    }
    else if(channel == 3)
    {
        for(y = 0; y < height; y++)
        {
            png_byte* row = row_pointers[y];
            locrow = 3*y*width;
            for (x = 0; x < width; x++)
            {
                png_byte* ptr = &(row[x*3]);
                ptr[2] = (int)floor(data[locrow+3*x]);    //Blue
                ptr[1] = (int)floor(data[locrow+3*x+1]);  //Green
                ptr[0] = (int)floor(data[locrow+3*x+2]);  //Red
            }
        }
    }

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    /** write bytes **/
    if(setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during writing bytes"<<endl;
    }

    png_write_image(png_ptr, row_pointers);

    /** end write */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        cerr<<"[write_png_file] Error during end of write"<<endl;
    }

    png_write_end(png_ptr, NULL);

    /** cleanup heap allocation */
    for (y = 0; y < height; y++)
    {
        free(row_pointers[y]);
    }
    free(row_pointers);
    png_destroy_info_struct(png_ptr, &info_ptr);
    fclose(fp);
}

void IOImage::testbmp()
{
    const char *fn = "";
    const char *fnw = "";
    int height = 0, width = 0;
    int maxg = 0,ch = 1;
    unsigned char *pixels = NULL;
    /**/
    //pixels = (unsigned char*)IOImage::read_bmp(fn, width, height, &maxg, 3);
    pixels = (unsigned char*)IOImage::read_jpg(fn, width, height, &maxg, ch);
    //IOImage::write_bmp(fnw, width, height, pixels, 3);
    cout<<width<<"\t"<<height<<endl;
    IOImage::write_jpg(fnw, pixels, width, height, 1, 60);
    //IOImage::pgm_pgmwrite(fnw, width, height, pixels, maxg, "", 3);
}

void IOImage::test()
{
    const char *pngimg = "/home/wlzhao/datasets/vgg/bt/bt1.bmp";
    const char *dstimg = "/home/wlzhao/sunflower4.pgm";
    /**/
    int w = 0, h = 0, val = 0;
    unsigned char *data = NULL;

    data = (unsigned char*)IOImage::read_bmp(pngimg, w, h, &val, 1);
    if(data == NULL)
        return;

    IOImage::write_pgm(dstimg, w, h, data, 255, "", 1);
    if(data != NULL)
        delete [] data;
    /**/
    //delete my;
    return ;
}
