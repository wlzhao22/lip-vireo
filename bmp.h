#ifndef _BMP_H_
#define _BMP_H_

// magic number "BM"
#define			BITMAP_ID			('B' + ('M'<<8))

// header byte type for RLE
const int		RLE_COMMAND		=	0;
const int		RLE_ENDOFLINE	=	0;
const int		RLE_ENDOFBITMAP	=	1;
const int		RLE_DELTA		=	2;

const int		BI_OS2			=	-1;

const int		BI_RGB			=	0;
const int		BI_RLE8			=	1;
const int		BI_RLE4			=	2;
const int		BI_BITFIELDS	=	3;

//#pragma warning( disable : 4103 ) //commented out by wanlei

// --------------------------------------------
// tagBITMAPFILEHEADER - bitmap file header.
// --------------------------------------------

#pragma pack(2)

typedef struct tagBITMAPFILEHEADER
{
        unsigned short	bfType;                 // magic number "BM"
        unsigned int	bfSize;                 // file size
        unsigned short	bfReserved1;            // reserved
        unsigned short	bfReserved2;            // reserved
        unsigned int	bfOffBits;                      // offset to bitmap data

} BITMAPFILEHEADER, *PBITMAPFILEHEADER;


#pragma pack(4)

// --------------------------------------------
// tagBITMAPCOREHEADER - bitmap core header.
// --------------------------------------------

typedef struct tagBITMAPCOREHEADER		// bmch
{
	unsigned int	bcSize;				// size of the structure
	unsigned int	bcWidth;			// image width
	unsigned int	bcHeight;			// image height
	unsigned int	bcPlanes;			// must be equal to 1
	unsigned int	bcBitCount;			// number of bits per pixel

} BITMAPCOREHEADER, *PBITMAPCOREHEADER;


// --------------------------------------------
// BITMAPFILEHEADER - bitmap info header.
// --------------------------------------------

typedef struct tagBITMAPINFOHEADER
{
	unsigned int	biSize;				// size of the structure
	unsigned int	biWidth;			// image width
	unsigned int	biHeight;			// image height
	unsigned short	biPlanes;			// must be equal to 1
	unsigned short	biBitCount;			// number of bits per pixel
	unsigned int	biCompression;		// compression type
	unsigned int	biSizeImage;		// size of data bitmap
	int		biXPelsPerMeter;	// number of pixels per meter on the X axis
	int		biYPelsPerMeter;	// number of pixels per meter on the Y axis
	unsigned int	biClrUsed;			// number of colors used
	unsigned int	biClrImportant;		// number of important colors

} BITMAPINFOHEADER, *PBITMAPINFOHEADER;

#endif

