#ifndef BOARD_H
#define BOARD_H

#include <cstring>

class Board
{
    private:
        unsigned int width, height;
        unsigned char *pix;
    public:
        Board(unsigned int w, unsigned h)
        {
            pix = new unsigned char[w*h];
            width = w; height = h;
            memset(pix, 0, w*h);
        }
        unsigned char getPixel(unsigned int x,  unsigned int y)
        {
             if(x < width && y < height)
                return pix[y*width+x];
             else
                return 1;
        }
        void setPixel(unsigned int x, unsigned int y, unsigned char val)
        {
             pix[y*width+x] = val;
        }
        unsigned int getWidth()
        {
            return width;
        }
        unsigned int getHeight()
        {
            return height;
        }
        virtual ~Board()
        {
           delete [] pix;
           pix = NULL;
        }
};

#endif // BOARD_H
