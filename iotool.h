#ifndef IOTOOL_H
#define IOTOOL_H

#include <fstream>
#include <vector>
#include "keypoint.h"

#ifndef ACCESS_MODE
#define ACCESS_MODE 0777
#endif

using namespace std;

class IOTool
{
public:
    static bool   IS_Dir(const char *dir);
    static bool   IS_FILE(const char *fn);
    static bool   getKeypoint(const char *keyfn, vector<KeyPoint *> &kplst);
    static float *load_pcaMat(const char *pcaMatFn, vector<float> &means,
                           vector<float> &variances, unsigned int &row, unsigned int &col);
    static float  *loadMat(const char *pcaMatFn, unsigned int &row, unsigned int &col);

    static bool creatDIR(const char* dirname);
};

#endif
