#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdio>
#include <cmath>
#include <list>

#include "vstring.h"

using namespace std;

VString::VString()
{}

bool VString::endWith(const char str[], const char end[])
{
    int lenend = strlen(end);
    int lenstr = strlen(str);
    if(lenend > lenstr)
        return false;

    int start = lenstr - lenend, i, j = 0;
    for(i = start; i < lenstr; i++,j++)
    {
        if(str[i]!= end[j])
            return false;
    }
    return true;
}

bool VString::startWith(const char str[], const char end[])
{
    int lenend  = strlen(end);
    int lenstr  = strlen(str), i;
    if(lenend > lenstr)
    {
        return false;
    }

    if(lenend == 0)
    {
        return true;
    }

    for(i = 0; i<lenend; i++)
    {
        if(str[i] != end[i])
            return false;
    }
    return true;
}

bool VString::existDir(const char *dir)
{
    if(dir == NULL || !strcmp(dir, ""))
      return false;

    DIR* mydir = opendir(dir);
    if(mydir == NULL)
    {
        return false;
    }
    else
    {
        closedir(mydir);
        return true;
    }
}


bool VString::validatePath(const char *path)
{
    unsigned int len = strlen(path);

    if(len == 0)
    return false;

    char *tmpdir = new char[len+1];
    bool  rslt   = VString::parseDIR(tmpdir, path);
    if(!rslt)
    {
        delete [] tmpdir;
        return false;
    }

    rslt = VString::existDir(tmpdir);


    if(!rslt)
    {
        delete [] tmpdir;
        return false;
    }
    else{
        char *tmpfn = new char[len+1];
        rslt = VString::parseFName(tmpfn, path);

        if(!rslt)
        {
            return false;
        }else{
            if(tmpfn[0] >= 'a' && tmpfn[0] <= 'z')
            {
                delete [] tmpfn;
                return true;
            }else if(tmpfn[0] >= 'A' && tmpfn[0] <= 'Z')
            {
                delete [] tmpfn;
                return true;
            }else if(tmpfn[0] >= '0' && tmpfn[0] <= '9')
            {
                delete [] tmpfn;
                return true;
            }else{
                delete [] tmpfn;
                return false;
            }
        }
    }
    return false;
}

bool VString::existFile(const char *fn)
{
    FILE* myfp = fopen(fn, "r");
    if(myfp == NULL)
    {
        return false;
    }
    else
    {
        fclose(myfp);
        return true;
    }
}

void VString::toUpper(char str[])
{
    int len = strlen(str), i;
    for(i = 0; i <len; i++)
    {
        str[i] = toupper(str[i]);
    }
}

void VString::toLower(char str[])
{
    int len = strlen(str), i;
    for(i = 0; i < len; i++)
    {
        str[i] = tolower(str[i]);
    }
    return ;
}

bool VString::parseDIR(char dirname[], const char path[])
{
    assert(dirname);
    assert(path);
    strcpy(dirname, "");

    if(strlen(path) == 0)
    {
        return false;
    }
    int  sidx  = 0;

    if(VString::startWith(path, "./"))
    {
        strcpy(dirname, "./");
        return true;
    }

    if(VString::startWith(path, ".\\"))
    {
        strcpy(dirname, ".\\");
        return true;
    }

    int  eidx  = VString::lastindexof(path, '/');
    int  _eidx = VString::lastindexof(path, '\\');

    if(eidx == -1 && _eidx == -1)
    {
        strcpy(dirname, "./");
        return true;
    }

    eidx = eidx>_eidx?eidx:_eidx;

    if(eidx == 0)
    {
        sprintf(dirname, "/");
    }

    if(sidx < eidx)
    {
        VString::subbtw(dirname, path, sidx, eidx);
        VString::replase(dirname, '\\', '/');
        return true;
    }
    else
    {
        return false;
    }
}

bool VString::parseFName(char fname[], const char path[])
{
    assert(fname);
    assert(path);
    strcpy(fname, "");

    if(strlen(path) == 0)
    {
        return false;
    }

    int  sidx  = VString::lastindexof(path, '\\');
    int  _sidx = VString::lastindexof(path, '/');
    int  eidx  = VString::lastindexof(path, '.');

    if(eidx == -1)
    {
        eidx = strlen(path) - 1;
    }
    else
    {
        if(path[0] == '.' && eidx == 0)
        {
            eidx = strlen(path) - 1;
        }else
        {
            eidx = eidx - 1;
        }
    }

    sidx = sidx>_sidx?sidx:_sidx;
    ///cout<<sidx<<"\t"<<eidx<<endl;
    if(sidx < eidx)
    {
        sidx = sidx + 1;
        VString::subbtw(fname, path, sidx, eidx);
        ///cout<<fname<<endl;
        return true;
    }
    else
    {
        return false;
    }
}

bool VString::validFName(const char fname[])
{
    int len = strlen(fname);
    if(len <= 0)
    {
        return false;
    }
    else if((fname[0] >= 'a' && fname[0] <= 'z')
            ||(fname[0] >= 'A' && fname[0] <= 'Z'))
    {
        return true;
    }else{
        return false;
    }
}

bool VString::parsePath(char fname[], const char path[])
{
    assert(fname);
    assert(path);
    strcpy(fname, "");

    if(strlen(path) == 0)
    {
        return false;
    }

    if(!strcmp(path, "./") ||!strcmp(path, ".\\"))
    {
        return false;
    }

    char *name = new char[strlen(path)+1];
    strcpy(name, "");

    int eidx = VString::lastindexof(path, '.');
    ///int vidx = eidx + 1;

    VString::parseFName(name, path);

    if(VString::validFName(name))
    {
        if(eidx <= 0)
        {
            strcpy(fname, path);
            VString::replase(fname, '\\', '/');
        }else{
            strncpy(fname, path, eidx);
            VString::replase(fname, '\\', '/');
            fname[eidx] = '\0';
        }
        delete [] name;
        return true;
    }else{
        delete [] name;
        return false;
    }
}

void VString::replase(char srcstr[], const char ch, const char nwch)
{
    unsigned int len = 0, i;
    if((len = strlen(srcstr)) == 0)
    {
        return ;
    }
    for(i = 0; i < len; i++)
    {
        if(srcstr[i] == ch)
        {
            srcstr[i] = nwch;
        }
    }
    return ;
}

void VString::subbtw(char dest[], const char src[],
                     const int start, const int end)
{
    int i = 0, r_end = 0, k = 0;
    int flen = 0;
    flen = strlen(src);
    assert(dest);

    if(start == flen)
    {
        return ;
    }
    r_end = end >= flen?(flen-1):end;
    for(i = start; i <= r_end; i++)
    {
        dest[k] = src[i];
        k++;
    }

    dest[k] = '\0';
    return ;
}

int VString::firstindexof(const char *src, const char ch)
{
    int i = 0, len = strlen(src);
    for(i = 0; i < len; i++)
    {
        if(src[i] == ch)
            return i;
    }
    return -1;
}

int VString::lastindexof(const char *src, char ch)
{
    int i, end = strlen(src) - 1;
    for(i = end; i >= 0; i--)
    {
        if(src[i] == ch)
            return i;
    }
    return -1;
}

void VString::split_twin(char *first, char *second, const char *src, const char ch)
{
    int i, j, len = strlen(src)-1;

    assert(first);
    assert(second);

    for(i = 0; i < len; i++)
    {
        if(src[i] != ch)
        {
            first[i] = src[i];
        }
        else
        {
            break;
        }
    }
    first[i] = '\0';
    i = i+1;

    for(j = i; j <= len; j++)
    {
        second[j-i] =  src[j];
    }
    second[j-i] = '\0';
    return ;
}

void VString::parsePair(const char src_str[], string &fst_str,
                     string &secd_str, const char splitter)
{
    int len, i, j;
    len =  strlen(src_str);

    for(i = 0; i < len; i++)
    {
        if(src_str[i] == splitter)
        {
            break;
        }
    }
    for(j=0; j < i; j++)
    {
        fst_str += src_str[j];
    }

    for(j = i+1; j < len; j++)
    {
        secd_str += src_str[j];
    }
    return ;
}

int VString::parse_words(const char *str, vector<int> &nums)
{
    char ch[2] = {' ', '\t'};
    int  i, len = strlen(str), count = 0, idx = 0;
    char *word  = new char[len+1];
    strcpy(word, "");
    bool FOUL  = false;
    for(i = 0; i < len; i++)
    {
        if(str[i] == ch[0] || str[i] == ch[1])
        {
            if(idx != 0)
            {
                if(!FOUL)
                {
                    word[idx]   = '\0';
                    idx         = 0;
                    nums.push_back(atoi(word));
                }
                else
                {
                    FOUL = false;
                    idx  = 0;
                }
            }
            if(i > 0)
            {
                if(str[i-1] != ch[0] && str[i-1] != ch[1])
                {
                    count++;
                }
            }
        }
        else
        {
            if(str[i] < '0' || str[i] > '9')
            {
                FOUL = true;
            }
            word[idx] = str[i];
            idx++;
        }
    }
    if(idx != 0)
    {
        word[idx]   = '\0';
        idx         = 0;
        if(!FOUL)
        {
            nums.push_back(atoi(word));
        }
    }

    delete [] word;
    return nums.size();
}

void VString::splitby(const char *line, const char splitter, vector<string> &words)
{
    assert(line != NULL);

    int   i = 0, len = strlen(line);
    char ch = ' ';
    string word = "";

    for(i = 0; i < len; i++)
    {
        ch = line[i];
        if(ch != splitter)
        {
            word = word + ch;
        }
        else
        {
            words.push_back(word);
            word = "";
        }
    }
    words.push_back(word);
    return ;
}

void VString::parseLine(const char *line, const char spliter, set<int> &regTab)
{
    assert(line != NULL);

    int len = strlen(line), i, size = 0, val;
    char word[30], ch =  ' ';
    bool PARSED = true;
    strcpy(word, "");

    for(i = 0; i < len; i++)
    {
        ch = line[i];

        if((ch >='0' && ch <= '9'))
        {
            PARSED     = false;
            word[size] = ch;
            size++;
        }
        else if(ch == '+'||ch=='-')
        {
            word[size] = ch;
            size++;
        }
        else if(ch == '.')
        {
            if(!PARSED)
            {
                word[size] = ch;
                size++;
            }
        }
        else
        {
            if(!PARSED)
            {
                word[size] ='\0';
                val = atoi(word);
                regTab.insert(val);

                PARSED = true;
                strcpy(word,"");
                size = 0;
            }
        }
    }
    word[size]='\0';
    if(!PARSED&& strcmp(word,""))
    {
        val = atoi(word);
        regTab.insert(val);
    }
}

void VString::parseLine(const char *line, const char spliter, vector<float> &parsedVals)
{
    assert(line != NULL);
    int len = strlen(line), i, size = 0;
    char word[30], ch;
    bool PARSED = true;
    strcpy(word, "");

    for(i = 0; i < len; i++)
    {
        ch = line[i];

        if(ch >='0' && ch <= '9')
        {
            PARSED     = false;
            word[size] = ch;
            size++;
        }
        else if(ch == '+'||ch=='-')
        {
            word[size] = ch;
            size++;
        }
        else if(ch == '.')
        {
            if(!PARSED)
            {
                word[size] = ch;
                size++;
            }
        }
        else
        {
            if(!PARSED)
            {
                word[size] = '\0';
                parsedVals.push_back(atof(word));
                PARSED = true;
                strcpy(word, "");
                size   = 0;
            }
        }
    }
    word[size] = '\0';
    if(!PARSED && strcmp(word, ""))
    {
        parsedVals.push_back(atof(word));
    }
}

set<unsigned int> VString::parseLine(const char *line, const char spliter)
{
    assert(line != NULL);
    set<unsigned int> mylist;

    int len = strlen(line), i, size = 0;
    char word[30], ch;
    unsigned int val = 0;
    bool PARSED = true;

    for(i = 0; i < len; i++)
    {
        ch = line[i];

        if(ch >= '0' && ch <= '9')
        {
            PARSED = false;
            word[size]=ch;
            size++;
        }
        else if(ch == '+'||ch == '-')
        {
            word[size]=ch;
            size++;
        }
        else if(ch == '.')
        {
            if(!PARSED)
            {
                word[size]=ch;
                size++;
            }
        }
        else
        {
            if(!PARSED)
            {
                word[size]='\0';
                val=atoi(word);
                mylist.insert(val);

                PARSED = true;
                strcpy(word,"");
                size   = 0;
            }
        }
    }
    word[size] = '\0';
    if(!PARSED && strcmp(word, ""))
    {
        val = atoi(word);
        mylist.insert(val);
    }
    return mylist;
}

void VString::trimEnd(char *line, const char ch)
{
    int len = strlen(line), i, j;

    if(len <= 0)
    {
        return ;
    }

    for(i = len-1; i > 0; i--)
    {
        if(line[i] != ch)
        {
            break;
        }
    }
    for(j = i+1; j < len; j++)
    {
        line[j] = '\0';
    }
    return ;
}

void VString::trimHead(char *line, const char ch)
{
    int i = 0, j = 0, len = strlen(line);

    if(len == 0)
        return ;

    for(i = 0; i < len; i++)
    {
        if(line[i] != ch)
        {
            break;
        }
    }
    for(j = 0; line[i] != '\0'; j++,i++)
    {
        line[j] = line[i];
    }
    return ;
}

void VString::trimStops(char *line)
{
    if(line == NULL)
    return ;

    if(strlen(line) == 0)
    return ;

    VString::trimHead(line, ' ');
    VString::trimEnd(line,  ' ');
    VString::trimEnd(line,  '\n');
    VString::trimEnd(line,  '\r');
    return ;
}

void VString::trimAfter(char *line, const char ch)
{
    int i, len = strlen(line);
    for(i = 0; i < len; i++)
    {
        if(line[i] == ch)
        {
            line[i] ='\0';
            break;
        }
    }
    return ;
}

void VString::time2Str(char *tmStr, const int seconds)
{
    int realtm = seconds;
    int tm_hour = realtm/3600;
    realtm = seconds%3600;
    int tm_min = realtm/60;
    int tm_sec = realtm%60;
    //char *tmstr = new char[12];
    char hour[3];
    char min[3];
    char sec[3];

    if(tm_hour>9)
    {
        sprintf(hour,"%d",tm_hour);
    }
    else if(tm_hour>0)
    {
        sprintf(hour,"0%d",tm_hour);
    }
    else
    {
        sprintf(hour,"0");
    }
    if(tm_min>9)
    {
        sprintf(min,"%d",tm_min);
    }
    else if(tm_min >0)
    {
        sprintf(min,"0%d",tm_min);
    }
    else
    {
        sprintf(min,"00");
    }

    if(tm_sec>9)
    {
        sprintf(sec,"%d",tm_sec);
    }
    else if(tm_sec > 0)
    {
        sprintf(sec,"0%d",tm_sec);
    }
    else
    {
        sprintf(sec,"00");
    }

    sprintf(tmStr,"%s:%s:%s",hour,min,sec);
}

const char *VString::time2Str(const int seconds)
{
    int realtm = seconds;
    int tm_hour = realtm/3600;
    realtm = seconds%3600;
    int tm_min = realtm/60;
    int tm_sec = realtm%60;
    char *tmstr = new char[12];
    char hour[3];
    char min[3];
    char sec[3];


    if(tm_hour > 9)
    {
        sprintf(hour,"%d",tm_hour);
    }
    else if(tm_hour > 0)
    {
        sprintf(hour, "0%d", tm_hour);
    }
    else
    {
        sprintf(hour, "0");
    }
    if(tm_min > 9)
    {
        sprintf(min, "%d", tm_min);
    }
    else if(tm_min > 0)
    {
        sprintf(min, "0%d", tm_min);
    }
    else
    {
        sprintf(min, "00");
    }

    if(tm_sec > 9)
    {
        sprintf(sec, "%d", tm_sec);
    }
    else if(tm_sec > 0)
    {
        sprintf(sec, "0%d", tm_sec);
    }
    else
    {
        sprintf(sec, "00");
    }

    sprintf(tmstr,"%s:%s:%s", hour, min, sec);
    return tmstr;
}

void VString::test()
{
    char dir[] = "2435sdfsrwese";

    if(VString::startWith(dir, "243"))
    {
        cout<<"It contains\n";
    }
    ///bool t = VString::validatePath("/home/1");
    ///cout<<t<<endl;
    bool t;
    char dest[64] = "";

    const char *hi[10] = {"/", "sdf/aa", "/mnt/a", "\\var/s1r", "s/er\1aa",
                         "/home/ss1.1", "/xe\\x", "12aass.1", "./a1", ".\\a"};
    for(int i = 0; i < 10; i++)
    {
        ///t = VString::parseFName(dest, hi[i]);
        ///t = VString::parseDIR(dest, hi[i]);
        ///t = VString::validatePath(hi[i]);
        t = VString::parsePath(dest, hi[i]);
        cout<<hi[i]<<"\t"<<t<<"\t"<<dest<<endl;
        //cout<<i<<endl;
    }
}

VString::~VString()
{

}
