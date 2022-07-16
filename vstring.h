#ifndef VSTRING_H
#define VSTRING_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <set>

using namespace std;

class VString
{
private:
    static void subbtw(char dest[], const char src[],
                         const int start, const int end);
    static void replase(char src[], const char ch, const char nwch);

public:
    VString();
    static void splitby(const char *line, const char spliter, vector<string> &words);
    static void parsePair(const char src_str[], string &fst_str,
                          string &secd_str, const char splitter);

    static bool startWith(const char *str, const char *pttn);
    static bool endWith(const char str[],  const char end[]);
    static void trimEnd(char *line,   const char ch);
    static void trimHead(char *line,  const char ch);
    static void trimAfter(char *line, const char ch);
    static void trimStops(char *line);

    static int  firstindexof(const char *src, const char ch);
    static int  lastindexof(const char *src,  const char ch);
    static int  count_words(const char *str);
    static int  parse_words(const char *str, vector<int> &nums);

    static void split_twin(char *first, char *second, const char *src, const char ch);

    static void parseLine(const char *line, const char spliter, set<int> &regTab);
    static set<unsigned int> parseLine(const char *line, const char spliter);
    static void parseLine(const char *line, const char spliter, vector<float> &vals);

    static bool existFile(const char *fn);
    static bool existDir(const char  *dir);
    static bool validatePath(const char *path);

    static bool parseFName(char fname[], const char path[]);
    static bool parseDIR(char dirname[], const char path[]);
    static bool parsePath(char fname[],  const char path[]);
    static bool validFName(const char fname[]);
    static void toUpper(char str[]);
    static void toLower(char str[]);

    static const char *time2Str(const int seconds);
    static void time2Str(char *tmStr, const int seconds);

    static void test();
    virtual ~VString();

};

#endif
