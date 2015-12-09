#ifndef IO_H
#define IO_H
#include <fstream>
class System;
class StatisticsSampler;
using std::ofstream;

class IO
{
private:
    ofstream movie;
    ofstream data;
public:
    IO();
    ~IO();

    void saveState(System *system);
    void saveStatistcalData(System *system, StatisticsSampler *statisticsSampler);
    void open(const char *filename, bool movie_flag);
    void close(bool movie_flag);

};
#endif
