#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include "statisticssampler.h"
#include <cstdlib>
using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close(0);
    close(1);
}

void IO::open(const char *filename, bool movie_flag) {
    if (movie_flag){
        if(movie.is_open()) {
            std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
            exit(1);
        }
        movie.open(filename);
    }

    else{
        if(data.is_open()) {
            std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
            exit(1);
        }
        data.open(filename);
    }
}

void IO::close(bool movie_flag) {
    if (movie_flag){
        if(movie.is_open())
            movie.close();
    }

    else{
        if(data.is_open())
            data.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
    if(movie.is_open()) {
        movie << system->atoms().size() << endl;
        movie << "The is an optional comment line that can be empty." << endl;
        for(Atom *atom : system->atoms()) {
            movie << "Ar " << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;
        }
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveStatistcalData(System *system, StatisticsSampler *statisticsSampler)
{
    if(data.is_open()) {
        data << system->steps() << "      "
             << system->time() << "      "
             << statisticsSampler->temperature() << "      "
             << statisticsSampler->kineticEnergy() << "      "
             << statisticsSampler->potentialEnergy() << "      "
             << statisticsSampler->totalEnergy()<< "      "
             << statisticsSampler->diffusionConstant()<< "      "
             << statisticsSampler->pressure()<< "      "
             << statisticsSampler->momentum()
             << endl;
    }
}
