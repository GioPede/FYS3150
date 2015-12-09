#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <omp.h>
#include <cmath>
using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    omp_set_num_threads(8);  // initialize OpenMP

    //DEFAULT PARAMETERS:
    int numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(1500.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26);   // measured in angstroms
    double dt = UnitConverter::timeFromSI(1e-15);                        // measured in seconds

    // USER DEFINED PARAMETERS:
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
    if(numberOfArguments > 4) dt = UnitConverter::timeFromSI(atof(argumentList[4]));

    // PRINT DIMENSION SCALING FACTORS:
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    for(double T =  UnitConverter::temperatureFromSI(100.0); T < UnitConverter::temperatureFromSI(3000.0); T +=  UnitConverter::temperatureFromSI(100.0)){
    // CREATE AND INITIALIZE THE LATTICE:
    System system;
    bool NoBoundaries = 0;      // flag to tigger boundary conditions
    system.createFCCLattice(numberOfUnitCells, latticeConstant, T);
    system.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405),
                                         UnitConverter::temperatureFromSI(119.8)));
    system.setIntegrator(new VelocityVerlet());
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    /*
    IO movie, data;
    movie.open("movie.xyz", 1);
    data.open("data.txt", 0);
    movie.saveState(&system);
    cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << endl;
    for(int timestep=0; timestep<20000; timestep++) {
        if(NoBoundaries) system.stepNoBoundaries(dt);
        else system.step(dt);
        statisticsSampler.sample(system);
        data.saveStatistcalData(&system, &statisticsSampler);
        if( !(timestep % 100) ) {
            cout << system.steps() << "      "
                 << UnitConverter::timeToSI(system.time()) << "      "
                 << UnitConverter::temperatureToSI(statisticsSampler.temperature()) << "      "
                 << statisticsSampler.kineticEnergy() << "      "
                 << statisticsSampler.potentialEnergy() << "      "
                 << statisticsSampler.totalEnergy()<< "      "
                 << statisticsSampler.diffusionConstant()<< "      "
                 << statisticsSampler.pressure()<< "      "
                 << statisticsSampler.momentum()
                 << endl;
            movie.saveState(&system);
        }
    }

    movie.close(1);
    data.close(0);

    */
    double N = 3000;
    for(int timestep=0; timestep<N; timestep++) {
        system.step(dt);
    }
    double D = 0.0;
    double D2 = 0.0;
    for(int timestep=0; timestep<N; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        D+=statisticsSampler.diffusionConstant();
        D2+=statisticsSampler.diffusionConstant()*statisticsSampler.diffusionConstant();
    }
    D2 = sqrt(D2/N - D/N*D/N);
    printf("%.f,    %.15f,    %.15f\n", UnitConverter::temperatureToSI(statisticsSampler.temperature()),
                                        UnitConverter::diffusionToSI(D/N),
                                        UnitConverter::diffusionToSI(D2));
    }

    return 0;
}
