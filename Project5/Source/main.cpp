#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "thermostat.h"
#include <iostream>
#include <omp.h>
#include <cmath>
#include <QElapsedTimer>
using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    //omp_set_num_threads(8);  // initialize OpenMP

    //DEFAULT PARAMETERS:
    int numberOfUnitCells = 10;
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
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

    // CREATE AND INITIALIZE THE LATTICE:
    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    double sigma = UnitConverter::lengthFromAngstroms(3.405);
    double epsilon = UnitConverter::temperatureFromSI(119.8);
    system.setPotential(new LennardJones(sigma, epsilon));
    system.potential()->rCut = 2.5*sigma;
    system.potential()->rShell = 0.3*sigma;
    system.setIntegrator(new VelocityVerlet());
    system.setThermostat(new Thermostat(UnitConverter::timeFromSI(1e-13),
                                        UnitConverter::temperatureFromSI(300.0)));
    system.removeTotalMomentum();
    system.thermostat()->setOff();

    int numberOfTimesteps = 2000;

    StatisticsSampler statisticsSampler;
    IO movie, data;
    movie.open("movie.xyz", 1);
    data.open("data.txt", 0);
    movie.saveState(&system);
    cout << "Timestep   Time  Temperature   KineticEnergy  PotentialEnergy"
         << "TotalEnergy   DiffusionK   Pressure    Momentum" << endl;


    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
            system.step(dt);
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
                     << UnitConverter::pressureToSI(statisticsSampler.pressure())<< "      "
                     << statisticsSampler.momentum()
                     << endl;
                movie.saveState(&system);
            }
        }

        movie.close(1);
        data.close(0);
    return 0;
}
