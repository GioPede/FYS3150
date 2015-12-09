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
    int numberOfUnitCells = 20;
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26*1.5);   // measured in angstroms
    double dt = UnitConverter::timeFromSI(1e-14);                        // measured in seconds

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


    for(double T = UnitConverter::temperatureFromSI(50.0); T > UnitConverter::temperatureFromSI(45.0); T -= UnitConverter::temperatureFromSI(10.0)){
        cout<< "New Temperature: "<<UnitConverter::temperatureToSI(T)<<endl;
        for(double vx = 2.0*pow(10,-26); vx < pow(10,-24); vx += pow(10,-26)){
            // CREATE AND INITIALIZE THE LATTICE:
            System system;
            system.createFCCLattice(numberOfUnitCells,UnitConverter::lengthFromSI(pow(vx,1.0/3.0)), T);
            double sigma = UnitConverter::lengthFromAngstroms(3.405);
            double epsilon = UnitConverter::temperatureFromSI(119.8);
            system.setPotential(new LennardJones(sigma, epsilon));
            system.potential()->rCut = 2.5*sigma;
            system.potential()->rShell = 0.3*sigma;
            system.setIntegrator(new VelocityVerlet());
            system.setThermostat(new Thermostat(UnitConverter::timeFromSI(1e-14),T));
            system.removeTotalMomentum();
            system.thermostat()->setOn();

            StatisticsSampler statisticsSampler;

            int numberOfTimesteps = 3000;
            double dT=0.05;
            statisticsSampler.sample(system);
            int count=0;
            while(count<1000){
                system.step(dt);
                statisticsSampler.sample(system);
               if((statisticsSampler.temperature() > T-dT)||(statisticsSampler.temperature()<T+dT)) count++;
               else count=0;
            }
            if (vx == 1.0*pow(10,-28))cout<< "Thermalized!"<<endl;
            double V = vx;
            double T_real = 0.0;
            double P = 0.0;
            for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
                system.step(dt);
                statisticsSampler.sample(system);
                P += statisticsSampler.pressure();
                T_real += statisticsSampler.temperature();
            }
            cout << UnitConverter::temperatureToSI(T_real/numberOfTimesteps) << ",     " << V << ",     "
                 << UnitConverter::pressureToSI( P/numberOfTimesteps) << endl;
        }
    }

    return 0;
}
