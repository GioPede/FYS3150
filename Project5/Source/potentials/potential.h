#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy = 0;
    double m_pressure = 0;
public:
    double rCut = 0;
    double rShell = 0;
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesNoBoundaries(System *system) = 0;
    double potentialEnergy();
    double pressure();
    void setPotentialEnergy(double potentialEnergy);
    void buildNeighborList(System &system);
};
#endif
