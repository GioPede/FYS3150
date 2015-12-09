#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "potential.h"

class LennardJones : public Potential
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;
public:
    double calculateMutualPotential(double norm);
    double calculateMutualForce(double dist);
    LennardJones(double sigma, double epsilon);
    virtual void calculateForces(System *system);
    virtual void calculateForcesNoBoundaries(System *system);
    void calculateForcesOld(System *system);
};
#endif
