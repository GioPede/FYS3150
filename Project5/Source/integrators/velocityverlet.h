#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H
#include "integrators/integrator.h"

class VelocityVerlet : public Integrator
{
public:
    bool firstStep = true;
    VelocityVerlet() { }
    virtual void integrate(System *system, double dt) override;
    virtual void integrateNoBoundaries(System *system, double dt) override;
};
#endif
