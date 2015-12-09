#include "velocityverlet.h"
#include "potentials/potential.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    if(firstStep) {
        system->potential()->buildNeighborList(*system);
        system->calculateForces();
        firstStep = false;
    }

    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / 2 / atom->mass();
        atom->position += atom->velocity*dt;
    }
    system->applyPeriodicBoundaryConditions();
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / 2 / atom->mass();
    }

}

void VelocityVerlet::integrateNoBoundaries(System *system, double dt)
{
    system->calculateForcesNoBoundaries();
    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / 2 / atom->mass();
        atom->position += atom->velocity*dt;
    }
    system->calculateForcesNoBoundaries();
    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / 2 / atom->mass();
    }
}
