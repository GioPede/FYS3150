#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "thermostat.h"
#include "math/random.h"
#include <algorithm>

CellList &System::cellList()
{
    return m_cellList;
}

System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for(Atom *atom : m_atoms) {
        if(atom->position.x() > m_systemSize.x()){
            atom->position.setX(atom->position.x() - m_systemSize.x());
            atom->initialPosition.setX(atom->initialPosition.x() - m_systemSize.x());
            atom->positionOnNeighborlistBuild.setX(atom->positionOnNeighborlistBuild.x() - m_systemSize.x());
        }
        if(atom->position.y() > m_systemSize.y()){
            atom->position.setY(atom->position.y() - m_systemSize.y());
            atom->initialPosition.setY(atom->initialPosition.y() - m_systemSize.y());
            atom->positionOnNeighborlistBuild.setY(atom->positionOnNeighborlistBuild.y() - m_systemSize.y());
        }
        if(atom->position.z() > m_systemSize.z()){
            atom->position.setZ(atom->position.z() - m_systemSize.z());
            atom->initialPosition.setZ(atom->initialPosition.z() - m_systemSize.z());
            atom->positionOnNeighborlistBuild.setZ(atom->positionOnNeighborlistBuild.z() - m_systemSize.z());
        }
        if(atom->position.x() < 0){
            atom->position.setX(atom->position.x() + m_systemSize.x());
            atom->initialPosition.setX(atom->initialPosition.x() + m_systemSize.x());
            atom->positionOnNeighborlistBuild.setX(atom->positionOnNeighborlistBuild.x() + m_systemSize.x());
        }
        if(atom->position.y() < 0){
            atom->position.setY(atom->position.y() + m_systemSize.y());
            atom->initialPosition.setY(atom->initialPosition.y() + m_systemSize.y());
            atom->positionOnNeighborlistBuild.setY(atom->positionOnNeighborlistBuild.y() + m_systemSize.y());
        }
        if(atom->position.z() < 0){
            atom->position.setZ(atom->position.z() + m_systemSize.z());
            atom->initialPosition.setZ(atom->initialPosition.z() + m_systemSize.z());
            atom->positionOnNeighborlistBuild.setZ(atom->positionOnNeighborlistBuild.z() + m_systemSize.z());
        }
    }
}

void System::removeTotalMomentum() {
    vec3 p;
    for(Atom *atom : m_atoms)
        p += atom->velocity;
    for(Atom *atom : m_atoms){
        atom->velocity.setX(atom->velocity.x() - p.x()/m_atoms.size());
        atom->velocity.setY(atom->velocity.y() - p.y()/m_atoms.size());
        atom->velocity.setZ(atom->velocity.z() - p.z()/m_atoms.size());
    }
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    for(int i = 0; i < numberOfUnitCellsEachDimension; i++)
        for(int j = 0; j < numberOfUnitCellsEachDimension; j++)
            for(int k = 0; k < numberOfUnitCellsEachDimension; k++)
                for(int l = 0; l < 4; l++)
                {
                    double x,y,z;
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    if (l == 0){
                        x = i * latticeConstant;
                        y = j * latticeConstant;
                        z = k * latticeConstant;
                    }
                    else if (l == 1){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant;
                    }
                    else if (l == 2){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant;
                        z = k * latticeConstant+latticeConstant/2.;
                    }
                    else if (l == 3){
                        x = i * latticeConstant;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant+latticeConstant/2.;
                    }

                    atom->position.set(x,y,z);
                    atom->initialPosition.set(x,y,z);
                    atom->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom);

                }
    setSystemSize(vec3(latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension));
}

void System::validateNeighborList() {
    for(Atom *atom : m_atoms) {
        double dr = (atom->position - atom->positionOnNeighborlistBuild).length();
        if(dr>0.5*m_potential->rShell) {
            m_potential->buildNeighborList(*this);
            break;
        }
    }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    validateNeighborList();
    m_potential->calculateForces(this);
}

void System::calculateForcesNoBoundaries() {
    resetForcesOnAllAtoms();
    m_potential->calculateForcesNoBoundaries(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_thermostat->changeVelocities(this, dt);
    m_steps++;
    m_time += dt;
}

void System::stepNoBoundaries(double dt) {
    m_integrator->integrateNoBoundaries(this, dt);
    m_steps++;
    m_time += dt;
}

