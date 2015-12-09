#ifndef SYSTEM_H
#define SYSTEM_H
#include "atom.h"
#include "math/vec3.h"
#include <vector>
#include "celllist.h"

class Potential; class Integrator; class Thermostat;
using std::vector;

class System
{
private:
    vec3 m_systemSize;
    vector<Atom*> m_atoms;
    Potential* m_potential = NULL;
    Integrator* m_integrator = NULL;
    Thermostat* m_thermostat = NULL;
    double m_time = 0;
    int m_steps = 0;
    double m_temperature = 1.0;
    CellList m_cellList;

    void validateNeighborList();
public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature);
    void applyPeriodicBoundaryConditions();
    void removeTotalMomentum();
    void calculateForces();
    void calculateForcesNoBoundaries();
    void step(double dt);
    void stepNoBoundaries(double dt);
    // Setters and getters
    vector<Atom *>& atoms() { return m_atoms; } // Returns a reference to the std::vector of atom pointers
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double time() { return m_time; }
    void setTime(double time) { m_time = time; }
    double temperature() { return m_temperature; }
    void setTemperature(double temperature) { m_temperature = temperature; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    Thermostat *thermostat() { return m_thermostat; }
    void setThermostat(Thermostat *thermostat) { m_thermostat = thermostat; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    double volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
    CellList &cellList();
};
#endif
