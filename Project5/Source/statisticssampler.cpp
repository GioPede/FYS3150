#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include "unitconverter.h"
#include <cmath>

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleMomentum(system);
    sampleDiffusionConstant(system);
    samplePressure(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    double v2 = 0;
    for(Atom *atom : system.atoms()) {
        v2 += atom->velocity.lengthSquared();
    }
    Atom* A = system.atoms()[1];
    m_kineticEnergy = 0.5 * A->mass() * v2;
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential()->potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    m_temperature = 2.0/3.0 * m_kineticEnergy / system.atoms().size();
    system.setTemperature(m_temperature);
}

void StatisticsSampler::sampleDensity(System &system)
{

}


void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    m_diffusionConstant = 0;
    double r2 = 0;
    for(Atom *atom : system.atoms()){
        r2 += (atom->position.x()-atom->initialPosition.x())*
              (atom->position.x()-atom->initialPosition.x());
        r2 += (atom->position.y()-atom->initialPosition.y())*
              (atom->position.y()-atom->initialPosition.y());
        r2 += (atom->position.z()-atom->initialPosition.z())*
              (atom->position.z()-atom->initialPosition.z());
    }
    m_diffusionConstant = r2 / system.atoms().size() / 6 / system.time();
}

void StatisticsSampler::sampleMomentum(System &system)
{
    vec3 p;
    m_momentum = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms())
        p += atom->velocity;
    p.set(p.x()/system.atoms().size(), p.y()/system.atoms().size(), p.z()/system.atoms().size());
    m_momentum = p.length();
}

void StatisticsSampler::samplePressure(System &system)
{
    //double N = system.atoms().size()/pow(system.systemSize().x(), 3);
    //m_pressure = N*8.314472*m_temperature*119.8/(pow(system.systemSize().x() * 1e-10, 3) - N * b_pressure)
    //           - a_pressure * N/pow(system.systemSize().x() * 1e-10, 6)*N;
    //m_pressure = N * 8.314472*m_temperature*119.8/pow(system.systemSize().x() * 1e-10, 3);
    double rho = system.atoms().size()/pow(system.systemSize().x(), 3);
    m_pressure = system.potential()->pressure() /( 3.0 * pow(system.systemSize().x(), 3))
               + rho * m_temperature;
}



