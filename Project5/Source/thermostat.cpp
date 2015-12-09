#include "thermostat.h"
#include <cmath>

Thermostat::Thermostat(double tau, double T)
{
    m_tau = tau;
    m_temperature = T;
}


void Thermostat::changeVelocities(System *system, double dt)
{
    if(is_on){
        double gamma = sqrt(1 + dt / m_tau * (m_temperature / system->temperature() - 1));
        for(Atom *atom : system->atoms()){
            atom->velocity *= gamma;
        }
    }
}
