#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "system.h"

class Thermostat
{
    private:
        double m_tau;
        double m_temperature;
        bool is_on;
    public:
        Thermostat(double tau, double T);
        void changeVelocities(System *system, double dt);

        void setTau(double tau) { m_tau = tau; }
        void setTemperature(double T) { m_temperature = T; }
        void setOn() { is_on = 1; }
        void setOff() { is_on = 0; }
};

#endif // THERMOSTAT_H
