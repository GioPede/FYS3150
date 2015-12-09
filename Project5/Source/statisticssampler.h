#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H

class System;
class StatisticsSampler
{
private:
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_momentum = 0;
    double m_diffusionConstant = 0;
    double m_pressure = 0;
    double a_pressure = 0.1355;
    double b_pressure = 0.00003201;


public:
    StatisticsSampler();
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleMomentum(System &system);
    void sampleDensity(System &system);
    void samplePressure(System &system);
    void sampleDiffusionConstant(System &system);
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double momentum() { return m_momentum; }
    double diffusionConstant(){ return m_diffusionConstant; }
    double pressure(){return m_pressure;}
};
#endif
