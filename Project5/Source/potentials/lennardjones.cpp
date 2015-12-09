#include "lennardjones.h"
#include <math.h>
#include <omp.h>
#include <algorithm>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

// SOME OLD STUFF
double LennardJones::calculateMutualPotential(double norm){
    double a6 = (m_sigma/norm);
    a6 = a6*a6*a6;
    a6=a6*a6;
    return 4*m_epsilon*(a6*a6-a6);

}

// STILL OLD STUFF
double LennardJones::calculateMutualForce(double dist){
    double a6 = (m_sigma/dist);
    a6 = a6*a6*a6;
    a6=a6*a6;
    return 4*m_epsilon*(a6 / dist * 6. -
                        a6*a6 / dist * 12.);

}

// MAIN LOOP OF THE PROGRAM. COMPUTES AND SETS ALL THE FORCES BETWEEN ATOMS
void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0;
    m_pressure = 0;

    double L = system->systemSize().x();
    double sigma6 = pow(m_sigma, 6);
    double energyAtRCut = 4*m_epsilon*pow(rCut,-6)*sigma6*(pow(rCut,-6)*sigma6 - 1.0);
    for(Atom *atom1 : system->atoms())
    {
        // set variables for the first atom, we don't want too many function calls
        double x1 = atom1->position.x();
        double y1 = atom1->position.y();
        double z1 = atom1->position.z();

        for(Atom *atom2 : atom1->neighbors) {
            double x2 = atom2->position.x();
            double x2x1 = x2-x1;
            double y2 = atom2->position.y();
            double y2y1 = y2-y1;
            double z2 = atom2->position.z();
            double z2z1 = z2-z1;
            if(x2x1 < -0.5*L) x2x1 += L;
            else if(x2x1 > 0.5*L) x2x1 -= L;

            if(y2y1 < -0.5*L) y2y1 += L;
            else if(y2y1 > 0.5*L) y2y1 -= L;

            if(z2z1 < -0.5*L) z2z1 += L;
            else if(z2z1 > 0.5*L) z2z1 -= L;
            double dr2 = x2x1*x2x1 + y2y1*y2y1 + z2z1*z2z1;
            if(dr2 > rCut*rCut) continue;

            double oneOverDr2 = 1.0/dr2;
            double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;
            double F = -24*m_epsilon*oneOverDr6*sigma6*(2.0*oneOverDr6*sigma6 - 1.0)*oneOverDr2;

            atom1->force[0] += F*x2x1;
            atom1->force[1] += F*y2y1;
            atom1->force[2] += F*z2z1;

            atom2->force[0] -= F*x2x1;
            atom2->force[1] -= F*y2y1;
            atom2->force[2] -= F*z2z1;

            m_potentialEnergy += 4*m_epsilon*oneOverDr6*sigma6*(oneOverDr6*sigma6 - 1.0) - energyAtRCut;
            m_pressure += atom1->force[0]*x2x1  +
                          atom1->force[1]*y2y1  +
                          atom1->force[2]*z2z1  ;
        }
    }
}

// FUNCTION TO COMPUTE THE FREE EXPANSION, NOT VERY RELIABLE
void LennardJones::calculateForcesNoBoundaries(System *system)
{
    m_potentialEnergy = 0;
    unsigned int i, j, N;
    double F, norm,x1,x2,y1,y2,z1,z2,x2x1,y2y1,z2z1,fx,fy,fz,ForceX,ForceY,ForceZ,cut_dist;
    Atom* M1;
    Atom* M2;
    N = system->atoms().size();
    cut_dist = m_sigma*5;
    for(i = 0; i < N; i++)
    {
        for(j = i+1 ; j < N; j++){

            M1 = system->atoms()[i];
            M2 = system->atoms()[j];

            x1 = M1->position.x();
            x2 = M2->position.x();
            x2x1 = x2-x1;
            fx = fabs(x2x1);
            if(fx > cut_dist) continue;

            y1 = M1->position.y();
            y2 = M2->position.y();
            y2y1 = y2-y1;
            fy = fabs(y2y1);
            if(fy > cut_dist) continue;

            z1 = M1->position.z();
            z2 = M2->position.z();
            z2z1 = z2-z1;
            fz = fabs(z2z1);
            if(fz > cut_dist) continue;


            norm = sqrt(fx*fx + fy*fy + fz*fz);

            m_potentialEnergy  += calculateMutualPotential(norm);
            F = calculateMutualForce(norm);

            ForceX = (F/norm*x2x1);
            ForceY = (F/norm*y2y1);
            ForceZ = (F/norm*z2z1);

            M1->force.setX(M1->force.x() + ForceX);
            M1->force.setY(M1->force.y() + ForceY);
            M1->force.setZ(M1->force.z() + ForceZ);
            M2->force.setX(M2->force.x() - ForceX);
            M2->force.setY(M2->force.y() - ForceY);
            M2->force.setZ(M2->force.z() - ForceZ);
        }
    }
}

