#include "potential.h"
#include <iostream>
using namespace std;
Potential::Potential()
{

}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}

void Potential::buildNeighborList(System &system)
{
    for(Atom *atom : system.atoms()) {
        atom->neighbors.clear();
        atom->positionOnNeighborlistBuild = atom->position;
    }
    system.cellList().cutoffDistance = rShell + rCut;
    system.cellList().build(&system);
    CellList &cellList = system.cellList();
    double L = system.systemSize().x();
    double neighborListRadiusSquared = (rShell+rCut)*(rShell+rCut);
    for(int cx = 0; cx<cellList.numberOfCellsPerDimension; cx++) {
        for(int cy = 0; cy<cellList.numberOfCellsPerDimension; cy++) {
            for(int cz = 0; cz<cellList.numberOfCellsPerDimension; cz++) {
                vector<Atom*> &cell1 = cellList.cell(cx,cy,cz);
                for(int dx = -1; dx<=1; dx++) {
                    for(int dy = -1; dy<=1; dy++) {
                        for(int dz = -1; dz<=1; dz++) {
                            int i = (cx + dx + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int j = (cy + dy + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int k = (cz + dz + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            vector<Atom*> &cell2 = cellList.cell(i,j,k);

                            for(Atom *atomi : cell1) {
                                for(Atom *atomj : cell2) {
                                    if(atomi<=atomj) continue;
                                    double dx = atomi->position.x() - atomj->position.x();
                                    double dy = atomi->position.y() - atomj->position.y();
                                    double dz = atomi->position.z() - atomj->position.z();
                                    if(dx < -0.5*L) dx += L;
                                    else if(dx > 0.5*L) dx -= L;

                                    if(dy < -0.5*L) dy += L;
                                    else if(dy > 0.5*L) dy -= L;

                                    if(dz < -0.5*L) dz += L;
                                    else if(dz > 0.5*L) dz -= L;

                                    double dr2 = dx*dx + dy*dy + dz*dz;
                                    if(dr2<neighborListRadiusSquared) {
                                        atomi->neighbors.push_back(atomj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}

double Potential::pressure()
{
    return m_pressure;
}
