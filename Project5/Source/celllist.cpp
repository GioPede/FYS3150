#include "celllist.h"
#include "system.h"
#include "atom.h"
#include <iostream>
using namespace std;

CellList::CellList()
{

}

void CellList::build(System *system)
{
    numberOfCellsPerDimension = system->systemSize().x() / cutoffDistance;
    cutoffDistance = system->systemSize().x() / numberOfCellsPerDimension;

    // Resize so our cell vectors have the correct size
    m_cells.resize(numberOfCellsPerDimension);
    for(int i=0; i<numberOfCellsPerDimension; i++) {
        m_cells[i].resize(numberOfCellsPerDimension);
        for(int j=0; j<numberOfCellsPerDimension; j++) {
            m_cells[i][j].resize(numberOfCellsPerDimension);
            for(int k=0; k<numberOfCellsPerDimension; k++) {
                m_cells[i][j][k].clear();
            }
        }
    }

    // Compute the cell dimensions
    double cellSize = system->systemSize().x()/numberOfCellsPerDimension;
    for(Atom *atom : system->atoms()) {
        int i = atom->position.x()/cellSize;
        int j = atom->position.y()/cellSize;
        int k = atom->position.z()/cellSize;
        m_cells.at(i).at(j).at(k).push_back(atom);
    }
}
