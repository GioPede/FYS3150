#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
using std::vector;
class Atom; class System;

class CellList
{
private:
    // This array will then have indices [cx][cy][cz][atomIndex], where the three first
    // are the 3-dimensional coordinates to the cell.
    vector<vector<vector<vector<Atom *> > > > m_cells;
public:
    double cutoffDistance = 2.5*3.405;
    int numberOfCellsPerDimension = 0;
    CellList();
    void build(System *system);
    vector<Atom*> & cell(int i, int j, int k) { return m_cells[i][j][k]; }
};
#endif // CELLLIST_H
