#include "neighbor.h"

void build_neighbor_list(struct MOL *mol, int numMol, double boxLength, double rm, struct NeighborList *neighborlist) {
    int mol1_id, mol2_id;
    double r;
    int pair_idx = 0;
    
    // check each pair
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id + 1; mol2_id < numMol; mol2_id++) {
            r = interAtomDistance(mol, mol1_id, mol2_id, boxLength);
            if (r < rm) {
                neighborlist->pairID[pair_idx][0] = mol1_id;
                neighborlist->pairID[pair_idx][1] = mol2_id;
                neighborlist->pairDistance[pair_idx] = r;
                pair_idx++;
            }
        }
    }
    neighborlist->numPair = pair_idx;
    //neighborlist.numPair = pair_idx;
}

void updateNeighborDistance(struct MOL *mol, int numMol, double boxLength, struct NeighborList *neighborlist) {
    int pair_idx, mol1_id, mol2_id;
    
    for (pair_idx = 0; pair_idx < neighborlist->numPair; pair_idx++) {
        mol1_id = neighborlist->pairID[pair_idx][0];
        mol2_id = neighborlist->pairID[pair_idx][1];
        neighborlist->pairDistance[pair_idx] = interAtomDistance(mol, mol1_id, mol2_id, boxLength);
    }
    
}