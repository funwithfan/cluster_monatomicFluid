#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "header.h"
#include "basic.h"
#include "readMDinformation.h"

void build_neighbor_list(struct MOL *mol, int numMol, double boxLength, double rm, struct NeighborList *neighborlist);
void updateNeighborDistance(struct MOL *mol, int numMol, double boxLength, struct NeighborList *neighborlist);