#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <stdbool.h>
#include "header.h"
#include "basic.h"
#include "readMDinformation.h"

void dfs(int node, bool *visited, int *adjMatrix, int numMol, struct MOL *mol, double boxLength, struct Cluster *currentCluster, int *clusterSize);
int findClusters(int *adjMatrix, int numMol, struct MOL *mol, double boxLength, struct Cluster *cluster);
void calculateMolDegree(int *adjMatrix, int numMol, struct MOL *mol);
void histogramClusterSize(struct Cluster *cluster, int numCluster, int *size_count_accumulated);
void histogramMolDegree(struct MOL *mol, int numMol, int *degree_count_accumulated);