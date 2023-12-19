#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "header.h"
#include "basic.h"

int getNumMols(char *filename);
void readAtomInformation(char *filename, int *numMol, double *boxLength, struct MOL *mol);
void findRelativeImageFlag(struct MOL *mol, int mol1_id, int mol2_id, double boxLength, int *relativeImageFlag);
double trueDistance(double distance, double box_length);
double interAtomDistance(struct MOL *mol, int mol1_id, int mol2_id, double boxLength);
//void getDistanceMatrix(struct MOL *mol, int numMol, double *distanceMatrix, double boxLength);
double intermolecularPE(char *forceField, double r_ij);
double calculatePE_perMol(char *forceField, struct MOL *mol, int numMol, double *peMatrix, struct NeighborList *neighborList);
double calculateKE_perMol(struct MOL *mol, int numMol);
double relativeKE(struct MOL *mol, int mol1_id, int mol2_id);
void getAdjacencyMatrix_Hill(struct MOL *mol, int *adjacencyMatrix, double *peMatrix, int numMol, struct NeighborList *neighborList);
void getAdjacencyMatrix_distance(struct MOL *mol, int *adjacencyMatrix, int numMol, double cutoff, struct NeighborList *neighborList);
void energy_perCluster_size(struct MOL *mol, int numClusters, struct Cluster *cluster, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated);
void inner_energy_perCluster_size(struct MOL *mol, int numMol, int numClusters, struct Cluster *cluster, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated, char *forceField, double boxLength);
void save_result(char *filename, int max_size, int numFrame, int *size_count_accumulated, int *degree_count_accumulated, double *PE_s_accumulated, double *KE_s_accumulated);
