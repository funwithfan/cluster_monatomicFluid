#include "graph.h"

// Depth-First Search
void dfs(int node, bool *visited, int *adjMatrix, int numMol, struct MOL *mol, double boxLength, struct Cluster *currentCluster, int *clusterSize) {
    int i, j;
    int relativeImageFlag[3];
    visited[node] = true;

    (*clusterSize)++;
    currentCluster->molIDList[*clusterSize - 1] = node;

    for (i = 0; i < numMol; i++) {
        if (getIntSymmetricMatrixValue(adjMatrix, numMol, node, i) == 1 && !visited[i]) {
            findRelativeImageFlag(mol, node, i, boxLength, relativeImageFlag);
            for (j = 0; j < 3; j++) {
                mol[i].imageFlag_hill[j] = mol[node].imageFlag_hill[j] + relativeImageFlag[j];
            }
            dfs(i, visited, adjMatrix, numMol, mol, boxLength, currentCluster, clusterSize);
        }
    }
}


// find clusters and return the number of clusters
int findClusters(int *adjMatrix, int numMol, struct MOL *mol, double boxLength, struct Cluster *cluster) {
    int i, j;
    bool visited[numMol];
    int numCluster = 0;

    // Initialize visited array
    initializeBoolArray(numMol, visited, false);
    
    // Traverse through each molecule
    for (i = 0; i < numMol; i++) {
        if (!visited[i]) {
            numCluster++;
            int clusterSize = 0;
            dfs(i, visited, adjMatrix, numMol, mol, boxLength, &cluster[numCluster - 1], &clusterSize);
            cluster[numCluster - 1].size = clusterSize;
            //printf("there are %d molecules in cluster #%d\n", clusterSize, numCluster - 1);
        }
    }
    return numCluster;
}

// Calculate the degree of each mol
void calculateMolDegree(int *adjMatrix, int numMol, struct MOL *mol) {
    int i, j, degree;
    for(i = 0; i < numMol; i++) {
        degree = 0;
        for(j = 0; j < numMol; j++) {
            degree += getIntSymmetricMatrixValue(adjMatrix, numMol, i, j);
        }
        mol[i].degree = degree;
    }
}

// Cluster size histogram
void histogramClusterSize(struct Cluster *cluster, int numCluster, int *size_count_accumulated) {
    int cluster_idx;
    for(cluster_idx = 0; cluster_idx < numCluster; cluster_idx++) {
        size_count_accumulated[cluster[cluster_idx].size - 1]++;
    }
}

// Cluster size histogram
void histogramMolDegree(struct MOL *mol, int numMol, int *degree_count_accumulated) {
    int mol_idx;
    for(mol_idx = 0; mol_idx < numMol; mol_idx++) {
        degree_count_accumulated[mol[mol_idx].degree]++;
    }
}