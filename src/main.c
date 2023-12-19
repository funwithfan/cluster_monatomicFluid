#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "readMDinformation.h"
#include "header.h"
#include "basic.h"
#include "graph.h"
#include "neighbor.h"
#include "time.h"


int main(int argc, char ** argv) {
    printf("Hello\n");
    int p ,t, numFrame, save_every;
    char forceField[100];
    char resultDir[FILENAME_MAX];
    char simulationDir[FILENAME_MAX];
    if(argc == 8) {
        p = atoi(argv[1]);
        t = atoi(argv[2]);
        numFrame = atoi(argv[3]);
        save_every = atoi(argv[4]);
        strcpy(resultDir, argv[5]);
        strcpy(simulationDir, argv[6]);
        strcpy(forceField, argv[7]);
    } else {
        // comment out these two lines if debugging locally
        //printf("Not enough input variables. Quiting...\n");
        //return 0; 

        printf("Local debugging:\n");
        p = 65;
        t = 155;
        numFrame = 1;
        save_every = 100;
        sprintf(resultDir, "/Users/jingcunfan/Documents/data_cp_network/Argon_localResult");
        sprintf(simulationDir, "/Users/jingcunfan/Documents/MD_simulations/Argon");
        sprintf(forceField, "ReaxFF"); // LJ or ReaxFF
    }
    printf("P = %d atm\n", p);
    printf("T = %d K\n", t);
    printf("Number of frames = %d\n", numFrame);
    printf("Saving frequency: every %d frames\n", save_every);
    printf("Result file directory: %s\n", resultDir);
    printf("Simulation file directory: %s\n", simulationDir);

    char resultDir_sub[FILENAME_MAX];
    char filename_dump[FILENAME_MAX];
    char filename_result_hill[FILENAME_MAX];
    //char filename_result_cutoff[FILENAME_MAX];

    // Create result directory
    char cmd[100];
    sprintf(resultDir_sub, "%s/%datm_fullEnergy", resultDir, p);
    sprintf(cmd, "mkdir -p %s", resultDir_sub);
    //printf("%s\n", cmd);
    system(cmd);
    sprintf(filename_result_hill, "%s/%s_%datm_%dK_hill.txt", resultDir_sub, forceField, p, t);
    //sprintf(filename_result_cutoff, "%s/%s_%datm_%dK_cutoff.txt", resultDir_sub, forceField, p, t);

    // Get number of molecules
    int numMol;
    int t0 = 250000;
    sprintf(filename_dump, "%s/%s_%datm/%dK/dump/%d.dump", simulationDir, forceField, p, t, t0);
    numMol = getNumMols(filename_dump);
    printf("Number of molecues = %d\n", numMol);

    // Create instance of structures
    struct MOL *mol = (struct MOL *)malloc(numMol * sizeof(struct MOL));
    struct Cluster *cluster = (struct Cluster *)malloc(numMol * sizeof(struct Cluster));
    struct NeighborList *neighborList = (struct NeighborList *)malloc(sizeof(struct NeighborList));

    // Define one-dimensional array to store diagonal and upper-right elements of a symmetric matrix. This is to save memory.
    int numPairs = numMol * (numMol + 1) / 2;
    double *peMatrix = (double *)malloc(numPairs * sizeof(double));
    int *adjacencyMatrix = (int *)malloc(numPairs * sizeof(int));

    int *size_count_accumulated_hill = (int *)malloc(numMol * sizeof(int));
    int *degree_count_accumulated_hill = (int *)malloc(numMol * sizeof(int));
    //int *size_count_accumulated_cutoff = (int *)malloc(numMol * sizeof(int));
    //int *degree_count_accumulated_cutoff = (int *)malloc(numMol * sizeof(int));
    initializeIntArray(numMol, size_count_accumulated_hill, 0);
    initializeIntArray(numMol, degree_count_accumulated_hill, 0);
    //initializeIntArray(numMol, size_count_accumulated_cutoff, 0);
    //initializeIntArray(numMol, degree_count_accumulated_cutoff, 0);
    
    double *PE_s_accumulated_hill = (double *)malloc(numMol * sizeof(double));
    double *KE_s_accumulated_hill = (double *)malloc(numMol * sizeof(double));
    //double *PE_s_accumulated_cutoff = (double *)malloc(numMol * sizeof(double));
    //double *KE_s_accumulated_cutoff = (double *)malloc(numMol * sizeof(double));
    initializeDoubleArray(numMol, PE_s_accumulated_hill, 0);
    initializeDoubleArray(numMol, KE_s_accumulated_hill, 0);
    //initializeDoubleArray(numMol, PE_s_accumulated_cutoff, 0);
    //initializeDoubleArray(numMol, KE_s_accumulated_cutoff, 0);


    int i, j, frame;
    double total_pe, total_ke;
    int numClusters_hill;
    //int numClusters_cutoff

    // Traverse through each snapshot
    printf("---------------------------------------------------------------------\n");
    printf("Traverse through each frame:\n");
    printf("frame\ttimestep\tke\tpe\tNc_hill\ttime\n");
    for (frame = 0; frame < numFrame; frame++) {
        tic();
        int timestep = t0 + frame * 100;

        // Load atom information  
        
        double boxLength;
        sprintf(filename_dump, "%s/%s_%datm/%dK/dump/%d.dump", simulationDir, forceField, p, t, timestep);
        readAtomInformation(filename_dump, &numMol, &boxLength, mol);
        //printf("Dump file read. Number of atoms = %d\n", numMol);

        // build or update neighbor 
        if (frame % 2 == 0) {
            build_neighbor_list(mol, numMol, boxLength, 20+3, neighborList);
            //printf("Neighbor list built, number of neighbors = %d\n", neighborList->numPair);
        } else {
            updateNeighborDistance(mol, numMol, boxLength, neighborList);
            //printf("Neighbor distance updated\n");
        }

        // Calculate energy of each molecule
        total_pe = calculatePE_perMol(forceField, mol, numMol, peMatrix, neighborList);
        //printf("pe matrix get\n");
        total_ke = calculateKE_perMol(mol, numMol);
        //printf("ke per mol get\n");


        /******************************************************************************
        * Find clusters based on Hill' energy criterion
        ******************************************************************************/
        getAdjacencyMatrix_Hill(mol, adjacencyMatrix, peMatrix, numMol, neighborList);
        //printf("Hill adjacency matrix get\n");
        numClusters_hill = findClusters(adjacencyMatrix, numMol, mol, boxLength, cluster);
        //printf("%d clusters found\n", numClusters_hill);
        calculateMolDegree(adjacencyMatrix, numMol, mol);
        //printf("nodes' degree get\n");
        

        histogramClusterSize(cluster, numClusters_hill, size_count_accumulated_hill);
        //printf("cluster size distribution get\n");
        histogramMolDegree(mol, numMol, degree_count_accumulated_hill);
        //printf("node degree distribution get\n");

        // Energy per cluster vs. cluster size
        energy_perCluster_size(mol, numClusters_hill, cluster, PE_s_accumulated_hill, KE_s_accumulated_hill);
        //inner_energy_perCluster_size(mol, numMol, numClusters_hill, cluster, PE_s_accumulated_hill, KE_s_accumulated_hill, forceField, boxLength);
        //printf("cluster energy get\n");

        
        printf("%d\t%d\t%e\t%e\t%d\t%.3f\n", frame, timestep, total_ke, total_pe, numClusters_hill, toc());

        // Save every certain frames
        if((frame + 1) % save_every == 0) {
            save_result(filename_result_hill, numMol, frame + 1, size_count_accumulated_hill, degree_count_accumulated_hill, PE_s_accumulated_hill, KE_s_accumulated_hill);
            //save_result(filename_result_cutoff, numMol, frame + 1, size_count_accumulated_cutoff, degree_count_accumulated_cutoff, PE_s_accumulated_cutoff, KE_s_accumulated_cutoff);
        }

    }

    // Save result
    save_result(filename_result_hill, numMol, numFrame, size_count_accumulated_hill, degree_count_accumulated_hill, PE_s_accumulated_hill, KE_s_accumulated_hill);
    //save_result(filename_result_cutoff, numMol, numFrame, size_count_accumulated_cutoff, degree_count_accumulated_cutoff, PE_s_accumulated_cutoff, KE_s_accumulated_cutoff);

    printf("Result averaged over all %d frames and saved to %s\n", numFrame, filename_result_hill);
    //printf("Result averaged over all %d frames and saved to %s\n", numFrame, filename_result_cutoff);

    // free memory
    free(mol);
    free(cluster);
    free(neighborList);
    free(peMatrix);
    free(adjacencyMatrix);
    free(size_count_accumulated_hill);
    free(degree_count_accumulated_hill);
    free(PE_s_accumulated_hill);
    free(KE_s_accumulated_hill);


    printf("All finished.\n");

    return 0;
}

    