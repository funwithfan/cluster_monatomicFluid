#include "readMDinformation.h"

// force field parameters
double mass = 39.948;
double charge = 0;

double lj_rc = 20;
double lj_eps = 0.224835086042065;
double lj_sig = 3.3677;

double reaxff_rc = 10;
double reaxff_p_vdw = 1.5591;
double reaxff_eps = 0.4060 ;
double reaxff_r_vdw = 1.9050;
double reaxff_alpha = 10.3000;
double reaxff_gamma_vdw = 40.0000;
double reaxff_gamma_coul = 0.5000;


// unit conversion
double NA = 6.02e23;
double kcal_mol2j = 6.9477e-21; 
double A2m = 1e-10;
double e2c = 1.60217663e-19;
double a_fs2m_s = 1e5;
double g_mol2kg = 1.6611e-27;


int getNumMols(char *filename) {
    FILE *fpr;
    int numMol;
    char str[100];
    
    fpr = fopen(filename, "r");
    
    // read lines of dump file
    fgets(str,100,fpr); // TIMESTEP
    fgets(str,100,fpr); // timestep
    fgets(str,100,fpr); // ITEM: NUMBER OF ATOMS
    fgets(str,100,fpr); // number of atoms
    numMol = atoi(str);

    fclose(fpr);

    return numMol;
}


void readAtomInformation(char *filename, int *numMol, double *boxLength, struct MOL *mol) {
    FILE *fpr;
    int i, j, index, mol_id;
    double lo, hi;
    char str[100];
    
    fpr = fopen(filename, "r");
    

    // read lines of dump file
    fgets(str,100,fpr); // TIMESTEP
    fgets(str,100,fpr); // timestep
    fgets(str,100,fpr); // ITEM: NUMBER OF ATOMS
    fgets(str,100,fpr); // number of atoms
    *numMol = atoi(str);
    fgets(str,100,fpr); // ITEM: BOX BOUNDS pp pp pp
    fscanf(fpr,"%s",str); // xlo
    lo = atof(str);
    fscanf(fpr,"%s",str); // xhi
    hi = atof(str);
    *boxLength = hi - lo;
    fgets(str,100,fpr); // skip line?
    fgets(str,100,fpr); // ylo yhi
    fgets(str,100,fpr); // zlo zhi
    fgets(str,100,fpr); // ITEM: ATOMS id mol type x y z ix iy iz vx vy vz element

    // read atom informations
    for(i = 0; i < *numMol; i++) {
        fscanf(fpr,"%s",str); // id
        mol_id = atof(str) - 1;
        fscanf(fpr, "%s", str); // type
        fscanf(fpr, "%s", str); // x
        mol[mol_id].position[0] = atof(str);
        fscanf(fpr, "%s", str); // y
        mol[mol_id].position[1] = atof(str);
        fscanf(fpr, "%s", str); // z
        mol[mol_id].position[2] = atof(str);
        fscanf(fpr, "%s", str); // vx
        mol[mol_id].velocity[0] = atof(str);
        fscanf(fpr, "%s", str); // vy
        mol[mol_id].velocity[1] = atof(str);
        fscanf(fpr, "%s", str); // vz
        mol[mol_id].velocity[2] = atof(str);
        fscanf(fpr, "%s", str); // element

        // initialze image flags
        mol[mol_id].imageFlag_hill[0] = 0;
        mol[mol_id].imageFlag_hill[1] = 0;
        mol[mol_id].imageFlag_hill[2] = 0;
        //mol[mol_id].imageFlag_cutoff[0] = 0;
        //mol[mol_id].imageFlag_cutoff[1] = 0;
        //mol[mol_id].imageFlag_cutoff[2] = 0;
    
    }
    fclose(fpr);
}


// Find the true distance considering PBC
double trueDistance(double distance, double box_length) {
    if (distance > box_length/2) {
        distance -= box_length;
    } else if (distance < -box_length/2) {
        distance += box_length;
    }
    return distance;
}

// Find relative image flag of two molecules
void findRelativeImageFlag(struct MOL *mol, int mol1_id, int mol2_id, double boxLength, int *relativeImageFlag) {
    int i;
    double d;
    for (i = 0; i < 3; i++) {
        relativeImageFlag[i] = 0;
        d = mol[mol1_id].position[i] - mol[mol2_id].position[i];
        if (d > boxLength/2) {
            relativeImageFlag[i] = 1;
        } else if (d < -boxLength/2) {
            relativeImageFlag[i] = -1;
        }
    }
}

// Return the distance between two atoms considering PBC
double interAtomDistance(struct MOL *mol, int mol1_id, int mol2_id, double boxLength) {
    int i;
    double d, r_sq;
    r_sq = 0;
    for (i = 0; i < 3; i++) {
        d = mol[mol1_id].position[i] - mol[mol2_id].position[i];
        if (d > boxLength/2) {
            d -= boxLength;
        } else if (d < -boxLength/2) {
            d += boxLength;
        }
        r_sq += d * d;
    }
    return sqrt(r_sq);
}





/*
void getDistanceMatrix(struct MOL *mol, int numMol, double *distanceMatrix, double boxLength) {
    int mol1_id, mol2_id;
    double dx, dy, dz, r;

    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id + 1; mol2_id < numMol; mol2_id++) {
            dx = mol[mol1_id].position[0] - mol[mol2_id].position[0];
            dy = mol[mol1_id].position[1] - mol[mol2_id].position[1];
            dz = mol[mol1_id].position[2] - mol[mol2_id].position[2];
            dx = trueDistance(dx, boxLength);
            dy = trueDistance(dy, boxLength);
            dz = trueDistance(dz, boxLength);
            r = sqrt(dx*dx + dy*dy + dz*dz);
            setDoubleSymmetricMatrixValue(distanceMatrix, numMol, mol1_id, mol2_id, r);
            //distanceMatrix[mol1_id][mol2_id] = r;
            //distanceMatrix[mol2_id][mol1_id] = r;
        }
    }
}
*/



// Calculate inter-molecular potential energy between two molecules, unit: kcal/mol
double intermolecularPE(char *forceField, double r_ij) {
    double pe, vdw, coul;
    int i, j;

    if (!strcmp(forceField, "LJ")) {
        // van der Waals
        vdw = 4 * lj_eps * (pow(lj_sig/r_ij, 12) - pow(lj_sig/r_ij, 6));

        // Coulomb
        coul = 0;
        //double Pi = 3.1415926;
        //double eps_0 = 8.85418782e-12; // F/m, dielectric permittivity
        //double C = 1 / (4 * Pi * eps_0); // energy-conversion constant
        //int dielectric = 1 ; // dielectric constant
        //coul = C * charge * charge / (dielectric * r_ij); 
        //coul = coul * e2c * e2c / A2m / kcal_mol2j
    }
    else if (!strcmp(forceField, "ReaxFF")) {
        // Taper correction
        double Tap = 20.0 * pow(r_ij / reaxff_rc, 7) - 70.0 * pow(r_ij / reaxff_rc, 6) + 84.0 * pow(r_ij / reaxff_rc, 5) - 35.0 * pow(r_ij / reaxff_rc, 4) + 1;
        
        // van der Waals
        double D_ij = reaxff_eps;
        double alpha_ij = reaxff_alpha;
        double r_vdw_ij = 2.0 * sqrt(reaxff_r_vdw * reaxff_r_vdw);
        //if(r_vdw[i] != r_vdw[j]) {
        //    r_vdw_ij = 1.0 * sqrt(r_vdw[i] * r_vdw[j]); // don't know why, but this seems to be consistent with MD results.
        //}
        double gamma_vdw_ij = reaxff_gamma_vdw;
        double f = pow(pow(r_ij, reaxff_p_vdw) + pow(1.0 / gamma_vdw_ij, reaxff_p_vdw), 1.0 / reaxff_p_vdw);
        double exp1 = exp(alpha_ij * (1 - f / r_vdw_ij));
        double exp2 = exp(0.5 * alpha_ij * (1 - f / r_vdw_ij));
        vdw = Tap * D_ij * (exp1 - 2 * exp2);

        // Coulomb
        coul = 0;
        //double gamma_coul_ij_inverse = 1.0 / sqrt(gamma_coul[i] * gamma_coul[j]);
        //coul = coul + Tap * C * mol[mol1].charge[i] * mol[mol2].charge[j] / pow((pow(r[i][j], 3) + pow(gamma_coul_ij_inverse, 3)), 0.333333333);
        //coul = coul * e2c * e2c / A2m / kcal_mol2j;
    }


    pe = vdw + coul ;
    //printf("pe = %f\n", pe);
    return pe;
}


// Calculate inter-molecular potential energy of each molecule and return total pe, unit: kcal/mol
double calculatePE_perMol(char *forceField, struct MOL *mol, int numMol, double *peMatrix, struct NeighborList *neighborList) {
    int pair_idx, mol1_id, mol2_id;
    double pe, total_pe;
    double rc, r;

    total_pe = 0;

    if (!strcmp(forceField, "LJ")) {
        rc = lj_rc;
    }
    else if (!strcmp(forceField, "ReaxFF")) {
        rc = reaxff_rc;
    }

    // Initialize PE of each molecule and initialize peMatrix
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        mol[mol1_id].PE = 0;
        for (mol2_id = mol1_id; mol2_id < numMol; mol2_id++) {
            setDoubleSymmetricMatrixValue(peMatrix, numMol, mol1_id, mol2_id, 0);
        }
    }
    
    // Traverse through each pair
    for (pair_idx = 0; pair_idx < neighborList->numPair; pair_idx++) {
        mol1_id = neighborList->pairID[pair_idx][0];
        mol2_id = neighborList->pairID[pair_idx][1];
        r = neighborList->pairDistance[pair_idx];
        if (r <= rc) {
            pe = intermolecularPE(forceField, r);
        }
        else {
            pe = 0;
        }
        setDoubleSymmetricMatrixValue(peMatrix, numMol, mol1_id, mol2_id, pe);
        mol[mol1_id].PE += 0.5 * pe;
        mol[mol2_id].PE += 0.5 * pe;
        total_pe += pe;
    }
    return total_pe;
}

// Calculate kinetic energy of each molecule and return total ke, unit: kcal/mol
double calculateKE_perMol(struct MOL *mol, int numMol) {
    int molId, direction;
    double mv2;
    double total_ke = 0;

    for (molId = 0; molId < numMol; molId++) {
        mv2 = 0;
        for (direction = 0; direction < 3; direction++) {
            mv2 += mass * mol[molId].velocity[direction] * mol[molId].velocity[direction];
        }
        mol[molId].KE = 0.5 * mv2 * g_mol2kg * a_fs2m_s * a_fs2m_s / kcal_mol2j;
        total_ke += mol[molId].KE;
    }
    return total_ke;
}


// Calculate relative kinetic energy between two molecules, unit: kcal/mol
double relativeKE(struct MOL *mol, int mol1_id, int mol2_id) {
    double v_relative_square, relative_ke;
    double v_relative[3];
    int i;
    v_relative_square = 0;

    for (i = 0; i < 3; i++) { // each direction
        v_relative[i] = mol[mol1_id].velocity[i] - mol[mol2_id].velocity[i];
        v_relative_square += v_relative[i] * v_relative[i];
    }
    relative_ke = 0.5 * mass * v_relative_square * g_mol2kg * a_fs2m_s * a_fs2m_s / kcal_mol2j;
    return relative_ke;
}

// Calculate adjacency matrix accoding to Hill's energy criterion
void getAdjacencyMatrix_Hill(struct MOL *mol, int *adjacencyMatrix, double *peMatrix, int numMol, struct NeighborList *neighborList) {
    int pair_idx, mol1_id, mol2_id;
    double ke, pe;
    double dx, dy, dz;

    // Initialize adjacencyMatrix
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id; mol2_id < numMol; mol2_id++) {
            setIntSymmetricMatrixValue(adjacencyMatrix, numMol, mol1_id, mol2_id, 0);
        }
    }

    // Traverse through each pair
    for (pair_idx = 0; pair_idx < neighborList->numPair; pair_idx++) {
        mol1_id = neighborList->pairID[pair_idx][0];
        mol2_id = neighborList->pairID[pair_idx][1];
        pe = getDoubleSymmetricMatrixValue(peMatrix, numMol, mol1_id, mol2_id);
        if (pe < 0) {
            ke = relativeKE(mol, mol1_id, mol2_id);
            if(ke + pe < 0) {
                setIntSymmetricMatrixValue(adjacencyMatrix, numMol, mol1_id, mol2_id, 1);
            }
        }
    }
}

// Calculate adjacency matrix accoding to distance cutoff
void getAdjacencyMatrix_distance(struct MOL *mol, int *adjacencyMatrix, int numMol, double cutoff, struct NeighborList *neighborList) {
    int pair_idx, mol1_id, mol2_id;
    double dx, dy, dz, r;

    // Initialize adjacencyMatrix
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id; mol2_id < numMol; mol2_id++) {
            setIntSymmetricMatrixValue(adjacencyMatrix, numMol, mol1_id, mol2_id, 0);
        }
    }

    // Traverse through each pair
    for (pair_idx = 0; pair_idx < neighborList->numPair; pair_idx++) {
        mol1_id = neighborList->pairID[pair_idx][0];
        mol2_id = neighborList->pairID[pair_idx][1];
        r = neighborList->pairDistance[pair_idx];
        if (r < cutoff) {
            setIntSymmetricMatrixValue(adjacencyMatrix, numMol, mol1_id, mol2_id, 1);
        }
    }
}

// Energy per cluster as a function of cluster size
void energy_perCluster_size(struct MOL *mol, int numClusters, struct Cluster *cluster, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated) {
    int i, cluster_id, size, molID;

    for (cluster_id = 0; cluster_id < numClusters; cluster_id++) {
        size = cluster[cluster_id].size;
        for (i = 0; i < size; i++) {
            molID = cluster[cluster_id].molIDList[i];
            mol[molID].clusterID = cluster_id;
            mol[molID].clusterSize = size;
            PE_perCluster_size_accumulated[size - 1] += mol[molID].PE;
            KE_perCluster_size_accumulated[size - 1] += mol[molID].KE;
        }
    }
}

void inner_energy_perCluster_size(struct MOL *mol, int numMol, int numClusters, struct Cluster *cluster, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated, char *forceField, double boxLength) {
    int i, j, cluster_id, size, molID, mol2ID;
    double inner_pe;
    double dx, dy, dz, r, rc;

    if (!strcmp(forceField, "LJ")) {
        rc = lj_rc;
    }
    else if (!strcmp(forceField, "ReaxFF")) {
        rc = reaxff_rc;
    }

    for (cluster_id = 0; cluster_id < numClusters; cluster_id++) {
        size = cluster[cluster_id].size;
        inner_pe = 0;
        for (i = 0; i < size; i++) {
            molID = cluster[cluster_id].molIDList[i];
            mol[molID].clusterID = cluster_id;
            mol[molID].clusterSize = size;
            KE_perCluster_size_accumulated[size - 1] += mol[molID].KE;
            for (j = i + 1; j < size; j++) {
                mol2ID = cluster[cluster_id].molIDList[j];
                dx = mol[molID].position[0] + boxLength * mol[molID].imageFlag_hill[0] - (mol[mol2ID].position[0] + boxLength * mol[mol2ID].imageFlag_hill[0]);
                dy = mol[molID].position[1] + boxLength * mol[molID].imageFlag_hill[1] - (mol[mol2ID].position[1] + boxLength * mol[mol2ID].imageFlag_hill[1]);
                dz = mol[molID].position[2] + boxLength * mol[molID].imageFlag_hill[2] - (mol[mol2ID].position[2] + boxLength * mol[mol2ID].imageFlag_hill[2]);
                r = sqrt(dx*dx + dy*dy + dz*dz);
                if (r < rc) {
                    inner_pe += intermolecularPE(forceField, r);
                }
            }
        }
        PE_perCluster_size_accumulated[size - 1] += inner_pe;
    }
}

void save_result(char *filename, int max_size, int numFrame, int *size_count_accumulated, int *degree_count_accumulated, double *PE_s_accumulated, double *KE_s_accumulated) {
    int i, size, degree;
    FILE *fpw;
    fpw = fopen(filename, "w");
    double size_count_mean, degree_count_mean, PE_s_mean, KE_s_mean;

    // Write header
    fprintf(fpw, "number of frames = %d\n", numFrame);
    fprintf(fpw, "s\tN_c\tPE(s)\tKE(s)\tk\tN_k\n");

    // Write data
    for (i = 0; i < max_size; i++) {
        size = i + 1;
        degree = i;
        size_count_mean = (1.0 * size_count_accumulated[i]) / (1.0 * numFrame);
        degree_count_mean = (1.0 * degree_count_accumulated[i]) / (1.0 * numFrame);
        PE_s_mean = PE_s_accumulated[i] / (1.0 * size_count_accumulated[i]);
        KE_s_mean = KE_s_accumulated[i] / (1.0 * size_count_accumulated[i]);

        fprintf(fpw, "%d\t%e\t%e\t%e\t%d\t%e\n", size, size_count_mean, PE_s_mean, KE_s_mean, degree, degree_count_mean);
    }

    fclose(fpw);
}






void unwrappPosition(struct MOL *mol, int *molIDs, int numMol) {

}