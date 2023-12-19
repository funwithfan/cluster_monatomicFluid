#ifndef __COMMON_H__
#define __COMMON_H__
#define N 25000 // maximum number of molecues
struct MOL {
    double position[3];
    double velocity[3];
    double PE;
    double KE;
    int clusterID;
    int clusterSize;
    int imageFlag_hill[3];
    //int imageFlag_cutoff[3];
    int degree;
};

struct Cluster {
    int molIDList[N];
    int size;
};

struct NeighborList {
    int numPair;
    int pairID[N*1000][2];
    double pairDistance[N*1000];
};




/*
// unit conversion
double NA;
double kcal_mol2j; 
double A2m;
double e2c;
double a_fs2m_s;
double g_mol2kg;

// force field parameters
double mass;
double charge;

double lj_rc;
double lj_eps;
double lj_sig;

double reaxff_rc;
double reaxff_p_vdw;
double reaxff_eps;
double reaxff_r_vdw;
double reaxff_alpha;
double reaxff_gamma_vdw;
double reaxff_gamma_coul;
*/

#endif


