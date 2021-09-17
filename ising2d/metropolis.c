#include"metropolis.h"

//Sweep through the square lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

//2D decomposition version of above function
void metropolis_sweep2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

//Triangular lattice metropolis sweep
void metropolis_sweep_triangular(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]  + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

//2D decomp version of above function
void metropolis_sweep_triangular2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

//Hex lattice metropolis sweep
void metropolis_sweep_hexagonal(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

//2D decomp version of above function
void metropolis_sweep_hexagonal2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

