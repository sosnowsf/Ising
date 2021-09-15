#include"metropolis.h"

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                        dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[i][modulo(j+1,n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

void metropolis_sweep_triangular(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                        dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r, g[i][modulo(j+1,n)]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(r,g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

void metropolis_sweep_hexagonal(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,s)+1;
                        dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep_2d(int g[m][n], int spin, int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                        dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[i][modulo(j+1,n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

void metropolis_sweep_triangular_2d(int g[m][n], int spin, int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                        dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r, g[i][modulo(j+1,n)]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(r,g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

void metropolis_sweep_hexagonal_2d(int g[m][n], int spin, int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,s)+1;
                        dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

