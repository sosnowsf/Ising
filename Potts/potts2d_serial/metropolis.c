#include"metropolis.h"

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int s, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(k=0; k<m*n; k++){
                i=gsl_rng_uniform_int(gsl_mt, m);
                j=gsl_rng_uniform_int(gsl_mt, n);
                r=(gsl_rng_uniform_int(gsl_mt, s))+1;
                dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[i][modulo(j+1,n)]));
                if(dE<=0) g[i][j] = r;
                else{
                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] = r;
                }
        }
}

//Sweep for triangular lattice
void metropolis_sweep_triangular(int g[m][n], int s, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(k=0; k<m*n; k++){
                        i=gsl_rng_uniform_int(gsl_mt, m);
                        j=gsl_rng_uniform_int(gsl_mt, n);
                        r=gsl_rng_uniform_int(gsl_mt, s)+1;
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r, g[i][modulo(j+1,n)]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(r,g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)])) ;
                        if(dE<=0) g[i][j] = r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] = r;
                        }
        }
}

//Sweep for triangular lattice
void metropolis_sweep_hexagonal(int g[m][n], int s, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(k=0; k<m*n; k++){
                        i=gsl_rng_uniform_int(gsl_mt, m);
                        j=gsl_rng_uniform_int(gsl_mt, n);
                        r=gsl_rng_uniform_int(gsl_mt, s)+1;
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        if(dE<=0) g[i][j] = r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] = r;
                        }
        }
}

