#include"metropolis.h"

void metropolis_sweep(int g[x][y][z], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                                dE = 2.0*J*(H(r,g[modulo(i+1,x)][j][k]) + H(r,g[modulo(i-1,x)][j][k]) + H(r,g[i][modulo(j+1,y)][k]) + H(r,g[i][modulo(j-1,y)][k]) + H(r,g[i][j][modulo(k+1,z)]) + H(r,g[i][j][modulo(k-1,z)]));
                                if(dE<=0) g[i][j][k]=r;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))) g[i][j][k]=r;
                                }
                        }
                }
        }
}

void metropolis_sweep_2d(int g[x][y][z], int spin, int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                                dE = 2.0*J*(H(r,g[modulo(i+1,x)][j][k]) + H(r,g[modulo(i-1,x)][j][k]) + H(r,g[i][modulo(j+1,y)][k]) + H(r,g[i][modulo(j-1,y)][k]) + H(r,g[i][j][modulo(k+1,z)]) + H(r,g[i][j][modulo(k-1,z)]));
                                if(dE<=0) g[i][j][k]=r;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))) g[i][j][k]=r;
                                }
                        }
                }
        }
}

void metropolis_sweep_3d(int g[x][y][z], int spin, int s, int e, int s2, int e2, int s3, int e3, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k,r;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                r=gsl_rng_uniform_int(gsl_mt,spin)+1;
                                dE = 2.0*J*(H(r,g[modulo(i+1,x)][j][k]) + H(r,g[modulo(i-1,x)][j][k]) + H(r,g[i][modulo(j+1,y)][k]) + H(r,g[i][modulo(j-1,y)][k]) + H(r,g[i][j][modulo(k+1,z)]) + H(r,g[i][j][modulo(k-1,z)]));
                                if(dE<=0) g[i][j][k]=r;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))) g[i][j][k]=r;
                                }
                        }
                }
        }
}

