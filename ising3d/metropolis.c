#include"metropolis.h"
#include"modulo.h"

//Metropolis sweep 1D decomposition
void metropolis_sweep(int g[x][y][z], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k;
        double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                if(dE<0) g[i][j][k]*=-1;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
                                                g[i][j][k]*=-1;
                                        }
                                }
                        }
                }
        }
}

//2D decomposition version
void metropolis_sweep2d(int g[x][y][z], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                if(dE<0) g[i][j][k]*=-1;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
                                                g[i][j][k]*=-1;
                                        }
                                }
                        }
                }
        }
}

//3D decomposition version
void metropolis_sweep3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                if(dE<0) g[i][j][k]*=-1;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
                                                g[i][j][k]*=-1;
                                        }
                                }
                        }
                }
        }
}


