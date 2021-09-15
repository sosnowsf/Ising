#include"observables.h"

double energy(int g[x][y][z], const double J){
        int i,j,k;
        double E=0;
        for(i=0; i<x; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                E += -J*(H(g[i][j][k],g[modulo(i+1,x)][j][k]) + H(g[i][j][k],g[modulo(i-1,x)][j][k]) + H(g[i][j][k],g[i][modulo(j+1,y)][k]) + H(g[i][j][k],g[i][modulo(j-1,y)][k]) + H(g[i][j][k],g[i][j][modulo(k+1,z)]) + H(g[i][j][k],g[i][j][modulo(k-1,z)]));
                        }
                }
        }
return E/(2.0*x*y*z);
}

double energy2(int g[x][y][z], const double J){
        int i,j,k;
        double dE;
        double E=0;
        for(i=0; i<x; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                dE = -J*(H(g[i][j][k],g[modulo(i+1,x)][j][k]) + H(g[i][j][k],g[modulo(i-1,x)][j][k]) + H(g[i][j][k],g[i][modulo(j+1,y)][k]) + H(g[i][j][k],g[i][modulo(j-1,y)][k]) + H(g[i][j][k],g[i][j][modulo(k+1,z)]) + H(g[i][j][k],g[i][j][modulo(k-1,z)]));
                                E+=dE*dE;
                        }
                }
        }
return E/(2.0*x*y*z);
}


double magnetisation(int g[x][y][z]){
        int i,j,k;
        double M=0;
        for(i=0; i<x; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k];
                        }
                }
        }
return fabs(M)/(x*y*z);
}

double magnetisation2(int g[x][y][z]){
        int i,j,k;
        double M=0;
        for(i=0; i<x; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k]*g[i][j][k];
                        }
                }
        }
return fabs(M)/(x*y*z*x*y*z);
}

