#include"observables.h"
#include"modulo.h"

//avg magnetisation
double magnetisation(int g[x][y][z], int s, int e){
        int i,j,k;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k];
                        }
                }
        }
return fabs(M/((e-s+1)*y*z));
}

double magnetisation2d(int g[x][y][z], int s, int e, int s2, int e2){
        int i,j,k;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k];
                        }
                }
        }
return fabs(M/((e-s+1)*(e2-s2+1)*z));
}

double magnetisation3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3){
        int i,j,k;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                M+=g[i][j][k];
                        }
                }
        }
return fabs(M/((e-s+1)*(e2-s2+1)*(e3-s3+1)));
}

//avg magnetisation variance
double var_mag(int g[x][y][z], int s, int e, double avg){
        int i,j,k;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
                        }
                }
        }
return var/((e-s+1)*y*z);
}

double var_mag2d(int g[x][y][z], int s, int e, int s2, int e2, double avg){
        int i,j,k;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*z);
}

double var_mag3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, double avg){
        int i,j,k;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*(e3-s3+1));
}

//avg energy
double energy(int g[x][y][z], const double J, int s, int e){
        int i,j,k;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                        }
                }
        }
return E/(2.0*(e-s+1)*y*z);
}

double energy2d(int g[x][y][z], const double J, int s, int e, int s2, int e2){
        int i,j,k;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                        }
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1)*z);
}

double energy3d(int g[x][y][z], const double J, int s, int e, int s2, int e2, int s3, int e3){
        int i,j,k;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                        }
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1)*(e3-s3+1));
}

//avg energy variance
double var_enrg(int g[x][y][z], int s, int e, double J, double avg){
        int i,j,k;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
                                var += (dE-avg)*(dE-avg);
                        }
                }
        }
return var/((e-s+1)*y*z);
}

double var_enrg2d(int g[x][y][z], int s, int e, int s2, int e2, double J, double avg){
        int i,j,k;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
                                var += (dE-avg)*(dE-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*z);
}

double var_enrg3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, double J, double avg){
        int i,j,k;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
                                var += (dE-avg)*(dE-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*(e3-s3+1));
}

