#include"observables.h"

//Calculate the energy per site
double energy(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}

double var_enrg(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


double energy_triangular(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_triangular2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}


double var_enrg_tri(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_tri2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


double energy_hexagonal(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_hexagonal2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}

double var_enrg_hex(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_hex2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


//Calculate magnetisation per site
double magnetisation(int g[m][n], int s, int e){
        int i,j;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        M+=g[i][j];
                }
        }
return fabs(M/((e-s+1)*n));
}

//Calculate magnetisation per site
double magnetisation2d(int g[m][n], int s, int e, int s2, int e2){
        int i,j;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        M+=g[i][j];
                }
        }
return fabs(M/((e-s+1)*(e2-s2+1)));
}

double var_mag(int g[m][n], int s, int e, double avg){
        int i,j;
        //double dM;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //dM = g[i][j]-avg;
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_mag2d(int g[m][n], int s, int e, int s2, int e2, double avg){
        int i,j;
        //double dM;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //dM = g[i][j]-avg;
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}
