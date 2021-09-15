#include"observables.h"

//Calculate the energy per site of the system
double energy(int g[m][n], const double J){
        int i,j;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E += -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
                }
        }
return E/(2.0*m*n);
}

double energy2(int g[m][n], const double J){
        int i,j;
        double dE;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
                        E+=dE*dE;
                }
        }
return E/(2.0*m*n);
}


//Calculate the energy per site of the system
double energy_triangular(int g[m][n], const double J){
        int i,j;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E+=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                }
        }
return E/(2.0*m*n);
}

double energy_triangular2(int g[m][n], const double J){
        int i,j;
        double dE;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        dE=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                E+=dE*dE;
                }
        }
return E/(2.0*m*n*m*n);
}


//Calculate the energy per site of the system
double energy_hexagonal(int g[m][n], const double J){
        int i,j;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E+=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                }
        }
return E/(2.0*m*n);
}

double energy_hexagonal2(int g[m][n], const double J){
        int i,j;
        double dE;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        dE=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        E+=dE*dE;
                }
        }
return E/(2.0*m*n*m*n);
}


//Calculate magnetisation persite
double magnetisation(int g[m][n]){
        int i,j;
        double M=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        M+=g[i][j];
                }
        }
return fabs(M/(m*n));
}

double magnetisation2(int g[m][n]){
        int i,j;
        double M=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        M+=g[i][j]*g[i][j];
                }
        }
return fabs(M/(m*n));
}


double var_mag(int g[m][n], double avg){
        int i,j;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/(m*n);
}

double var_enrg(int g[m][n], double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/(m*n);
}

double var_enrg_tri(int g[m][n], double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/(m*n);
}

double var_enrg_hex(int g[m][n], double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/(m*n);
}

