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

//average energy variance
double var_enrg(int g[x][y][z], const double J, double avg){
	int i,j,k;
	double var=0;
	double dE;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				dE = -J*(H(g[i][j][k],g[modulo(i+1,x)][j][k]) + H(g[i][j][k],g[modulo(i-1,x)][j][k]) + H(g[i][j][k],g[i][modulo(j+1,y)][k]) + H(g[i][j][k],g[i][modulo(j-1,y)][k]) + H(g[i][j][k],g[i][j][modulo(k+1,z)]) + H(g[i][j][k],g[i][j][modulo(k-1,z)]));
                                dE/=2.0;
				var += (dE-avg)*(dE-avg);
			}
		}
	}
return var/(x*y*z);
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

//average magnetisation variance
double var_mag(int g[x][y][z], double avg){
	int i,j,k;
	double var=0;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				var += (g[i][j][k] - avg)*(g[i][j][k] - avg);
			}
		}
	}
return var/(x*y*z);
}
