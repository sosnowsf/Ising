#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<mpi.h>

#define m 50
#define n 50

void init_grid(int g[m][n]);
void metropolis_sweep(int g[m][n], const double , double, const double);
void metropolis_sweep_triangular(int g[m][n], const double , double, const double);
double energy(int g[m][n], const double);
double energy_triangular(int g[m][n], const double);
double magnetisation(int g[m][n]);
double specific_heat(int g[m][n], const double, double);
double specific_heat_triangular(int g[m][n], const double, double);
double susceptibility(int g[m][n], double);
void print_matrix(int g[m][n]);
void write_stats(double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
int modulo(int, int);

int main(){
	double T = 0.0; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	const double k = pow(1.38064852, -23); //Boltzmann constant
	
	double r1 = 100;
	int r2 = 100;
	int r3 = 100;

	int seed = 1997;
	srand48(seed);
	int g[m][n]; //define grid
	
	init_grid(g);

	//print_matrix(g);
	printf("%f\n", magnetisation(g));

	double *mag, *enrgy, *spec, *sus, *temp;
	mag = malloc(100*sizeof(double));
	enrgy = malloc(100*sizeof(double));
	temp = malloc(100*sizeof(double));
	spec = malloc(100*sizeof(double));
	sus = malloc(100*sizeof(double));	

	for(int k=0; k<r1; k++){
		T=0.0;
		double M=0.0;
		double E=0.0;
		for(int j=0; j<r2; j++){
			for(int i=0; i<r3; i++){
				metropolis_sweep_triangular(g,J,T,k);
			}
			double m1 = magnetisation(g)/r1;
			double e1 = energy_triangular(g,J)/r1;
			mag[j] += m1;//magnetisation(g)/r1;
			enrgy[j] += e1;//energy(g,J)/r1;
			M += m1*m1;//(magnetisation(g)/r1)*(magnetisation(g)/r1);
			E += e1*e1;//(energy(g,J)/r1)*(energy(g,J)/r1);
			spec[j] += (E-enrgy[j]*enrgy[j])*(T*T)/r1;//specific_heat(g,J,T)/r1;
			sus[j] += (M-mag[j]*mag[j])*T/r1;//susceptibility(g,T)/r1;
			if(k==0) temp[j] = T;
			T+=0.05;
			//print_matrix(g);
			init_grid(g);
		}
	}
	write_stats(mag,enrgy,spec,sus,temp,r1);
	//printf("%f\n", magnetisation(g));
	//g[3][5]=9;
	//print_matrix(g);
	free(mag);
	free(enrgy);
	free(temp);
return 0;
}

//Initialise grid for random start
void init_grid(int g[m][n]){
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			double U = drand48();
			if(U<0.5) g[i][j]=-1;
			else{ 
				g[i][j]=1;
			}
		}
	}
	//print_matrix(g);
}

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], const double J, double T, const double k){
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			//Calculate the energy difference and check if spin has to be flipped
			double dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
			//double dE = J*g[i][j]*(g[(i+1)%m][j]+g[(i-1)&m][j]+g[i][(j+1)%n]+g[i][(j-1)%n]);
			if(dE<0) g[i][j] *= -1;
			else{
				if(drand48()<(exp(-dE/T))){
					//printf("%f %f %f\n", dE, T, exp(-dE/T));
					g[i][j] *= -1;
				}
			}
		}
	}
}

//Sweep for triangular lattice
void metropolis_sweep_triangular(int g[m][n], const double J, double T, const double k){
        int i,j;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        double dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        //double dE = J*g[i][j]*(g[(i+1)%m][j]+g[(i-1)&m][j]+g[i][(j+1)%n]+g[i][(j-1)%n]);
                        if(dE<0) g[i][j] *= -1;
                        else{
                                if(drand48()<(exp(-dE/T))){
                                        //printf("%f %f %f\n", dE, T, exp(-dE/T));
                                        g[i][j] *= -1;
                                }
                        }
                }
        }
}

//Calculate the energy per site of the system
double energy(int g[m][n], const double J){
	int i,j;
	double E=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
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
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*m*n);
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
return fabs(M)/(m*n);
}

double specific_heat(int g[m][n], const double J, double T){
	int i,j;
	double E;
	double E1=0;
	double E2=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			E=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
			E1+=E;
			E2+=E*E;
		}
	}
	E1/=(2.0*m*n);
	E2/=(2.0*m*n);
	return (E2-(E1*E1));//*T*T;
}

double specific_heat_triangular(int g[m][n], const double J, double T){
        int i,j;
        double E;
        double E1=0;
        double E2=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[i+1][j+((int)pow(-1,i))] + g[i-1][j+((int)pow(-1,i))]);
                        E1+=E;
                        E2+=E*E;
                }
        }
        E1/=(2.0*m*n);
        E2/=(2.0*m*n);
        return (E2-(E1*E1));//*T*T;
}


double susceptibility(int g[m][n], double T){
	int i,j;
	double M1=0;
	double M2=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			M1+=g[i][j];
			M2+=g[i][j]*g[i][j];
		}
	}
	M1/=(2.0*m*n);
	M2/=(2.0*m*n);
	return (M2-(M1*M1))/T;

}

void print_matrix(int g[m][n]){
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%d ", g[i][j]);
		}
	printf("\n");
	}
}

//Reset grid after every temp increase
void reset_grid(int G[m][n], int g[m][n]){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			g[i][j]=G[i][j];
		}
	}
}

void write_stats(double *mag, double *energy, double *sus, double *spec, double *T, int r1){
	FILE *fp = fopen("stats_serial.txt", "w");
	for(int i=0; i<r1; i++){
		fprintf(fp,"%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}

int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}
