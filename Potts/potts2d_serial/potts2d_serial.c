/*#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define m 50
#define n 50
*/
#include"grid_size.h"
#include"metropolis.h"
#include"observables.h"
#include"metropolis.h"

void init_grid(int g[m][n], int, int,  gsl_rng *);
/*
void metropolis_sweep(int g[m][n], int, const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], int, const double , double, gsl_rng *);
*/
/*
double energy(int g[m][n], const double);
double energy2(int g[m][n], const double);
double energy_triangular(int g[m][n], const double);
double energy_triangular2(int g[m][n], const double);
double energy_hexagonal(int g[m][n], const double);
double energy_hexagonal2(int g[m][n], const double);
double magnetisation(int g[m][n]);
double magnetisation2(int g[m][n]);
*/
void print_matrix(int g[m][n]);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
/*
int modulo(int, int);
*/
/*
int delta(int, int);
int H(int, int);
*/
/*
double var_mag(int g[m][n], double);
double var_enrg(int g[m][n], double, double);
double var_enrg_tri(int g[m][n], double, double);
double var_enrg_hex(int g[m][n], double, double);
*/

int main(int argc, char *argv[]){

	int spin =  5;
	//Parse variables
	int sq=0;
	int triangle=0;
	int hex=0;
	int c=0;
	int opt;
	while((opt = getopt(argc, argv, "123ch")) != -1){
		switch(opt){
			case '1':
				sq = 1;
				break;
			case '2':
				triangle = 1;
				break;
			case '3':
				hex = 1;
				break;
			case 'h':
				printf("Command lines:\n -1 for square grid\n -2 for triangular grid \n -3 for hexagonal grid\n");
				return 0;
			case 'c':
				c=1;
				break;
		}
	}
	//Error checking
	if((sq+triangle+hex)==0){
		perror("No parse command given, use -h for help");
		_exit(EXIT_FAILURE);
	}
	if((sq+triangle+hex)>1){
		perror("Please provide only one parse command at a time");
		_exit(EXIT_FAILURE);
	}
	double T;
	const double J = 1.0; //Coupling constant
	//const double k = pow(1.38064852, -23); //Boltzmann constant
	
	double r1 = 50.0; //Number of times simulation is ran to ge the average
	int r2 = 100; //Number of temperatures
	int r3 = 1000; //Number of sweeps

	//int seed = 1997;
	//srand48(seed);
	int g[m][n]; //define grid
	
	/*Set up rng*/
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);

	init_grid(g, c, spin, gsl_mt);

	//printf("%f %f\n", gsl_ran_gaussian(gsl_mt,1), gsl_ran_gaussian(gsl_mt,1));

	//print_matrix(g);
	//printf("%f\n", magnetisation(g));

	double *mag, *enrgy, *spec, *sus, *temp;
	mag = malloc(r2*sizeof(double));
	enrgy = malloc(r2*sizeof(double));
	temp = malloc(r2*sizeof(double));
	spec = malloc(r2*sizeof(double));
	sus = malloc(r2*sizeof(double));	

        if(triangle==1){
        for(int k=0; k<r1; k++){
                T=0.0001;
                double M=0.0;
                double E=0.0;
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_triangular(g,spin,J,t,gsl_mt);
                        }
                        double m1 = magnetisation(g);
                        double e1 = energy_triangular(g,J);
                        mag[j] += m1/r1;//magnetisation(g)/r1;
                        enrgy[j] += e1/r1;//energy(g,J)/r1;
                        M = var_mag(g,m1);//(magnetisation(g)/r1)*(magnetisation(g)/r1);
                        E = var_enrg_tri(g,J,e1);//(energy(g,J)/r1)*(energy(g,J)/r1);
                        spec[j] += E*t*t/r1;// fabs(E-(e1*e1))*(T*T)/r1;//specific_heat(g,J,T)/r1;
                        sus[j] += M*t/r1;//fabs(M-(m1*m1))*T/r1;//susceptibility(g,T)/r1;
                        if(k==0) temp[j] = T;
                        T+=0.1;
                        init_grid(g,c,spin,gsl_mt);
                }
        }
        }


        if(hex==1){
        for(int k=0; k<r1; k++){
                T=0.0001;
                double M=0.0;
                double E=0.0;
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_hexagonal(g,spin,J,t,gsl_mt);
                        }
                        double m1 = magnetisation(g);
                        double e1 = energy_hexagonal(g,J);
                        mag[j] += m1/r1;//magnetisation(g)/r1;
                        enrgy[j] += e1/r1;//energy(g,J)/r1;
                        M = var_mag(g,m1);//magnetisation2(g);//(magnetisation(g)/r1)*(magnetisation(g)/r1);
                        E = var_enrg_hex(g,J,e1);//energy2(g,J);//(energy(g,J)/r1)*(energy(g,J)/r1);
                        spec[j] += E*t*t/r1;//fabs(E-(e1*e1))*(T*T)/r1;//specific_heat(g,J,T)/r1;
                        sus[j] += M*t/r1;//fabs(M-(m1*m1))*T/r1;//susceptibility(g,T)/r1;
                        if(k==0) temp[j] = T;
                        T+=0.1;
                        init_grid(g,c,spin,gsl_mt);
                }
        }
        }


	if(sq==1){
	for(int k=0; k<r1; k++){
		T=0.0001;
		double M=0.0;
		double E=0.0;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int i=0; i<r3; i++){
				metropolis_sweep(g,spin,J,t,gsl_mt);
			}
			double m1 = magnetisation(g);
			double e1 = energy(g,J);
			mag[j] += m1/r1;//magnetisation(g)/r1;
			enrgy[j] += e1/r1;//energy(g,J)/r1;
			M = magnetisation2(g);//(magnetisation(g)/r1)*(magnetisation(g)/r1);
			E = energy2(g,J);//(energy(g,J)/r1)*(energy(g,J)/r1);
			//printf("%f %f %f %f\n", M, m1, E, e1);
			spec[j] += var_enrg(g,J,e1)*t*t/r1;//(E-(e1*e1))*(T*T)/(r1*m*n);//specific_heat(g,J,T)/r1;
			sus[j] += var_mag(g,m1)*t/r1;//(M-(m1*m1))*T/(r1*m*n);//susceptibility(g,T)/r1;
			if(k==0) temp[j] = T;
			T+=0.1;
			init_grid(g,c,spin,gsl_mt);
		}
	}
	}
	//Save observables and free memory
	char title[100];
	if(hex==1 && c==0) snprintf(title, 100, "potts2d_stats_hexagonal_serial_hot_q%d.txt", spin);
	if(triangle==1 && c==0) snprintf(title, 100, "potts2d_stats_triangular_serial_hot_q%d.txt", spin);
	if(sq==1 && c==0) snprintf(title, 100, "potts2d_stats_square_serial_hot_q%d.txt", spin);
	if(hex==1 && c==1) snprintf(title, 100, "potts2d_stats_hexagonal_serial_cold_q%d.txt", spin);
        if(triangle==1 && c==1) snprintf(title, 100, "potts2d_stats_triangular_serial_cold_q%d.txt", spin);
        if(sq==1 && c==1) snprintf(title, 100, "potts2d_stats_square_serial_cold_q%d.txt", spin);

	write_stats(title,mag,enrgy,sus,spec,temp,r2);
	free(mag);
	free(enrgy);
	free(temp);
	free(spec);
	free(sus);
return 0;
}

//Initialise grid for random start
void init_grid(int g[m][n], int c, int s, gsl_rng *gsl_mt){
	int i,j;
	if(c==1){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				g[i][j]=s;
			}
		}
	}
	if(c==0){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				g[i][j] = gsl_rng_uniform_int(gsl_mt,s)+1;
			}
		}
	}
}
			//g[i][j]=(gsl_rng_uniform_int(gsl_mt, s))+1;
			//if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j] *= -1;
			//g[i][j]=(2*gsl_rng_uniform_int(gsl_mt, s))-s;
			/*double U = gsl_rng_uniform(gsl_mt);//drand48();
			(U<0.5) g[i][j]=-1;
			else{ 
				g[i][j]=1;
			}*/
/*
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
*/
/*
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
*/
/*
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
*/

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

void write_stats(char *title, double *mag, double *energy, double *sus, double *spec, double *T, int r2){
	FILE *fp = fopen(title, "w");
	//FILE *fp = fopen("stats_hexagonal_serial.txt", "w");
	for(int i=0; i<r2; i++){
		fprintf(fp,"%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}
/*
int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}
*/
/*
int delta(int a, int b){
if(a==b) return 1;
else{ return 0;}
}

int H(int a, int b){
return 1-delta(a,b);
}
*/
/*
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
*/
