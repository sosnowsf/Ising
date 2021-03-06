#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define x 20
#define y 20
#define z 20

void init_grid(int g[x][y][z], int, gsl_rng *);
void metropolis_sweep(int g[x][y][z], const double, double, gsl_rng *);
double energy(int g[x][y][z], const double);
double energy_var(int g[x][y][z], const double, double);
double energy2(int g[x][y][z], const double);
double magnetisation(int g[x][y][z]);
double mag_var(int g[x][y][z], double);
double magnetisation2(int g[x][y][z]);
void write_stats(char *, double *, double *, double *, double *, double *, int);
int modulo(int, int);

int main(int argc, char *argv[]){
	double T; //Temperature
	const double J = 1.0; //Coupling Temperature
	int c=0;
	int opt;
	while((opt=getopt(argc, argv, "ch")) != -1){
		switch(opt){
			case 'c':
				c=1;
				break;
			case 'h':
				printf("-c for cold start(hot start default)");
				break;
		}
	}

	double r1 = 50.0; //Number of simulations
	int r2 = 60; //Number of tempreatures;
	int r3 = 1000; //Number of sweeps

	//Seed rng
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);
	//int seed = 1997;
	//srand48(seed);

	//Initialise grid
	int g[x][y][z];
	init_grid(g,c,gsl_mt);

	//Allocate memory for observables
	double *mag, *enrgy, *spec, *sus, *temp;
	mag = malloc(r2*sizeof(double));
	enrgy = malloc(r2*sizeof(double));
	spec = malloc(r2*sizeof(double));
	sus = malloc(r2*sizeof(double));
	temp = malloc(r2*(sizeof(double)));

	for(int i=0; i<r1; i++){
		//Reset variables;
		T=0.00001;
		double m2;
		double e2;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,J,t,gsl_mt);
			}
			double m = magnetisation(g);
			double e = energy(g,J);
			m2 = mag_var(g,m);
			e2 = energy_var(g,J,e);
			mag[j] += m/r1;
			enrgy[j] += e/r1;
			spec[j] += e2*t*t/r1;
			sus[j] += m2*t/r1;
			if(i==0) temp[j]=T;
			T+=0.1;
			init_grid(g,c,gsl_mt);
		}
	}
	char title[100];
	if(c==0) snprintf(title, 100, "stats_cube_hot.txt");
	if(c==1) snprintf(title, 100, "stats_cube_cold.txt");
	write_stats(title,mag,enrgy,sus,spec,temp,r2);
	free(mag);
	free(enrgy);
	free(spec);
	free(sus);
	free(temp);
return 0;
}

//Initialise grid
void init_grid(int g[x][y][z], int c, gsl_rng *gsl_mt){
	int i,j,k;
	if(c==1){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					g[i][j][k]=1;
				}
			}
		}
	}
	if(c==0){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j][k]=-1;
					else{ g[i][j][k]=1;}
				}
			}
		}
	}
}

//Cubic metropolis sweep
void metropolis_sweep(int g[x][y][z], const double J, double T, gsl_rng *gsl_mt){
	int i,j,k;
	double dE;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				dE = 2*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
				if(dE<=0) g[i][j][k] *= -1;
				else{
					if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j][k]*=-1;
				}
			}
		}
	}
}

//average energy
double energy(int g[x][y][z], const double J){
	int i,j,k;
	double E=0;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				E += -J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
			}
		}
	}
return E/(2.0*x*y*z);
}

//average energy variance
double energy_var(int g[x][y][z], const double J, double avg){
	int i,j,k;
	double var=0;
	double dE;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				dE = -J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
				var += (dE-avg)*(dE-avg);
			}
		}
	}
return var/(x*y*z);
}

//average magnetisation
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
double mag_var(int g[x][y][z], double avg){
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

void write_stats(char *title, double *mag, double *energy, double *sus, double *spec, double *T, int r2){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r2; i++){
		fprintf(fp, "%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}

//modulo function
int modulo(int a, int b){
	int r=a%b;
	if(r<0) r+=b;
return r;
}
