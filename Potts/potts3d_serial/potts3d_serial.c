#include"grid_size.h"
#include"metropolis.h"
#include"observables.h"
#include"modulo.h"
#include"hamiltonian.h"

void init_grid(int g[x][y][z], int, int, gsl_rng *);
void write_stats(char *, double *, double *, double *, double *, double *, int);

int main(int argc, char **argv){
	double T; //Temperature
	const double J = 1.0; //Coupling Temperature
	int spin=2;	
	double r1 = 50.0; //Number of simulations
	int r2 = 100; //Number of tempreatures;
	int r3 = 1000; //Number of sweeps
	int c = 0;
	int opt;
	while((opt = getopt(argc, argv, "hc")) != -1){
		switch(opt){
			case 'c':
				c=1;
				break;
			case 'h':
				printf("Use -c for cold start (hot start default)");
				break;
		}
	}


	//Seed rng
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
        unsigned long seed = 1999;
        gsl_rng_set(gsl_mt, seed);

	//int seed = 1997;
	//srand48(seed);

	//Initialise grid
	int g[x][y][z];
	init_grid(g,c,spin,gsl_mt);

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
		//double M=0.0;
		//double E=0.0;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,spin,J,T,gsl_mt);
			}
			double m = magnetisation(g);
			double e = energy(g,J);
			double m2 = var_mag(g,m);
			double e2 = var_enrg(g,J,e);
			mag[j] += m/r1;
			enrgy[j] += e/r1;
			//M = magnetisation2(g)/r1;
			//E = energy2(g,J)/r1;
			spec[j] += e2*t*t/r1;//fabs(E-(e*e))*(T*T)/r1;
			sus[j] += m2*t/r1;//fabs(M-(m*m))*T/r1;
			if(i==0) temp[j]=T;
			T+=0.06;
			init_grid(g,c,spin,gsl_mt);
		}
	}
	char title[100];
	if(c==0) snprintf(title, 100, "potts3d_stats_cube_hot_q%d.txt", spin);
	if(c==1) snprintf(title, 100, "potts3d_stats_cube_cold_q%d.txt", spin);
	write_stats(title,mag,enrgy,spec,sus,temp,r2);
	free(mag);
	free(enrgy);
	free(spec);
	free(sus);
	free(temp);
return 0;
}

void init_grid(int g[x][y][z], int c, int s, gsl_rng *gsl_mt){
	int i,j,k;
	if(c==1){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					g[i][j][k]=s;
				/*if(drand48()<0.5) g[i][j][k]=-1;
				else{
					g[i][j][k]=1;
				}*/
				}
			}
		}
	}
	if(c==0){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					g[i][j][k]=gsl_rng_uniform_int(gsl_mt,s)+1;
				}
			}
		}
	}
}

void write_stats(char *title, double *mag, double *energy, double *sus, double *spec, double *T, int r2){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r2; i++){
		fprintf(fp, "%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}
