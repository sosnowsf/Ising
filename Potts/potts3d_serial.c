#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define x 10
#define y 10
#define z 10

void init_grid(int g[x][y][z], int, gsl_rng *);
void metropolis_sweep(int g[x][y][z], int, const double, double, gsl_rng *);
double energy(int g[x][y][z], const double);
double magnetisation(int g[x][y][z]);
void write_stats(char *, double *, double *, double *, double *, double *, int);
int modulo(int, int);
int delta(int, int);
int H(int, int);

int main(){
	double T; //Temperature
	const double J = 1.0; //Coupling Temperature
	
	double r1 = 1.0; //Number of simulations
	int r2 = 100; //Number of tempreatures;
	int r3 = 1000; //Number of sweeps

	//Seed rng
	int seed = 1997;
	srand48(seed);

	//Initialise grid
	int g[x][y][z];
	init_grid(g);

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
		double M=0.0;
		double E=0.0;
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,J,T);
			}
			double m = magnetisation(g)/r1;
			double e = energy(g,J)/r1;
			mag[j] += m;
			enrgy[j] += e;
			M += m*m;
			E += e*e;
			spec[j] += (E-enrgy[j]*enrgy[j])*(T*T)/r1;
			sus[j] += (M-mag[j]*mag[j])*T/r1;
			if(i==0) temp[j]=T;
			T+=0.05;
			init_grid(g);
		}
	}
	char title[100];
	snprintf(title, 100, "potts3d_stats_cube.txt");
	write_stats(title,mag,enrgy,spec,sus,temp,r2);
	free(mag);
	free(enrgy);
	free(spec);
	free(sus);
	free(temp);
return 0;
}

void init_grid(int g[x][y][z], int s, gsl_rng *gsl_mt){
	int i,j,k;
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

void metropolis_sweep(int g[x][y][z], int s, const double J, double T, gsl_rng *gsl_mt){
	int i,j,k,r;
	double dE;
	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				r=gsl_rng_uniform(gsl_mt,s)+1;
				dE = 2.0*J*(H(r,g[modulo(i+1,x)][j][k]) + H(r,g[modulo(i-1,x)][j][k]) + H(r,g[i][modulo(j+1,y)][k]) + H(r,g[i][modulo(j-1,y)][k]) + H(r,g[i][j][modulo(k+1,z)]) + H(r,g[i][j][modulo(k-1,z)]));
				if(dE<=0) g[i][j][k] = r;
				else{
					if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))) g[i][j][k] = r;
				}
			}
		}
	}
}

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

void write_stats(char *title, double *mag, double *energy, double *sus, double *spec, double *T, int r2){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r2; i++){
		fprintf(fp, "%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}

int modulo(int a, int b){
	int r=a%b;
	if(r<0) r+=b;
return r;
}

int delta(int a, int b){
if(a==b) return 1;
else{ return 0;}
}

int H(int a, int b){
return 1-delta(a,b);
}
