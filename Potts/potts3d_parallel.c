#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include"decomp1d.h"

#define x 10
#define y 10
#define z 10

void init_grid(int g[x][y][z], int, gsl_rng *);
void metropolis_sweep(int g[x][y][z], int, int, const double , double, gsl_rng *);
double energy(int g[z][y][z], const double, int, int);
double magnetisation(int g[x][y][z], int, int);
void exchange(int g[x][y][z], int, int, int, int, MPI_Comm comm);
void write_stats(char *, double *, double *, double *, double *, double *, int);
int modulo(int, int);
void reset_grid(int G[x][y][z], int g[x][y][z]);
int delta(int, int);
int H(int, int);

int main(int argc, char **argv){
	int rank, size, s, e, nbrtop, nbrbot; //variables needed for parallelisation
	double T; //Temperature
	const double J = 1.0; //Coupling constant

	int g[x][y][z];

	double r1 = 1.0; //Number of simulations
	int r2 = 100; //Number of temperatures
	int r3 = 1000; //Number of sweeps per simulation per temp

	//Seed prng and initiate grid
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);

	init_grid(g, gsl_mt);
	//int seed = 1997;
	//srand48(1997);

	//Variables for MPI Cart
	int ndims;
	int dims[1];
	int periods[1];
	int reorder;
	int source1, source2;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm cart;
	ndims=1;
	dims[0]=size;
	periods[0]=1;
	reorder=0;

	//Create MPI_Cart for finding neighbouring processes
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart);
	MPI_Cart_shift(cart, 0, -1, &source1, &nbrtop);
	MPI_Cart_shift(cart, 0, 1, &source2, &nbrbot);

	//Find the size of each subarray
	decomp1d(x, size, rank, &s, &e);

	//Change prng for each process 
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);

	double *mag, *eng, *specs, *sus, *temps;
	if(rank==0){
		mag = (double *)malloc(r2*sizeof(double));
		eng = (double *)malloc(r2*sizeof(double));
		temps = (double *)malloc(r2*sizeof(double));
		specs = malloc(r2*sizeof(double));
		sus = malloc(r2*sizeof(double));
	}


	for(int i=0; i<r1; i++){
		T=0.0001;
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,s,e,J,T,gsl_mt);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
			}
			double M,E;
			double M1 = magnetisation(g,s,e);
			double E1 = energy(g,J,s,e);
			MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				M/=size;
				E/=size;
				//printf("%f %f %f\n", mag, enrg, T);
				mag[j]+=M/r1;
				eng[j]+=E/r1;
				specs[j] += ((E/r1)*(E/r1) - (eng[j])*(eng[j]))*(T*T)/r1;
				sus[j] += ((M/r1)*(M/r1) - (mag[j])*(mag[j]))*T/r1;
				if(i==0) temps[j] = T;
			}
		T+=0.05;
		}
	init_grid(g,gsl_mt);
	//exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
	}

	if(rank==0){
		char title[100];
		snprintf(title, 100, "potts3d_stats_cube_parallel.txt");
		write_stats(title,mag,eng,specs,sus,temps,r2);
	}
	MPI_Finalize();
return 0;
}

//Initiate grid for random start
void init_grid(int g[x][y][z], int s, gsl_rng *gsl_mt){
	int i,j,k;
	for(i=0;i<x;i++){
		for(j=0;j<y;j++){
			for(k=0;k<z;k++){
				g[i][j][k]= s;//2*gsl_rng_uniform_int(gsl_mt, 2)-1;
				/*if(drand48()<0.5) g[i][j][k] = -1;
				else{
					g[i][j][k] = 1;
				}*/
			}
		}
	}
}

void metropolis_sweep(int g[x][y][z], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
	int i,j,k,r;
	double dE;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				r=gsl_rng_uniform_int(gsl_mt,spin)+1;
				dE = 2.0*J*(H(r,g[modulo(i+1,x)][j][k]) + H(r,g[modulo(i-1,x)][j][k]) + H(r,g[i][modulo(j+1,y)][k]) + H(r,g[i][modulo(j-1,y)][k]) + H(r,g[i][j][modulo(k+1,z)]) + H(r,g[i][j][modulo(k-1,z)]));
				if(dE<=0) g[i][j][k]=r;
				else{
					if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))) g[i][j][k]=r;
				}	
			}
		}
	}
}

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

double energy(int g[x][y][z], const double J, int s, int e){
	int i,j,k;
	double E=0;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				E+=-J*(H(g[i][j][k],g[modulo(i+1,x)][j][k]) + H(g[i][j][k],g[modulo(i-1,x)][j][k]) + H(g[i][j][k],g[i][modulo(j+1,y)][k]) + H(g[i][j][k],g[i][modulo(j-1,y)][k]) + H(g[i][j][k],g[i][j][modulo(k+1,z)]) + H(g[i][j][k],g[i][j][modulo(k-1,z)]));
			}
		}
	}
return E/(2.0*(e-s+1)*y*z);
}

void exchange(int g[x][y][z], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
	MPI_Request reqs[4];
	MPI_Datatype slice;
	MPI_Type_vector(y*z, 1, 1, MPI_INT, &slice);
	MPI_Type_commit(&slice);

	MPI_Irecv(&g[modulo(s-1,x)][0][0], 1, slice, nbrtop, 0, comm, &reqs[0]);
	MPI_Irecv(&g[modulo(e+1,x)][0][0], 1, slice, nbrbot, 0, comm, &reqs[1]);
	MPI_Isend(&g[modulo(e,x)][0][0], 1, slice, nbrbot, 0, comm, &reqs[2]);
	MPI_Isend(&g[modulo(s,x)][0][0], 1, slice, nbrtop, 0, comm, &reqs[3]);


//	MPI_Irecv(&g[modulo(s-1,x)][0][0], 100, MPI_INT, nbrtop, 0, comm, &reqs[0]);
//	MPI_Irecv(&g[modulo(e+1,x)][0][0], 100, MPI_INT, nbrbot, 0, comm, &reqs[1]);
//	MPI_Isend(&g[modulo(e,x)][0][0], 100, MPI_INT, nbrbot, 0, comm, &reqs[2]);
//	MPI_Isend(&g[modulo(s,x)][0][0], 100, MPI_INT, nbrtop, 0, comm, &reqs[3]);

	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&slice);
}

void write_stats(char *title, double *mag, double *energy, double *specs, double *sus, double *T, int r1){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r1; i++){
		fprintf(fp,"%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], specs[i]);
	}
	fclose(fp);
}

//Reset grid after every temp increase
void reset_grid(int G[x][y][z], int g[x][y][z]){
	for(int i=0; i<x; i++){
		for(int j=0; j<y; j++){
			for(int k=0; k<z; k++){
				g[i][j][k]=G[i][j][k];
			}
		}
	}
}
int modulo(int a, int b){
	int r = a%b;
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
