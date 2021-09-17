#include"grid_size.h"
#include"modulo.h"
#include"hamiltonian.h"
#include"exchange.h"
#include"metropolis.h"
#include"observables.h"
#include"mpi_io.h"

void init_grid(int g[x][y][z], int, int, gsl_rng *);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[x][y][z], int g[x][y][z]);

int main(int argc, char **argv){
	int rank, size, s, e, s2, e2, nbrtop, nbrbot, nbrleft, nbrright; //variables needed for parallelisation
	double T; //Temperature
	const double J = 1.0; //Coupling constant
	int spin = 2;

	int g[x][y][z];

	double r1 = 10.0; //Number of simulations
	int r2 = 50; //Number of temperatures
	int r3 = 1000; //Number of sweeps per simulation per temp
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


	//Seed prng and initiate grid
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);

	init_grid(g, c, spin, gsl_mt);
	//int seed = 1997;
	//srand48(1997);

	//Variables for MPI Cart
	int ndims;
	//int dims[1];
	//int periods[1];
	int reorder;
	int source1, source2, source3, source4;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//Variables for 2D decomposition
	int dims2[2];
	int periods2[2];

	int ndims2=2;
	int ddims2[2];
	ddims2[0]=ddims2[1]=0;


	MPI_Comm cart;
	ndims=1;
	dims2[0]=size;
	dims2[1]=size;
	periods2[0]=1;
	periods2[1]=1;
	reorder=0;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims2, periods2, reorder, &cart);

	//Find the neighbouring processes
	MPI_Dims_create(size, ndims2, ddims2);

	MPI_Cart_shift(cart, 0, -1, &source1, &nbrleft);
	MPI_Cart_shift(cart, 0, 1, &source2, &nbrright);
	MPI_Cart_shift(cart, 0, -ddims2[0], &source3, &nbrtop);
	MPI_Cart_shift(cart, 0, ddims2[0], &source4, &nbrbot);

	if((rank%ddims2[0])==0) nbrleft = (nbrleft+ddims2[0])%size;
	if((rank%ddims2[0])==(ddims2[0]-1)) nbrright = modulo(nbrright-ddims2[0],size);

	decomp1d(y, ddims2[0], rank%ddims2[0], &s2, &e2); //Decompose cols
	decomp1d(x, ddims2[1], (rank-(rank%ddims2[0]))/ddims2[0], &s, &e); //Decompose rows

	//printf("%d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, s, e, s2, e2);

	//srand48(seed*rank+1);

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
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int k=0; k<r3; k++){
				metropolis_sweep2d(g,spin,s,e,s2,e2,J,T,gsl_mt);
				exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
			}
			double M,E,M2,E2;
			double M1 = magnetisation2d(g,s,e,s2,e2);
			double E1 = energy2d(g,J,s,e,s2,e2);
			double M22 = var_mag2d(g,s,e,s2,e2,M1);
			double E22 = var_enrg2d(g,s,e,s2,e2,J,E1);
			MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M22, &M2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E22, &E2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				M/=size;
				E/=size;
				M2/=size;
				E2/=size;
				//printf("%f %f %f\n", mag, enrg, T);
				mag[j]+=M/r1;
				eng[j]+=E/r1;
				specs[j] += E2*t*t/r1;
				sus[j] += M2*t/r1;
				if(i==0) temps[j] = T;
			}
		T+=0.1;
		}
	init_grid(g,c,spin,gsl_mt);
	exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
	}

	if(rank==0){
		char title[100];
		if(c==0) snprintf(title, 100, "potts3d_stats_cube_parallel_hot_2d.txt");
		if(c==1) snprintf(title, 100, "potts3d_stats_cube_parallel_cold_2d.txt");
		write_stats(title,mag,eng,specs,sus,temps,r2);
	}
	MPI_Finalize();
return 0;
}

//Initiate grid for random start
void init_grid(int g[x][y][z], int c, int s, gsl_rng *gsl_mt){
	int i,j,k;
	if(c==1){
		for(i=0;i<x;i++){
			for(j=0;j<y;j++){
				for(k=0;k<z;k++){
					g[i][j][k] = s;//2*gsl_rng_uniform_int(gsl_mt, 2)-1;
					/*if(drand48()<0.5) g[i][j][k] = -1;
					else{
						g[i][j][k] = 1;
					}*/
				}
			}
		}
	}
	if(c==0){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					g[i][j][k] = gsl_rng_uniform_int(gsl_mt,s)+1;
				}
			}
		}
	}
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
