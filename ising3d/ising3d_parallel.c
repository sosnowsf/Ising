#include"grid_size.h"
#include"modulo.h"
#include"decomp1d.h"
#include"observables.h"
#include"exchange.h"
#include"metropolis.h"
#include"mpi_io.h"

void init_grid(int g[x][y][z], int, gsl_rng *);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[x][y][z], int g[x][y][z]);
void init_grid2(int g[x][y][z], int, int);
void init_grid2_2d(int g[x][y][z], int, int, int, int);
void init_grid2_3d(int g[x][y][z], int, int, int, int, int, int);
//void init_grid2d2(int g[x][y][z], int, int, int, int);

int main(int argc, char **argv){
	int rank, size, s, e, s2, e2, s3, e3, nbrtop, nbrbot, nbrleft, nbrright, nbrfront, nbrback; //variables needed for parallelisation
	double T; //Temperature
	const double J = 1.0; //Coupling constant

        int c=0;
	int d=1;
	char *dval;
        int opt;
	//Parsing
        while((opt=getopt(argc, argv, "chd:")) != -1){
                switch(opt){
                        case 'c':
                                c=1;
                                break;
			case 'd':
				dval = optarg;
				d = atoi(dval);
				break;
                        case 'h':
                                printf("-c for cold start(hot start default)\n -d for decomposition dimensions (1,2 or 3)\n");
                                break;
                }
        }
	//Error check
	if(d!=1 && d!=2 && d!=3){
		perror("Invalid value for decomposition (can only be 1,2 or 3)");
		exit(EXIT_FAILURE);
	}

	int g[x][y][z];

	double r1 = 50.0; //Number of simulations
	int r2 = 60; //Number of temperatures
	int r3 = 1000; //Number of sweeps per simulation per temp

	//Seed prng and initiate grid
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);

	init_grid(g, c, gsl_mt);
	//int seed = 1997;
	//srand48(1997);

	//Variables for MPI Cart
	int ndims;
	int dims[1];
	int periods[1];
	int reorder;
	int source1, source2, source3, source4, source5, source6;

	//Variables for 2D decomposition
	int dims2[2];
	int periods2[2];

	int ndims2=2;
	int ddims2[2];
	ddims2[0]=ddims2[1]=0;

	//Variables for 3D decomposition
	int dims3[3];
	int periods3[3];

	int ndims3=3;
	int ddims3[3];
	ddims3[0]=ddims3[1]=ddims3[2]=0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	//1D decomposition
	if(d==1){
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

	//printf("%d %d %d\n", rank, s, e);

	//Change prng for each process 
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);

	double *mag, *eng, *specs, *sus, *temps;
	if(rank==0){
		mag = (double *)malloc(r2*sizeof(double));
		eng = (double *)malloc(r2*sizeof(double));
		temps = (double *)malloc(r2*sizeof(double));
		specs = (double *)malloc(r2*sizeof(double));
		sus = (double *)malloc(r2*sizeof(double));
	}

	for(int i=0; i<r1; i++){
		T=0.0001;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,s,e,J,T,gsl_mt);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
			}
			//calculate the observables at each rank and reduce to root process
			double M,E,m2,e2;
			double M1 = magnetisation(g,s,e);
			double E1 = energy(g,J,s,e);
			double M2 = var_mag(g,s,e,M1);
			double E2 = var_enrg(g,s,e,J,E1);
			MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &m2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				M/=size;
				E/=size;
				m2/=size;
				e2/=size;
				mag[j]+=M/r1;
				eng[j]+=E/r1;
				specs[j] += e2*t*t/r1;
				sus[j] += m2*t/r1;
				if(i==0) temps[j] = T;
			}
		T+=0.1;
		}
	init_grid(g,c,gsl_mt);
	exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
	}

        //Test MPI-IO
/*	init_grid2(g,s,e);
	printf("%d %f\n", rank, magnetisation(g,s,e));
	save_grid(g,s,e);
	init_grid(g,1,gsl_mt);
	printf("%d %f\n", rank, magnetisation(g,s,e));
	MPI_Barrier(MPI_COMM_WORLD);
	load_grid(g,s,e);
	printf("%d %f\n", rank, magnetisation(g,s,e));
*/
	

	if(rank==0){
		char title[100];
		if(c==0) snprintf(title, 100, "stats_cube_parallel_hot.txt");
		if(c==1) snprintf(title, 100, "stats_cube_parallel_cold.txt");
		write_stats(title,mag,eng,specs,sus,temps,r2);
	}
	}


	//2D decomposition
	if(d==2){
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

	//Change prng for each process
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);
	//srand48(seed*rank+1);
	

	//Allocate memory to store the required stats
	double *mag, *eng, *specs, *sus, *temps;
	if(rank==0){
		mag = (double *)malloc(r2*sizeof(double));
		eng = (double *)malloc(r2*sizeof(double));
		temps = (double *)malloc(r2*sizeof(double));
		specs = (double *)malloc(r2*sizeof(double));
		sus = (double *)malloc(r2*sizeof(double));
	}		

        for(int i=0; i<r1; i++){
                T=0.0001;
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int k=0; k<r3; k++){
                                metropolis_sweep2d(g,s,e,s2,e2,J,T,gsl_mt);
                                exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                        }
                        double M,E,m2,ee2;
                        double M1 = magnetisation2d(g,s,e,s2,e2);
                        double E1 = energy2d(g,J,s,e,s2,e2);
                        double M2 = var_mag2d(g,s,e,s2,e2,M1);
                        double E2 = var_enrg2d(g,s,e,s2,e2,J,E1);
                        MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&M2, &m2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E2, &ee2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                M/=size;
                                E/=size;
                                m2/=size;
                                ee2/=size;
                                mag[j]+=M/r1;
                                eng[j]+=E/r1;
                                specs[j] += ee2*t*t/r1;
                                sus[j] += m2*t/r1;
                                if(i==0) temps[j] = T;
                        }
                T+=0.1;
                }
        init_grid(g,c,gsl_mt);
        exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
        }

        if(rank==0){
                char title[100];
                if(c==0) snprintf(title, 100, "stats_cube_parallel_hot_2d.txt");
                if(c==1) snprintf(title, 100, "stats_cube_parallel_cold_2d.txt");
                write_stats(title,mag,eng,specs,sus,temps,r2);
        }

        //Test MPI-IO
/*        init_grid2_2d(g,s,e,s2,e2);
        printf("%d %f\n", rank, magnetisation2d(g,s,e,s2,e2));
        save_grid2d(g,s,e,s2,e2);
        init_grid(g,1,gsl_mt);
        printf("%d %f\n", rank, magnetisation2d(g,s,e,s2,e2));
        MPI_Barrier(MPI_COMM_WORLD);
        load_grid2d(g,s,e,s2,e2);
        printf("%d %f\n", rank, magnetisation2d(g,s,e,s2,e2));
*/

	}


	//3D Decomposition
	if(d==3){

        MPI_Comm cart;
        ndims=1;
        dims3[0]=size;
        dims3[1]=size;
	dims3[2]=size;
        periods3[0]=1;
        periods3[1]=1;
	periods3[2]=2;
        reorder=0;
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims3, periods3, reorder, &cart);

        //Find the neighbouring processes
        MPI_Dims_create(size, ndims3, ddims3);

        MPI_Cart_shift(cart, 0, -1, &source1, &nbrleft);
        MPI_Cart_shift(cart, 0, 1, &source2, &nbrright);
        MPI_Cart_shift(cart, 0, -ddims3[0], &source3, &nbrtop);
        MPI_Cart_shift(cart, 0, ddims3[0], &source4, &nbrbot);
	MPI_Cart_shift(cart, 0, ddims3[1]*ddims3[2], &source5, &nbrfront);
	MPI_Cart_shift(cart, 0, -ddims3[1]*ddims3[2], &source6, &nbrback);


	if(rank==0) printf("%d %d %d\n", ddims3[0], ddims3[1], ddims3[2]);

	//Find correct left and right neighbours
        if((rank%ddims3[1])==0) nbrleft = (nbrleft+ddims3[1])%size;
        if((rank%ddims3[1])==(ddims3[1]-1)) nbrright = modulo(nbrright-ddims3[1],size);

	//Find the correct top and bottom neighbours
	int r = (rank - (rank%(ddims3[1]*ddims3[2]))); //Position along the dimension decomposition ddims3[0]

	nbrtop = rank-ddims3[2];
	nbrbot = rank+ddims3[2];

	if(nbrtop < r) nbrtop = modulo(nbrtop+ddims3[1]*ddims3[2],size);
	if(nbrbot >= (r+(ddims3[1]*ddims3[2])) ) nbrbot -= ddims3[1]*ddims3[2];

        //printf("%d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot);

	
	decomp1d(z, ddims3[2], rank%ddims3[2], &s, &e); //Decompose along y axis
	decomp1d(y, ddims3[1], ((rank-rank%ddims3[2])/ddims3[2])%ddims3[2], &s2, &e2); //Decompose along z axis
	decomp1d(z, ddims3[0], ((rank-rank%(ddims3[1]*ddims3[2]))/(ddims3[1]*ddims3[2])), &s3, &e3); //Decompose along x axis

	//printf("%d\n", ((rank-rank%(ddims3[1]*ddims3[2]))/(ddims3[1]*ddims3[2])));


        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, nbrfront, nbrback, s, e, s2, e2, s3, e3);

        //Change prng for each process
        unsigned long seed2 = (1999*rank);
        gsl_rng_set(gsl_mt, seed2);
        //srand48(seed*rank+1);

        //Allocate memory to store the required stats
        double *mag, *eng, *specs, *sus, *temps;
        if(rank==0){
                mag = (double *)malloc(r2*sizeof(double));
                eng = (double *)malloc(r2*sizeof(double));
                temps = (double *)malloc(r2*sizeof(double));
                specs = (double *)malloc(r2*sizeof(double));
                sus = (double *)malloc(r2*sizeof(double));
        }
	
        for(int i=0; i<r1; i++){
                T=0.0001;
                double t;
                for(int j=0; j<r2; j++){
                        t=1/T;
                        for(int k=0; k<r3; k++){
                                metropolis_sweep3d(g,s,e,s2,e2,s3,e3,J,T,gsl_mt);
                                exchange3d(g,s,e,s2,e2,s3,e3,nbrtop,nbrbot,nbrleft,nbrright,nbrfront,nbrback,MPI_COMM_WORLD);
                        }
                        double M,E,m2,ee2;
                        double M1 = magnetisation3d(g,s,e,s2,e2,s3,e3);
                        double E1 = energy3d(g,J,s,e,s2,e2,s3,e3);
                        double M2 = var_mag3d(g,s,e,s2,e2,s3,e3,M1);
                        double E2 = var_enrg3d(g,s,e,s2,e2,s3,e3,J,E1);
                        MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&M2, &m2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E2, &ee2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                M/=size;
                                E/=size;
                                m2/=size;
                                ee2/=size;
                                mag[j]+=M/r1;
                                eng[j]+=E/r1;
                                specs[j] += ee2*t*t/r1;
                                sus[j] += m2*t/r1;
                                if(i==0) temps[j] = T;
                        }
                T+=0.1;
                }
        init_grid(g,c,gsl_mt);
        exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
        }

        if(rank==0){
                char title[100];
                if(c==0) snprintf(title, 100, "stats_cube_parallel_hot_3d.txt");
                if(c==1) snprintf(title, 100, "stats_cube_parallel_cold_3d.txt");
                write_stats(title,mag,eng,specs,sus,temps,r2);
        }
/*
        //Test MPI-IO
        init_grid2_3d(g,s,e,s2,e2,s3,e3);
        printf("%d %f\n", rank, magnetisation3d(g,s,e,s2,e2,s3,e3));
        save_grid3d(g,s,e,s2,e2,s3,e3);
        init_grid(g,1,gsl_mt);
        printf("%d %f\n", rank, magnetisation3d(g,s,e,s2,e2,s3,e3));
        MPI_Barrier(MPI_COMM_WORLD);
        load_grid3d(g,s,e,s2,e2,s3,e3);
        printf("%d %f\n", rank, magnetisation3d(g,s,e,s2,e2,s3,e3));
*/

	}


	MPI_Finalize();
return 0;
}

//Initiate grid for random start
void init_grid(int g[x][y][z], int c, gsl_rng *gsl_mt){
	int i,j,k;
	if(c==1){
		for(i=0;i<x;i++){
			for(j=0;j<y;j++){
				for(k=0;k<z;k++){
					g[i][j][k]= 1;
				}
			}
		}
	}
	if(c==0){
		for(i=0; i<x; i++){
			for(j=0; j<y; j++){
				for(k=0; k<z; k++){
					if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j][k]=-1;
					else{ g[i][j][k]=1; }
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


void init_grid2(int g[x][y][z], int s, int e){
	int i,j,k;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				g[i][j][k]= i*y*z + j*z + k;
			}
		}
	}
}

void init_grid2_2d(int g[x][y][z], int s, int e, int s2, int e2){
        int i,j,k;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                g[i][j][k]= i*y*z + j*z + k;
                        }
                }
        }
}

void init_grid2_3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3){
        int i,j,k;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                g[i][j][k]= i*y*z + j*z + k;
                        }
                }
        }
}

