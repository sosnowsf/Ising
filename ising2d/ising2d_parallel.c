#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include"decomp1d.h"

#define m 100
#define n 100

void init_grid(int g[m][n], int, gsl_rng *);
void metropolis_sweep(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep_triangular2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
double energy(int g[m][n], const double, int, int);
double energy2d(int g[m][n], const double, int, int, int, int);
double energy2(int g[m][n], const double, int, int);
double energy_triangular(int g[m][n], const double, int, int);
double energy_triangular2d(int g[m][n], const double, int, int, int, int);
double energy_hexagonal(int g[m][n], const double, int, int);
double energy_hexagonal2d(int g[m][n], const double, int, int, int, int);
double magnetisation(int g[m][n], int, int);
double magnetisation2d(int g[m][n], int, int, int, int);
double magnetisation2(int g[m][n], int, int);
void print_matrix(int g[m][n]);
void exchange(int g[m][n], int, int, int, int, MPI_Comm);
void exchange2(int g[m][n], int, int, int, int, int, int, MPI_Comm);
void exchange2d(int g[m][n], int, int, int, int, int, int, int, int, MPI_Comm);
//void gather(int g[m][n], int, int, int, int);
void print_in_order(int g[m][n], MPI_Comm comm);
int modulo(int, int);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
double var_mag(int g[m][n], int, int, double);
double var_mag2d(int g[m][n], int, int, int, int, double);
double var_enrg(int g[m][n], int, int, double, double);
double var_enrg2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_tri(int g[m][n], int, int, double, double);
double var_enrg_tri2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_hex(int g[m][n], int, int, double, double);
double var_enrg_hex2d(int g[m][n], int, int, int, int, double, double);
void save_grid(int g[m][n], int, int);
void save_grid2d(int g[m][n], int, int, int, int);
void load_grid(int g[m][n], int, int);
void load_grid2d(int g[m][n], int, int, int, int);
void init_grid2(int g[m][n], int, int);
void init_grid2d2(int g[m][n], int, int, int, int);

int main(int argc, char **argv){
	int rank, size, s, e, s2, e2, nbrtop, nbrbot, nbrleft, nbrright;
	double T = 0.0001; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	//const double k = pow(1.38064852, -23); //Boltzmann constant
	int c=0;

	double r1 = 10.0;//50.0; //Number of simulations
	int r2 = 100; //Number of temperatures
	int r3 = 10000; //Number of sweeps

	//Seed prng and inititate grid
	//Seed prng and initiate grid
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);
	//int seed = 1997;
	//srand48(seed);

	int g[m][n]; //define grid
	//int G[m][n];
	//init_grid(g,c,gsl_mt);
	//reset_grid(G,g);


        //Parse variables
        int sq=0;
        int triangle=0;
        int hex=0;
	int d=1;
        int opt;
        while((opt = getopt(argc, argv, "123hcd")) != -1){
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
			case 'c':
				c=1;
				break;
			case 'd':
				d=2;
				break;
			case 'h':
				printf("Command lines:\n -1 for square grid\n -2 for triangular grid \n -3 for hexagonal grid\n -c for cold start (hot start default)\n -d for two dimensional exchange");
				return 0;
                }
        }
	init_grid(g,c,gsl_mt);
	//Error checking
        if((sq+triangle+hex)==0){
                perror("No parse command given, use -h for help");
                _exit(EXIT_FAILURE);
        }
        if((sq+triangle+hex)>1){
                perror("Please provide only one parse command at a time");
                _exit(EXIT_FAILURE);
        }

	//Variables for MPI Cart
	int ndims;
	int dims[1];
	int periods[1];
	int reorder;
	int source1, source2, source3, source4;

	//Variables for 2D decomposition
	int dims2[2];
	int periods2[2];

	int ndims2=2;
	int ddims2[2];
	ddims2[0]=ddims2[1]=0;
	
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

	//Create new communicator for finding neighbours
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart);
	MPI_Cart_shift(cart, 0, -1, &source1, &nbrtop);
	MPI_Cart_shift(cart, 0, 1, &source2, &nbrbot);

	//Find size of each subarray
	decomp1d(n, size, rank, &s, &e);
	
	//Change prng for each process
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);
	//srand48(seed*rank+1);
	
	//printf("%d %d %d %d %d\n", rank, s, e, nbrtop, nbrbot);

	//Allocate memory to store the required stats
	double *magnetisations, *energies, *specs, *sus, *temps;
	if(rank==0){
		magnetisations = (double *)malloc(r2*sizeof(double));
		energies = (double *)malloc(r2*sizeof(double));
		temps = (double *)malloc(r2*sizeof(double));
		specs = malloc(r2*sizeof(double));
		sus = malloc(r2*sizeof(double));
	}

	//Perform the simulations and calculate an average for each observables at each temperature
	if(sq==1){
	for(int k=0; k<r1; k++){
		T=0.0001; //reset the temperature for each simulation
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int i=0; i<r3; i++){
				metropolis_sweep(g,s,e,J,t,gsl_mt);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
				/*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);*/
				//if(rank%2==i%2) metropolis_sweep(g,s,e,J,T,k);
				//exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
			}
			//calculate the magnetisation per site for each process then get the average on the root process
			double mag, enrg, mag2, enrg2;
			double M = magnetisation(g,s,e);
			double E = energy(g,J,s,e);
			double M2 = var_mag(g,s,e,M);//magnetisation2(g,s,e);
			double E2 = var_enrg(g,s,e,J,E);//energy2(g,J,s,e);
			//printf("%d %f %f %f %f\n", rank, M, E, M2, E2);
			MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				mag/=size;
				enrg/=size;
				mag2/=size;
				enrg2/=size;
				magnetisations[j]+=mag/r1;
				energies[j]+=enrg/r1;
				specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
				sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
				if(k==0) temps[j] = T;
			}
			T+=0.05;
			//Save configuration
			//if(k==0 && j==0) save_grid(g,s,e);

			//Reset grid to original arrangement
			//reset_grid(G,g);
			init_grid(g,c,gsl_mt);
			exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
		}
	}
	}

	if(triangle==1){
        for(int k=0; k<r1; k++){
                T=0.0001; //reset the temperature for each simulation
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_triangular(g,s,e,J,t,gsl_mt);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                /*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);*/
                                //if(rank%2==i%2) metropolis_sweep_triangular(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the magnetisation per site for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation(g,s,e);
                        double E = energy_triangular(g,J,s,e);
			double M2 = var_mag(g,s,e,M);//magnetisation2(g,s,e);
			double E2 = var_enrg_tri(g,s,e,J,E);//energy_triangular2(g,J,s,e);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
                                sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.07;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,gsl_mt);
			exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                }
        }
	}
	if(hex==1){
        for(int k=0; k<r1; k++){
                T=0.0001; //reset the temperature for each simulation
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_hexagonal(g,s,e,J,t,gsl_mt);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                /*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);*/
                                //if(rank%2==i%2) metropolis_sweep_hexagonal(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the magnetisation per site for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation(g,s,e);
                        double E = energy_hexagonal(g,J,s,e);
			double M2 = var_mag(g,s,e,M);//magnetisation2(g,s,e);
			double E2 = var_enrg_hex(g,s,e,J,E);//energy_hexagonal2(g,J,s,e);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
                                sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.05;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,gsl_mt);
			exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                }
        }
	}
	//Test MPI-IO
/*	print_matrix(g);
	save_grid(g,s,e);
	MPI_Barrier(MPI_COMM_WORLD);
	init_grid2(g,0,n-1);//s,e);
	print_matrix(g);
	load_grid(g,s,e);
	MPI_Barrier(MPI_COMM_WORLD);
	print_matrix(g);
*/
    //    init_grid2d2(g,s,e,0,n);//0,n-1, 0, m-1);//s,e);
     //   print_in_order(g,MPI_COMM_WORLD);//matrix(g);
//	exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
   //     load_grid2d(g,s,e,s2,e2);
    //    MPI_Barrier(MPI_COMM_WORLD);
//        print_in_order(g, MPI_COMM_WORLD);
	
	//Save observables and free memory
	if(rank==0){
		char title[100];
	      	if(hex==1 && c==0) snprintf(title, 100, "stats_hexagonal_parallel_hot.txt");
        	if(triangle==1 && c==0) snprintf(title, 100, "stats_triangular_parallel_hot.txt");
        	if(sq==1 && c==0) snprintf(title, 100, "stats_square_parallel_hot.txt");
                if(hex==1 && c==1) snprintf(title, 100, "stats_hexagonal_parallel_cold.txt");
                if(triangle==1 && c==1) snprintf(title, 100, "stats_triangular_parallel_cold.txt");
                if(sq==1 && c==1) snprintf(title, 100, "stats_square_parallel_cold.txt");

		write_stats(title, magnetisations, energies, specs, sus, temps, r2);
		free(magnetisations);
		free(energies);
		free(temps);
		free(specs);
		free(sus);
	}
	}
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
	//printf("%d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot);

	decomp1d(m, ddims2[0], rank%ddims2[0], &s2, &e2); //Decompose cols
	decomp1d(n, ddims2[1], (rank-(rank%ddims2[0]))/ddims2[0], &s, &e); //Decompose rows

	printf("%d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, s, e, s2, e2);

	//Change prng for each process
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);
	//srand48(seed*rank+1);
	
	//printf("%d %d %d %d %d\n", rank, s, e, nbrtop, nbrbot);

	//Allocate memory to store the required stats
	double *magnetisations, *energies, *specs, *sus, *temps;
	if(rank==0){
		magnetisations = (double *)malloc(r2*sizeof(double));
		energies = (double *)malloc(r2*sizeof(double));
		temps = (double *)malloc(r2*sizeof(double));
		specs = malloc(r2*sizeof(double));
		sus = malloc(r2*sizeof(double));
	}

	//Perform the simulations and calculate an average for each observables at each temperature
	if(sq==1){
	for(int k=0; k<r1; k++){
		T=0.0001; //reset the temperature for each simulation
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int i=0; i<r3; i++){
				metropolis_sweep2d(g,s,e,s2,e2,J,t,gsl_mt);
				exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
				/*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);*/
				//if(rank%2==i%2) metropolis_sweep(g,s,e,J,T,k);
				//exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
			}
			//calculate the magnetisation per site for each process then get the average on the root process
			double mag, enrg, mag2, enrg2;
			double M = magnetisation2d(g,s,e,s2,e2);
			double E = energy2d(g,J,s,e,s2,e2);
			double M2 = var_mag2d(g,s,e,s2,e2,M);//magnetisation2(g,s,e);
			double E2 = var_enrg2d(g,s,e,s2,e2,J,E);//energy2(g,J,s,e); 
			MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				mag/=size;
				enrg/=size;
				mag2/=size;
				enrg2/=size;
				magnetisations[j]+=mag/r1;
				energies[j]+=enrg/r1;
				specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
				sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
				if(k==0) temps[j] = T;
			}
			T+=0.05;
			//Reset grid to original arrangement
			//reset_grid(G,g);
			init_grid(g,c,gsl_mt);
			exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
		}
	}
	}

	if(triangle==1){
        for(int k=0; k<r1; k++){
                T=0.0001; //reset the temperature for each simulation
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_triangular2d(g,s,e,s2,e2,J,t,gsl_mt);
                                exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                                /*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);*/
                                //if(rank%2==i%2) metropolis_sweep_triangular(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the magnetisation per site for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation2d(g,s,e,s2,e2);
                        double E = energy_triangular2d(g,J,s,e,s2,e2);
			double M2 = var_mag2d(g,s,e,s2,e2,M);//magnetisation2(g,s,e);
			double E2 = var_enrg_tri2d(g,s,e,s2,e2,J,E);//energy_triangular2(g,J,s,e);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
                                sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.07;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,gsl_mt);
			exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                }
        }
	}

	if(hex==1){
        for(int k=0; k<r1; k++){
                T=0.0001; //reset the temperature for each simulation
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_hexagonal2d(g,s,e,s2,e2,J,t,gsl_mt);
                                exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                                /*if(rank%2==0) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(rank%2==1) metropolis_sweep_triangular(g,s,e,J,T,k);
                                exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                                MPI_Barrier(MPI_COMM_WORLD);*/
                                //if(rank%2==i%2) metropolis_sweep_hexagonal(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the magnetisation per site for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation2d(g,s,e,s2,e2);
                        double E = energy_hexagonal2d(g,J,s,e,s2,e2);
			double M2 = var_mag2d(g,s,e,s2,e2,M);//magnetisation2(g,s,e);
			double E2 = var_enrg_hex2d(g,s,e,s2,e2,J,E);//energy_hexagonal2(g,J,s,e);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t/r1;//fabs(E2 - ((enrg)*(enrg)))*T*T/r1;
                                sus[j] += mag2*t/r1;//fabs(M2 - ((mag)*(mag)))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.05;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,gsl_mt);
			exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                }
        }
	}

        //Test MPI-IO
        //print_matrix(g);
//	init_grid2d2(g,s,e,s2,e2);
 //       save_grid2d(g,s,e,s2,e2);
//	init_grid(g,c,gsl_mt);
 //       MPI_Barrier(MPI_COMM_WORLD);


/*
        init_grid2d2(g,s,e,s2,e2);//0,n-1, 0, m-1);//s,e);
        print_in_order(g,MPI_COMM_WORLD);//matrix(g);
	exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
   //     load_grid2d(g,s,e,s2,e2);
    //    MPI_Barrier(MPI_COMM_WORLD);
        print_in_order(g, MPI_COMM_WORLD);

*/

	//Save observables and free memory
	if(rank==0){
		char title[100];
        	if(hex==1 && c==0) snprintf(title, 100, "stats_hexagonal_parallel_hot_2d.txt");
        	if(triangle==1 && c==0) snprintf(title, 100, "stats_triangular_parallel_hot_2d.txt");
        	if(sq==1 && c==0) snprintf(title, 100, "stats_square_parallel_hot_2d.txt");
                if(hex==1 && c==1) snprintf(title, 100, "stats_hexagonal_parallel_cold_2d.txt");
                if(triangle==1 && c==1) snprintf(title, 100, "stats_triangular_parallel_cold_2d.txt");
                if(sq==1 && c==1) snprintf(title, 100, "stats_square_parallel_cold_2d.txt");

		write_stats(title, magnetisations, energies, specs, sus, temps, r2);
		free(magnetisations);
		free(energies);
		free(temps);
		free(specs);
		free(sus);
	}	
	}
	MPI_Finalize();
return 0;
}

//Initialise grid for random start
void init_grid(int g[m][n], int c, gsl_rng *gsl_mt){
	int i,j;
	if(c==1){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				/*double U = drand48();
				if(U<0.5) g[i][j]=-1;
				else{ 
					g[i][j]=1;
				}*/
				g[i][j]=1;
			}
		}
	}
	if(c==0){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j]=-1;
				else{ g[i][j]=1;}
			}
		}
	}
}

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
	int i,j;
	double dE;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			//Calculate the energy difference and check if spin has to be flipped
			dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
			if(dE<=0) g[i][j]*=-1;
			else{
				if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
			}
		}
	}
}

void metropolis_sweep2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

void exchange(int g[m][n], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
	MPI_Request reqs[4];

/*	MPI_Irecv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[1]);
        MPI_Isend(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[3]);
*/
        MPI_Irecv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[1]);
        MPI_Isend(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[3]);


	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}

void exchange2d(int g[m][n], int s, int e, int s2, int e2, int nbrtop, int nbrbot, int nbrleft, int nbrright, MPI_Comm comm){
	MPI_Request reqs[8];

	//Send columns
	MPI_Datatype col;
	MPI_Type_vector(e-s, 1, n, MPI_INT, &col);//e-s
	MPI_Type_commit(&col);
	MPI_Irecv(&g[s][modulo(e2+1,n)], 1, col, nbrright, 0, comm, &reqs[0]);
	MPI_Irecv(&g[s][modulo(s2-1,n)], 1, col, nbrleft, 1, comm, &reqs[1]);
	MPI_Isend(&g[s][s2], 1, col, nbrleft, 0, comm, &reqs[2]);
	MPI_Isend(&g[s][e2], 1, col, nbrright, 1, comm, &reqs[3]);

	//Send rows
	MPI_Irecv(&g[modulo(e+1,m)][s2], e2-s2+1, MPI_INT, nbrbot, 2, comm, &reqs[4]);
	MPI_Irecv(&g[modulo(s-1,m)][s2], e2-s2+1, MPI_INT, nbrtop, 3, comm, &reqs[5]);
	MPI_Isend(&g[modulo(s,m)][s2], e2-s2+1, MPI_INT, nbrtop, 2, comm, &reqs[6]);
	MPI_Isend(&g[modulo(e,m)][s2], e2-s2+1, MPI_INT, nbrbot, 3, comm, &reqs[7]);

	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&col);
}

void exchange2(int g[m][n], int s, int e, int nbrtop, int nbrbot, int rank, int turn, MPI_Comm comm){
	if(turn%2==0){
		if(rank%2==0){
			//printf("a\n");
			MPI_Send(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm);
			//printf("b\n");
			MPI_Send(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm);
		}
		if(rank%2==1){
			//printf("c\n");
			MPI_Recv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrtop, 0, comm, MPI_STATUS_IGNORE);
			//printf("d\n");
			MPI_Recv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrbot, 0, comm, MPI_STATUS_IGNORE);
		}
	}
        if(turn%2==1){
                if(rank%2==1){
			//printf("e\n");
                        MPI_Send(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm);
			//printf("f\n");
			MPI_Send(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm);
                }
                if(rank%2==0){
			//printf("g\n");
			MPI_Recv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, MPI_STATUS_IGNORE);
                        //printf("h\n");
			MPI_Recv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, MPI_STATUS_IGNORE);
                }
        }
}

void metropolis_sweep_triangular(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
	double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]  + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

void metropolis_sweep_triangular2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]  + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}


void metropolis_sweep_hexagonal(int g[m][n], int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
	double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}

void metropolis_sweep_hexagonal2d(int g[m][n], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        if(dE<=0) g[i][j]*=-1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE*T)) g[i][j]*=-1;
                        }
                }
        }
}


/*void gather(int g[m][n], int s, int e, int rank, int size){
	int *recvptr = &(g[0][0]);
	int recvcounts[size];
	//recvcounts[0]=e-s+1;
	int displs[size];
	//displs[0]=s*n;
	//int *recvcounts, *displs;
	//recvcounts = (int *)malloc((e-s+1)*sizeof(int));// *sizeof(int);
	//displs = (int *)malloc(s*n*sizeof(int));
	//printf("%d %d %d %d\n", recvcounts[0], displs[0], e-s+1, s*n);
	int rc = e-s+1;
	if(rank==0){
		
	}
	MPI_Gatherv(&(g[0][0]), n*(e-s+1), MPI_INT, recvptr, &recvcounts[rank], &displs[rank], MPI_INT, 0, MPI_COMM_WORLD);
}*/

/*void IO(int g[m][n], int s, int e, int rank, int size){
	int dimx = m;
	int dimy = n;
	int subdimx = e-s+1;
	int subdimy = n;
	int sizes[2] = {dimx, dimy};
	int starts[2];
	

}*/

//Calculate the energy per site
double energy(int g[m][n], const double J, int s, int e){
	int i,j;
	double E=0;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
		}
	}
return E/(2.0*(e-s+1)*n);
}

double energy2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}


//Calculate the energy per site of the system
/*double energy2(int g[m][n], const double J, int s, int e){
        int i,j;
        double dE;
        double E=0;
        for(i=s; i<e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
                        E+=dE*dE;
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
}*/

double var_enrg(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


double energy_triangular(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_triangular2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}


/*double energy_triangular2(int g[m][n], const double J, int s, int e){
        int i,j;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
			E+=dE*dE;
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
}*/

double var_enrg_tri(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
			dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_tri2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


double energy_hexagonal(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_hexagonal2d(int g[m][n], const double J, int s, int e, int s2, int e2){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1));
}


/*double energy_hexagonal2(int g[m][n], const double J, int s, int e){
        int i,j;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
			E+=dE*dE;
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
}*/

double var_enrg_hex(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
			dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_hex2d(int g[m][n], int s, int e, int s2, int e2, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


//Calculate magnetisation per site
double magnetisation(int g[m][n], int s, int e){
	int i,j;
	double M=0;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			M+=g[i][j];
		}
	}
return fabs(M/((e-s+1)*n));
}

//Calculate magnetisation per site
double magnetisation2d(int g[m][n], int s, int e, int s2, int e2){
        int i,j;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        M+=g[i][j];
                }
        }
return fabs(M/((e-s+1)*(e2-s2+1)));
}


//Calculate magnetisation per site
/*double magnetisation2(int g[m][n], int s, int e){
        int i,j;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        M+=g[i][j]*g[i][j];
                }
        }
return fabs(M/((e-s+1)*n*(e-s+1)*n));
}*/

double var_mag(int g[m][n], int s, int e, double avg){
        int i,j;
        //double dM;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //dM = g[i][j]-avg;
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_mag2d(int g[m][n], int s, int e, int s2, int e2, double avg){
        int i,j;
        //double dM;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        //dM = g[i][j]-avg;
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/((e-s+1)*(e2-s2+1));
}


void print_matrix(int g[m][n]){
	int i,j;
	printf("\n");
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%d ", g[i][j]);
		}
	printf("\n");
	}
}


void print_in_order(int g[m][n], MPI_Comm comm)
{
  int myid, size;
  int i;

  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &size);
  MPI_Barrier(comm);
  printf("Attempting to print in order\n");
  sleep(1);
  MPI_Barrier(comm);

  for(i=0; i<size; i++){
    if( i == myid ){
      printf("proc %d\n",myid);
      print_matrix(g);
    }
    fflush(stdout);
    usleep(500);	
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}

void write_stats(char *title, double *mag, double *energy, double *specs, double *sus, double *T, int r1){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r1; i++){
		fprintf(fp,"%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], specs[i]);
	}
	fclose(fp);
}

//Reset grid after every temp increase
void reset_grid(int G[m][n], int g[m][n]){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			g[i][j]=G[i][j];
		}
	}
}

void save_grid(int g[m][n], int s, int e){
	MPI_Datatype subarray;
	int ndims=2;
	int sizes[2] = {m,n};
	int subsizes[2] = {(e-s+1),n};
	int starts[2] = {s,0};
	MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
	MPI_Type_commit(&subarray);

	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, /*s*n*/0, MPI_INT, subarray, "native", MPI_INFO_NULL);
	MPI_File_write_all(fh, &g[s][0]/*[s][0]*/, (e-s+1)*n, /*subarray*/ MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	MPI_Type_free(&subarray);
}

void save_grid2d(int g[m][n], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),(e2-s2+1)};
        int starts[2] = {s,s2};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, /*s*n*/0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &g[s][s2], /*1*/(e-s+1)*n, /*subarray*/ MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void load_grid(int g[m][n], int s, int e){
	MPI_Datatype subarray;
	int ndims=2;
	int sizes[2] = {m,n};
	int subsizes[2] = {(e-s+1),n};
	int starts[2] = {s,0};
	MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
	MPI_Type_commit(&subarray);

	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, 0, subarray, subarray, "native", MPI_INFO_NULL);
	MPI_File_read_all(fh, &g[0][0], 1, subarray, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	MPI_Type_free(&subarray);



//	Alternative code to load full grid instead of just relevant part

//	MPI_File fh1;
//	MPI_File_open(MPI_COMM_WORLD, "configuration", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh1);
//	MPI_File_set_view(fh1, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//	MPI_File_read_all(fh1, &g[0][0], m*n, MPI_INT, MPI_STATUS_IGNORE);
//	MPI_File_close(&fh1);

}

void load_grid2d(int g[m][n], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
	//printf("%d %d %d %d\n", s, e, s2, e2);
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),(e2-s2+1)};
        int starts[2] = {s,s2};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, /*((s*n)+(e2-s2+1))*sizeof(int)*/0, /*MPI_INT*/subarray, /*MPI_INT*/subarray, "native", MPI_INFO_NULL);
	MPI_File_read_all(fh, &g[s][s2], (e-s+1)*n/*(e2-s2+1)*sizeof(int)*/, MPI_INT, MPI_STATUS_IGNORE);
        //MPI_File_read_all(fh, &g[s][s2], (e-s+1)*(e2-s2+1)*sizeof(int), /*subarray*/MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);

//      Alternative code to load full grid instead of just relevant part

//      MPI_File fh1;
//      MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh1);
//      MPI_File_set_view(fh1, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//      MPI_File_read_all(fh1, &g[0][0], m*n, MPI_INT, MPI_STATUS_IGNORE);
//      MPI_File_close(&fh1);
}


void init_grid2(int g[m][n], int s, int e){
	for(int i=s; i<=e; i++){
		for(int j=0; j<n; j++){
			g[i][j]=i*m+j;
		}
	}
}

void init_grid2d2(int g[m][n], int s, int e, int s2, int e2){
        for(int i=s; i<=e; i++){
                for(int j=s2; j<=e2; j++){
                        g[i][j]=i*m+j;
                }
        }
}
