/*#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include"decomp1d.h"

#define m 50
#define n 50
*/
//#include"grid_size.h"
#include"decomp1d.h"
#include"exchange.h"
#include"metropolis.h"
#include"observables.h"

void init_grid(int g[m][n], int, int, gsl_rng *);
/*
void metropolis_sweep(int g[m][n], int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], int, int, int, const double , double, gsl_rng *);
*/
/*
double energy(int g[m][n], const double, int, int);
double energy2(int g[m][n], const double, int, int);
double energy_triangular(int g[m][n], const double, int, int);
double energy_triangular2(int g[m][n], const double, int, int);
double energy_hexagonal(int g[m][n], const double, int, int);
double energy_hexagonal2(int g[m][n], const double, int, int);
double magnetisation(int g[m][n], int, int);
double magnetisation2(int g[m][n], int, int);
*/
void print_matrix(int g[m][n]);
/*
void exchange(int g[m][n], int, int, int, int, MPI_Comm);
void exchange2(int g[m][n], int, int, int, int, int, int, MPI_Comm);
*/
//void gather(int g[m][n], int, int, int, int);
void print_in_order(int g[m][n], MPI_Comm comm);
//int modulo(int, int);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
//int delta(int, int);
//int H(int, int);
/*
double var_mag(int g[m][n], int, int, double);
double var_enrg(int g[m][n], int, int, double, double);
double var_enrg_tri(int g[m][n], int, int, double, double);
double var_enrg_hex(int g[m][n], int, int, double, double);
*/
int main(int argc, char **argv){
	int rank, size, s, e, nbrtop, nbrbot;
	int spin = 2;
	double T = 0.0001; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	int c=0;

	//const double k = pow(1.38064852, -23); //Boltzmann constant
	double r1 = 50.0; //Number of simulations
	int r2 = 100; //Number of temperatures
	int r3 = 1000; //Number of sweeps

	//Seed prng and inititate grid
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);
	//int seed = 1997;
	//srand48(seed);

	int g[m][n]; //define grid
	//int G[m][n];
	init_grid(g, c, spin, gsl_mt);
	//reset_grid(G,g);


        //Parse variables
        int sq=0;
        int triangle=0;
        int hex=0;
        int opt;
        while((opt = getopt(argc, argv, "123hc")) != -1){
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

	//Create new communicator for finding neighbours
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart);
	MPI_Cart_shift(cart, 0, -1, &source1, &nbrtop);
	MPI_Cart_shift(cart, 0, 1, &source2, &nbrbot);

	//Find size of each subarray
	decomp1d(n, size, rank, &s, &e);
	
	//Change prng for each process
	srand48(seed*rank+1);
	
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
		double t;
		T=0.0001; //reset the temperature for each simulation
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int i=0; i<2*r3; i++){
				metropolis_sweep(g,spin,s,e,J,T,gsl_mt);
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
			double M2 = var_mag(g,s,e,M);
			double E2 = var_enrg(g,s,e,J,E);
			MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				mag/=size;
				enrg/=size;
				mag2/=size;
				enrg2/=size;
				//printf("%f %f %f\n", mag, enrg, T);
				magnetisations[j]+=mag/r1;
				energies[j]+=enrg/r1;
				specs[j] += enrg2*t*t/r1;//fabs(enrg2 - (enrg/r1)*(enrg/r1))*(T*T)/r1;
				sus[j] += mag2*t/r1;//fabs(mag2 - (mag/r1)*(mag/r1))*T/r1;
				if(k==0) temps[j] = T;
			}
			T+=0.1;
			//Reset grid to original arrangement
			//reset_grid(G,g);
			init_grid(g,c,spin,gsl_mt);
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
                        for(int i=0; i<2*r3; i++){
                                metropolis_sweep_triangular(g,spin,s,e,J,T,gsl_mt);
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
			double M2 = var_mag(g,s,e,M);
			double E2 = var_enrg_tri(g,s,e,J,E);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                //printf("%f %f %f\n", mag, enrg, T);
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t;//fabs(enrg2 - (enrg/r1)*(enrg/r1))*(T*T)/r1;
                                sus[j] += mag2*t/r1;//fabs(mag2 - (mag/r1)*(mag/r1))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.1;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,spin,gsl_mt);
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
                        for(int i=0; i<2*r3; i++){
                                metropolis_sweep_triangular(g,spin,s,e,J,T,gsl_mt);
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
			double M2 = var_mag(g,s,e,M);
			double E2 = var_enrg_hex(g,s,e,J,E);
                        MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &mag2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &enrg2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                mag/=size;
                                enrg/=size;
				mag2/=size;
				enrg2/=size;
                                //printf("%f %f %f\n", mag, enrg, T);
                                magnetisations[j]+=mag/r1;
                                energies[j]+=enrg/r1;
                                specs[j] += enrg2*t*t/r1;//fabs(enrg2 - (enrg/r1)*(enrg/r1))*(T*T)/r1;
                                sus[j] += mag2*t/r1;//fabs(mag2 - (mag/r1)*(mag/r1))*T/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.1;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,spin,gsl_mt);
			exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
                }
        }
	}
	//Save observables and free memory
	if(rank==0){
		char title[100];
        	if(hex==1 && c==0) snprintf(title, 100, "potts2d_stats_hexagonal_parallel_hot.txt");
        	if(triangle==1 && c==0) snprintf(title, 100, "potts2d_stats_triangular_parallel_hot.txt");
        	if(sq==1 && c==0) snprintf(title, 100, "potts2d_stats_square_parallel_hot.txt");
                if(hex==1 && c==1) snprintf(title, 100, "potts2d_stats_hexagonal_parallel_cold.txt");
                if(triangle==1 && c==1) snprintf(title, 100, "potts2d_stats_triangular_parallel_cold.txt");
                if(sq==1 && c==1) snprintf(title, 100, "potts2d_stats_square_parallel_cold.txt");

		write_stats(title, magnetisations, energies, specs, sus, temps, r2);
		free(magnetisations);
		free(energies);
		free(temps);
		free(specs);
		free(sus);
	}
	MPI_Finalize();
return 0;
}

//Initialise grid for random start
void init_grid(int g[m][n], int c, int s, gsl_rng *gsl_mt){
	int i,j;
	if(c==1){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				/*double U = drand48();
				if(U<0.5) g[i][j]=-1;
				else{ 
					g[i][j]=1;
				}*/
				g[i][j]=s;
			}
		}
	}
	if(c==0){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				g[i][j]=gsl_rng_uniform_int(gsl_mt,s)+1;
			}
		}
	}
}
/*
//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
	int i,j,r;
	double dE;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			//Calculate the energy difference and check if spin has to be flipped
			r=gsl_rng_uniform_int(gsl_mt,spin)+1;
			dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[i][modulo(j+1,n)]));
			if(dE<=0) g[i][j]=r;
			else{
				if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
			}
		}
	}
}
*/
/*
void exchange(int g[m][n], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
	MPI_Request reqs[4];

	MPI_Irecv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[1]);
        MPI_Isend(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[3]);

	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
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
*/
/*
void metropolis_sweep_triangular(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
	double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,spin)+1;
			dE = 2.0*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r, g[i][modulo(j+1,n)]) + H(r,g[i][modulo(j-1,n)]) + H(r,g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(r,g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}

void metropolis_sweep_hexagonal(int g[m][n], int spin, int s, int e, const double J, double T, gsl_rng *gsl_mt){
        int i,j,r;
	double dE;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        r=gsl_rng_uniform_int(gsl_mt,s)+1;
			dE = 2*J*(H(r,g[modulo(i+1,m)][j]) + H(r,g[modulo(i-1,m)][j]) + H(r,g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        if(dE<=0) g[i][j]=r;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<exp(-dE/T)) g[i][j]=r;
                        }
                }
        }
}
*/

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
/*
//Calculate the energy per site
double energy(int g[m][n], const double J, int s, int e){
	int i,j;
	double E=0;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			E += -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
			//E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
		}
	}
return E/(2.0*(e-s+1)*n);
}

double energy2(int g[m][n], const double J, int s, int e){
        int i,j;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
			E+=dE*dE;
                        //E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
}


double energy_triangular(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
			E+=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
                        //E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_triangular2(int g[m][n], const double J, int s, int e){
        int i,j;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]));
			E+=dE*dE;
                        //E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
}
*/
/*
double energy_hexagonal(int g[m][n], const double J, int s, int e){
        int i,j;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
			E+=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        //E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n);
}

double energy_hexagonal2(int g[m][n], const double J, int s, int e){
        int i,j;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        dE=-J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
			E+=dE*dE;
                        //E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*(e-s+1)*n*(e-s+1)*n);
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

double magnetisation2(int g[m][n], int s, int e){
        int i,j;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        M+=g[i][j]*g[i][j];
                }
        }
return fabs(M/((e-s+1)*n*(e-s+1)*n));
}
*/

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
/*
int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}
*/
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
double var_mag(int g[m][n], int s, int e, double avg){
        int i,j;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<n; j++){
                        var += (g[i][j] - avg)*(g[i][j] - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]));
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_tri(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+1,n)]) + H(g[i][j],g[i][modulo(j-1,n)]) + H(g[i][j],g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)]) + H(g[i][j],g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)])); 
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}

double var_enrg_hex(int g[m][n], int s, int e, double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=0; j<m; j++){
                        dE = -J*(H(g[i][j],g[modulo(i+1,m)][j]) + H(g[i][j],g[modulo(i-1,m)][j]) + H(g[i][j],g[i][modulo(j+((int)pow(-1,i+j)),n)]));
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/((e-s+1)*n);
}
*/
