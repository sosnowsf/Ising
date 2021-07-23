#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include"decomp1d.h"

#define m 10
#define n 10

void init_grid(int g[m][n]);
void metropolis_sweep(int g[m][n], int, int, const double , double, const double);
double energy(int g[m][n], const double, int, int);
double magnetisation(int g[m][n], int, int);
void print_matrix(int g[m][n]);
void exchange(int g[m][n], int, int, int, int, MPI_Comm comm);
//void gather(int g[m][n], int, int, int, int);
void print_in_order(int g[m][n], MPI_Comm comm);
int modulo(int, int);
void write_stats(double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);

int main(int argc, char **argv){
	int rank, size, s, e, nbrtop, nbrbot;
	double T = 0; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	const double k = pow(1.38064852, -23); //Boltzmann constant

	double r1 = 100;
	int r2 = 100;
	int r3 = 100;

	//Seed prng and inititate grid
	int seed = 1997;
	srand48(seed);

	int g[m][n]; //define grid
	int G[m][n];
	init_grid(G);
	reset_grid(G,g);

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
	
	printf("%d %d %d %d %d\n", rank, s, e, nbrtop, nbrbot);
	/*for(int i=s; i<e+1; i++){
		for(int j=0; j<n; j++){
			g[i][j]=i*n+j;
		}
	}*/
	double *magnetisations, *energies, *specs, *sus, *temps;
	if(rank==0){
		magnetisations = (double *)malloc(r1*sizeof(double));
		energies = (double *)malloc(r1*sizeof(double));
		temps = (double *)malloc(r1*sizeof(double));
	}
	for(int k=0; k<r1; k++){
		T=0;
		for(int j=0; j<r2; j++){
			for(int i=0; i<r3; i++){
				metropolis_sweep(g,s,e,J,T,k);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
			}

			double mag, enrg;
			double M = magnetisation(g,s,e);
			double E = energy(g,J,s,e);
			MPI_Reduce(&M, &mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E, &enrg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				mag/=size;
				enrg/=size;
				//printf("%f %f %f\n", mag, enrg, T);
				magnetisations[j]+=mag/r1;
				energies[j]+=enrg/r1;
				if(k==0) temps[j] = T;
			}
			T+=0.05;
			//Reset grid to original arrangement
			reset_grid(G,g);
		}
	}
	if(rank==0){
		write_stats(magnetisations, energies, temps, r1);
	}
	//printf("%f %f\n", magnetisation(g,s,e), energy(g,J,s,e));
	MPI_Finalize();
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
}

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], int s, int e, const double J, double T, const double k){
	int i,j;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			//Calculate the energy difference and check if spin has to be flipped
			double dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
			if(dE<0) g[i][j]*=-1;
			else{
				double r = drand48();
				if(r<exp(-(1.0/T)*dE)){
					g[i][j]*=-1;
				}
			}
		}
	}
}

void exchange(int g[m][n], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
	MPI_Request reqs[4];

	MPI_Irecv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[1]);
        MPI_Isend(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[3]);

	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
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

//Calculate the energy of the system
double energy(int g[m][n], const double J, int s, int e){
	int i,j;
	double E=0;
	for(i=s; i<=e; i++){
		for(j=0; j<n; j++){
			E+=-J*g[i][j]*(g[(i+1)%m][j]+g[(i-1)&m][j]+g[i][(j+1)%n]+g[i][(j-1)%n]);
		}
	}
return E/(2.0*m*n);
}

//Calculate magnetisation persite
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

void write_stats(double *mag, double *energy, double *T, int r1){
	FILE *fp = fopen("stats.txt", "w");
	for(int i=0; i<r1; i++){
		fprintf(fp,"%f %f %f\n", T[i], mag[i], energy[i]);
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
