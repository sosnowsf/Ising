//#include"grid_size.h"
#include"decomp1d.h"
#include"exchange.h"
#include"metropolis.h"
#include"observables.h"
#include"mpi_io.h"

void init_grid(int g[m][n], int, int, gsl_rng *);
void print_matrix(int g[m][n]);
void print_in_order(int g[m][n], MPI_Comm comm);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);

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
