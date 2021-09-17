#include"decomp1d.h"
#include"modulo.h"
#include"mpi_io.h"
#include"observables.h"
#include"metropolis.h"
#include"exchange.h"

void init_grid(int g[m][n], int, int, gsl_rng *);
void print_matrix(int g[m][n]);
void print_in_order(int g[m][n], MPI_Comm comm);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
void init_grid2(int g[m][n], int, int);
void init_grid2d2(int g[m][n], int, int, int, int);

int main(int argc, char **argv){
	int rank, size, s, e, s2, e2, nbrtop, nbrbot, nbrleft, nbrright, nbrd1, nbrd2, nbrd3, nbrd4;
	double T = 0.0001; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	//const double k = pow(1.38064852, -23); //Boltzmann constant
	int c=0;

	int spin=5;

	double r1 = 50.0; //Number of simulations
	int r2 = 100; //Number of temperatures
	int r3 = 1000; //Number of sweeps

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
			case 'c':
				c=1;
				break;
			case 'h':
				printf("Command lines:\n -1 for square grid\n -2 for triangular grid \n -3 for hexagonal grid\n -c for cold start (hot start default)\n -d for two dimensional exchange");
				return 0;
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

	//Initiate grid
	init_grid(g,c,spin,gsl_mt);

	//Variables for MPI Cart
	int ndims;
	//int dims[1];
	//int periods[1];
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

	decomp1d(m, ddims2[0], rank%ddims2[0], &s2, &e2); //Decompose cols
	decomp1d(n, ddims2[1], (rank-(rank%ddims2[0]))/ddims2[0], &s, &e); //Decompose rows

	//printf("%d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, s, e, s2, e2);

	//Change prng for each process
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);
	//srand48(seed*rank+1);
	

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
				metropolis_sweep_2d(g,spin,s,e,s2,e2,J,t,gsl_mt);
				exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
				//if(rank%2==i%2) metropolis_sweep(g,s,e,J,T,k);
				//exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
			}
			//calculate the observables for each process then get the average on the root process
			double mag, enrg, mag2, enrg2;
			double M = magnetisation_2d(g,s,e,s2,e2);
			double E = energy_2d(g,J,s,e,s2,e2);
			double M2 = var_mag_2d(g,s,e,s2,e2,M);
			double E2 = var_enrg_2d(g,s,e,s2,e2,J,E); 
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
				specs[j] += enrg2*t*t/r1;
				sus[j] += mag2*t/r1;
				if(k==0) temps[j] = T;
			}
			T+=0.05;
			//Reset grid to original arrangement
			init_grid(g,c,spin,gsl_mt);
			exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
		}
	}
	}

	if(triangle==1){

	//Find diagonal neighbours
	nbrd1 = modulo(nbrleft-ddims2[0],size);
	nbrd2 = modulo(nbrright-ddims2[0],size);
	nbrd3 = modulo(nbrleft+ddims2[0],size);
	nbrd4 = modulo(nbrright+ddims2[0],size);
        printf("%d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, nbrd1, nbrd2, nbrd3, nbrd4);

        for(int k=0; k<r1; k++){
                T=0.0001; //reset the temperature for each simulation
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_triangular_2d(g,spin,s,e,s2,e2,J,t,gsl_mt);
                                exchange2d_tri(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,nbrd1,nbrd2,nbrd3,nbrd4,MPI_COMM_WORLD);


                                //if(rank%2==i%2) metropolis_sweep_triangular(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the observables for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation_2d(g,s,e,s2,e2);
                        double E = energy_triangular_2d(g,J,s,e,s2,e2);
			double M2 = var_mag_2d(g,s,e,s2,e2,M);
			double E2 = var_enrg_tri_2d(g,s,e,s2,e2,J,E);
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
                                specs[j] += enrg2*t*t/r1;
                                sus[j] += mag2*t/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.07;
                        //Reset grid to original arrangement
                        init_grid(g,c,spin,gsl_mt);
			exchange2d_tri(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,nbrd1,nbrd2,nbrd3,nbrd4,MPI_COMM_WORLD);
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
                                metropolis_sweep_hexagonal_2d(g,spin,s,e,s2,e2,J,t,gsl_mt);
                                exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
                                
				//if(rank%2==i%2) metropolis_sweep_hexagonal(g,s,e,J,T,k);
                                //exchange2(g,s,e,nbrtop,nbrbot,rank,i,MPI_COMM_WORLD);
                        }
                        //calculate the observables for each process then get the average on the root process
                        double mag, enrg, mag2, enrg2;
                        double M = magnetisation_2d(g,s,e,s2,e2);
                        double E = energy_hexagonal_2d(g,J,s,e,s2,e2);
			double M2 = var_mag_2d(g,s,e,s2,e2,M);
			double E2 = var_enrg_hex_2d(g,s,e,s2,e2,J,E);
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
                                specs[j] += enrg2*t*t/r1;
                                sus[j] += mag2*t/r1;
                                if(k==0) temps[j] = T;
                        }
                        T+=0.05;
                        //Reset grid to original arrangement
                        //reset_grid(G,g);
                        init_grid(g,c,spin,gsl_mt);
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



	//Test exchanges
    //    init_grid2d2(g,s,e,s2,e2);//0,n-1, 0, m-1);//s,e);
  //      print_in_order(g,MPI_COMM_WORLD);//matrix(g);
//	exchange2d(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,MPI_COMM_WORLD);
   //     load_grid2d(g,s,e,s2,e2);
    //    MPI_Barrier(MPI_COMM_WORLD);
//        print_in_order(g, MPI_COMM_WORLD);

	//Test exchanges for triangular lattice
//        nbrd1 = modulo(nbrleft-ddims2[0],size);
//        nbrd2 = modulo(nbrright-ddims2[0],size);
//        nbrd3 = modulo(nbrleft+ddims2[0],size);
//        nbrd4 = modulo(nbrright+ddims2[0],size);
//	init_grid2d2(g,s,e,s2,e2);
//	print_in_order(g,MPI_COMM_WORLD);
//	exchange2d_tri(g,s,e,s2,e2,nbrtop,nbrbot,nbrleft,nbrright,nbrd1,nbrd2,nbrd3,nbrd4,MPI_COMM_WORLD);
//	print_in_order(g, MPI_COMM_WORLD);	

	//Save observables and free memory
	if(rank==0){
		char title[100];
        	if(hex==1 && c==0) snprintf(title, 100, "potts2d_stats_hexagonal_parallel_hot_2d.txt");
        	if(triangle==1 && c==0) snprintf(title, 100, "potts2d_stats_triangular_parallel_hot_2d.txt");
        	if(sq==1 && c==0) snprintf(title, 100, "potts2d_stats_square_parallel_hot_2d.txt");
                if(hex==1 && c==1) snprintf(title, 100, "potts2d_stats_hexagonal_parallel_cold_2d.txt");
                if(triangle==1 && c==1) snprintf(title, 100, "potts2d_stats_triangular_parallel_cold_2d.txt");
                if(sq==1 && c==1) snprintf(title, 100, "potts2d_stats_square_parallel_cold_2d.txt");

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

//write stats to file
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

//grid initiations for testing code
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
