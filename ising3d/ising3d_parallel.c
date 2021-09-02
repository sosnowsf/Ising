#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include"decomp1d.h"

#define x 24
#define y 24
#define z 24

void init_grid(int g[x][y][z], int, gsl_rng *);
void metropolis_sweep(int g[x][y][z], int, int, const double , double, gsl_rng *);
void metropolis_sweep2d(int g[x][y][z], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep3d(int g[x][y][z], int, int, int, int, int, int, const double , double, gsl_rng *);
double energy(int g[z][y][z], const double, int, int);
double energy2d(int g[z][y][z], const double, int, int, int, int);
double energy3d(int g[z][y][z], const double, int, int, int, int, int, int);
double energy2(int g[z][y][z], const double, int, int);
double var_enrg(int g[x][y][z], int, int, double, double);
double var_enrg2d(int g[x][y][z], int, int, int, int, double, double);
double var_enrg3d(int g[x][y][z], int, int, int, int, int, int, double, double);
double magnetisation(int g[x][y][z], int, int);
double magnetisation2d(int g[x][y][z], int, int, int, int);
double magnetisation3d(int g[x][y][z], int, int, int, int, int, int);
double magnetisation2(int g[x][y][z], int, int);
double var_mag(int g[x][y][z], int, int, double);
double var_mag2d(int g[x][y][z], int, int, int, int, double);
double var_mag3d(int g[x][y][z], int, int, int, int, int, int, double);
void exchange(int g[x][y][z], int, int, int, int, MPI_Comm);
void exchange2d(int g[x][y][z], int, int, int, int, int, int, int, int, MPI_Comm);
void exchange3d(int g[x][y][z], int, int, int, int, int, int, int, int, int, int, int, int, MPI_Comm);
void write_stats(char *, double *, double *, double *, double *, double *, int);
int modulo(int, int);
void reset_grid(int G[x][y][z], int g[x][y][z]);

int main(int argc, char **argv){
	int rank, size, s, e, s2, e2, s3, e3, nbrtop, nbrbot, nbrleft, nbrright, nbrfront, nbrback; //variables needed for parallelisation
	double T; //Temperature
	const double J = 1.0; //Coupling constant

        int c=0;
	int d=1;
	char *dval;
        int opt;
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
	if(d!=1 && d!=2 && d!=3){
		perror("Invalid value for decomposition (can only be 1,2 or 3)");
		exit(EXIT_FAILURE);
	}

	int g[x][y][z];

	double r1 = 10.0; //Number of simulations
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

	//printf("%d %d %d\n", rank, s, e);

	for(int i=0; i<r1; i++){
		T=0.0001;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int k=0; k<r3; k++){
				metropolis_sweep(g,s,e,J,T,gsl_mt);
				exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
			}
			double M,E,m2,e2;
			double M1 = magnetisation(g,s,e);
			double E1 = energy(g,J,s,e);
			double M2 = var_mag(g,s,e,M1);
			double E2 = var_enrg(g,s,e,J,E1);
			//if(E1>0) printf("%d %d %d %f %f\n", rank, s, e, T, E1);
			MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M2, &m2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&E2, &e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank==0){
				M/=size;
				E/=size;
				//printf("%f %f\n", M, E);
				m2/=size;
				e2/=size;
				mag[j]+=M/r1;
				eng[j]+=E/r1;
				specs[j] += e2*t*t/r1;//fabs(e2 - (E*E))*T*T/r1;
				sus[j] += m2*t/r1;//fabs(m2 - (M*M))*T/r1;
				if(i==0) temps[j] = T;
			}
		T+=0.1;
		}
	init_grid(g,c,gsl_mt);
	//exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
	}

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
	//printf("%d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot);

	decomp1d(y, ddims2[0], rank%ddims2[0], &s2, &e2); //Decompose cols
	decomp1d(x, ddims2[1], (rank-(rank%ddims2[0]))/ddims2[0], &s, &e); //Decompose rows

	printf("%d %d %d %d %d %d %d %d %d\n", rank, nbrleft, nbrright, nbrtop, nbrbot, s, e, s2, e2);

	//Change prng for each process
	unsigned long seed2 = (1999*rank);
	gsl_rng_set(gsl_mt, seed2);
	//srand48(seed*rank+1);
	
	//printf("%d %d %d %d %d\n", rank, s, e, nbrtop, nbrbot);

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
			//printf("%d %d %d %d %d\n", rank, s, e, s2, e2);
                        double M1 = magnetisation2d(g,s,e,s2,e2);
                        double E1 = energy2d(g,J,s,e,s2,e2);
                        double M2 = var_mag2d(g,s,e,s2,e2,M1);
                        double E2 = var_enrg2d(g,s,e,s2,e2,J,E1);
			//printf("%d %f %f %f %f\n", rank, M1, E1, M2, E2);
                        MPI_Reduce(&M1, &M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E1, &E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&M2, &m2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&E2, &ee2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if(rank==0){
                                M/=size;
                                E/=size;
                                m2/=size;
                                ee2/=size;
				//printf("%f %f %f %f\n", M, E, m2, e2);
                                mag[j]+=M/r1;
                                eng[j]+=E/r1;
                                specs[j] += ee2*t*t/r1;//fabs(e2 - (E*E))*T*T/r1;
                                sus[j] += m2*t/r1;//fabs(m2 - (M*M))*T/r1;
                                if(i==0) temps[j] = T;
                        }
                T+=0.1;
                }
        init_grid(g,c,gsl_mt);
        //exchange(g,s,e,nbrtop,nbrbot,MPI_COMM_WORLD);
        }

        if(rank==0){
                char title[100];
                if(c==0) snprintf(title, 100, "stats_cube_parallel_hot_2d.txt");
                if(c==1) snprintf(title, 100, "stats_cube_parallel_cold_2d.txt");
                write_stats(title,mag,eng,specs,sus,temps,r2);
        }


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
					g[i][j][k]= 1;//2*gsl_rng_uniform_int(gsl_mt, 2)-1;
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
					if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j][k]=-1;
					else{ g[i][j][k]=1; }
				}
			}
		}
	}
}

void metropolis_sweep(int g[x][y][z], int s, int e, const double J, double T, gsl_rng *gsl_mt){
	int i,j,k;
	double dE;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
				if(dE<0) g[i][j][k]*=-1;
				else{
					if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
						g[i][j][k]*=-1;
					}
				}	
			}
		}
	}
}

void metropolis_sweep2d(int g[x][y][z], int s, int e, int s2, int e2, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                if(dE<0) g[i][j][k]*=-1;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
                                                g[i][j][k]*=-1;
                                        }
                                }
                        }
                }
        }
}

void metropolis_sweep3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, const double J, double T, gsl_rng *gsl_mt){
        int i,j,k;
        double dE;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                dE = 2.0*J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                if(dE<0) g[i][j][k]*=-1;
                                else{
                                        if(gsl_rng_uniform(gsl_mt)<(exp(-dE/T))){
                                                g[i][j][k]*=-1;
                                        }
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

double magnetisation2d(int g[x][y][z], int s, int e, int s2, int e2){
        int i,j,k;
        double M=0;
	//printf("%d %d %d %d\n", s, e, s2, e2);
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k];
                        }
                }
        }
	//printf("%f %d %d %d %d\n", M, s, e, s2, e2);
return fabs(M/((e-s+1)*(e2-s2+1)*z));
}

double magnetisation3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3){
        int i,j,k;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                M+=g[i][j][k];
                        }
                }
        }
return fabs(M/((e-s+1)*(e2-s2+1)*(e3-s2+1)));
}


double magnetisation2(int g[x][y][z], int s, int e){
        int i,j,k;
        double M=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                M+=g[i][j][k]*g[i][j][k];
                        }
                }
        }
return fabs(M/(((e-s+1)*y*z)*((e-s+1)*y*z)));
}

double var_mag(int g[x][y][z], int s, int e, double avg){
	int i,j,k;
	double var=0;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
			}
		}
	}
return var/((e-s+1)*y*z);
}

double var_mag2d(int g[x][y][z], int s, int e, int s2, int e2, double avg){
        int i,j,k;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*z);
}

double var_mag3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, double avg){
        int i,j,k;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                var+=(g[i][j][k]-avg)*(g[i][j][k]-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*(e3-s3+1));
}

double energy(int g[x][y][z], const double J, int s, int e){
	int i,j,k;
	double E=0;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
			}
		}
	}
return E/(2.0*(e-s+1)*y*z);
}

double energy2d(int g[x][y][z], const double J, int s, int e, int s2, int e2){
        int i,j,k;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                        }
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1)*z);
}

double energy3d(int g[x][y][z], const double J, int s, int e, int s2, int e2, int s3, int e3){
        int i,j,k;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                E+=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                        }
                }
        }
return E/(2.0*(e-s+1)*(e2-s2+1)*(e3-s3+1));
}

double var_enrg(int g[x][y][z], int s, int e, double J, double avg){
	int i,j,k;
	double dE;
	double var=0;
	for(i=s; i<=e; i++){
		for(j=0; j<y; j++){
			for(k=0; k<z; k++){
				dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
				dE/=2.0;
				var += (dE-avg)*(dE-avg);
			}
		}
	}
return var/((e-s+1)*y*z);
}

double var_enrg2d(int g[x][y][z], int s, int e, int s2, int e2, double J, double avg){
        int i,j,k;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=0; k<z; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
                                var += (dE-avg)*(dE-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*z);
}

double var_enrg3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, double J, double avg){
        int i,j,k;
        double dE;
        double var=0;
        for(i=s; i<=e; i++){
                for(j=s2; j<=e2; j++){
                        for(k=s3; k<=e3; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
                                dE/=2.0;
                                var += (dE-avg)*(dE-avg);
                        }
                }
        }
return var/((e-s+1)*(e2-s2+1)*(e3-s3+1));
}


double energy2(int g[x][y][z], const double J, int s, int e){
        int i,j,k;
	double dE;
        double E=0;
        for(i=s; i<=e; i++){
                for(j=0; j<y; j++){
                        for(k=0; k<z; k++){
                                dE=-J*g[i][j][k]*(g[modulo(i+1,x)][j][k] + g[modulo(i-1,x)][j][k] + g[i][modulo(j+1,y)][k] +g[i][modulo(j-1,y)][k] + g[i][j][modulo(k+1,z)] + g[i][j][modulo(k-1,z)]);
				E+=dE*dE;
                        }
                }
        }
return E/(2.0*(e-s+1)*y*z*(e-s+1)*y*z);
}


void exchange(int g[x][y][z], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
	MPI_Request reqs[4];
	MPI_Datatype slice;
	MPI_Type_vector(y*z, 1, 1, MPI_INT, &slice);
	MPI_Type_commit(&slice);

	MPI_Irecv(&g[modulo(e+1,x)][0][0], 1, slice, nbrtop, 0, comm, &reqs[0]);
	MPI_Irecv(&g[modulo(s-1,x)][0][0], 1, slice, nbrbot, 1, comm, &reqs[1]);
	MPI_Isend(&g[modulo(s,x)][0][0], 1, slice, nbrbot, 0, comm, &reqs[2]);
	MPI_Isend(&g[modulo(e,x)][0][0], 1, slice, nbrtop, 1, comm, &reqs[3]);


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

void exchange2d(int g[x][y][z], int s, int e, int s2, int e2, int nbrtop, int nbrbot, int nbrleft, int nbrright, MPI_Comm comm){
	MPI_Request reqs[8];
	MPI_Datatype slice1;
	//MPI_Type_vector((e-s)*z, 1, z, MPI_INT, &slice1);
	MPI_Type_vector((e-s-1/*+1*/), z, y*z, MPI_INT, &slice1);
	MPI_Type_commit(&slice1);

	//Send to left and right nbrs
        MPI_Irecv(&g[s+1][modulo(e2+1,x)][0], 1, slice1, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s+1][modulo(s2-1,x)][0], 1, slice1, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s+1][modulo(s2,x)][0], 1, slice1, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s+1][modulo(e2,x)][0], 1, slice1, nbrright, 1, comm, &reqs[3]);

	//Send to top and bot nbrs
        MPI_Datatype slice2;
        MPI_Type_vector((e2-s2+1)*z, 1, 1, MPI_INT, &slice2);
        MPI_Type_commit(&slice2);

        MPI_Irecv(&g[modulo(e+1,x)][s2][0], 1, slice2, nbrtop, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,x)][s2][0], 1, slice2, nbrbot, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,x)][s2][0], 1, slice2, nbrbot, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,x)][s2][0], 1, slice2, nbrtop, 3, comm, &reqs[7]);
	
	MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	MPI_Type_free(&slice1);
	MPI_Type_free(&slice2);
}

void exchange3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, int nbrtop, int nbrbot, int nbrleft, int nbrright, int nbrfront, int nbrback, MPI_Comm comm){
        MPI_Request reqs[12];
        MPI_Datatype slice1;
        //MPI_Type_vector((e-s)*z, 1, z, MPI_INT, &slice1);
        MPI_Type_vector((e-s-1/*+1*/), z, y*z, MPI_INT, &slice1);
        MPI_Type_commit(&slice1);

        //Send to left and right nbrs
        MPI_Irecv(&g[s+1][modulo(e2+1,x)][s3], 1, slice1, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s+1][modulo(s2-1,x)][s3], 1, slice1, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s+1][modulo(s2,x)][s3], 1, slice1, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s+1][modulo(e2,x)][s3], 1, slice1, nbrright, 1, comm, &reqs[3]);

        //Send to top and bot nbrs
        MPI_Datatype slice2;
        MPI_Type_vector((e2-s2+1)*z, 1, 1, MPI_INT, &slice2);
        MPI_Type_commit(&slice2);

        MPI_Irecv(&g[modulo(e+1,x)][s2][s3], 1, slice2, nbrtop, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,x)][s2][s3], 1, slice2, nbrbot, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,x)][s2][s3], 1, slice2, nbrbot, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,x)][s2][s3], 1, slice2, nbrtop, 3, comm, &reqs[7]);

	//Send front and back nbrs
	MPI_Datatype slice3;
	MPI_Datatype row1;
	MPI_Type_vector((e-s-1)/*(e2-s2-1)*/, 1, z, MPI_INT, &row1);
	MPI_Type_commit(&row1);
	MPI_Type_vector((e2-s2-1), 1, y*z, row1, &slice3);
	MPI_Type_commit(&slice3);

	MPI_Irecv(&g[s][s2][modulo(e3+1,x)], 1, slice2, nbrfront, 2, comm, &reqs[8]);
	MPI_Irecv(&g[s][s2][modulo(s3-1,x)], 1, slice2, nbrback, 2, comm, &reqs[9]);
	MPI_Isend(&g[s][s2][s3], 1, slice2, nbrback, 2, comm, &reqs[10]);
	MPI_Isend(&g[s][s2][e3], 1, slice2, nbrfront, 2, comm, &reqs[11]);

        MPI_Waitall(12, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&slice1);
        MPI_Type_free(&slice2);
	MPI_Type_free(&slice3);
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
