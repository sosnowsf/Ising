#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define m 100
#define n 100

void init_grid(int g[m][n], int, gsl_rng *);
void metropolis_sweep(int g[m][n], const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], const double , double, gsl_rng *);
double energy(int g[m][n], const double);
double energy2(int g[m][n], const double);
double energy_triangular(int g[m][n], const double);
double energy_triangular2(int g[m][n], const double);
double energy_hexagonal(int g[m][n], const double);
double energy_hexagonal2(int g[m][n], const double);
double magnetisation(int g[m][n]);
double magnetisation2(int g[m][n]);
void print_matrix(int g[m][n]);
void write_stats(char *, double *, double *, double *, double *, double *, int);
void save_configuration(char *, int g[m][n]);
void reset_grid(int G[m][n], int g[m][n]);
int modulo(int, int);
double var_mag(int g[m][n], double);
double var_enrg(int g[m][n], double, double);
double var_enrg_tri(int g[m][n], double, double);
double var_enrg_hex(int g[m][n], double, double);
void write_stats_autocorr(char *title, double *autocorr, int r2);
double autocorrelation(int k, int N, double *x);

int main(int argc, char *argv[]){

	//Parse variables
	int sq=0;
	int triangle=0;
	int hex=0;
	int c=0;
	int opt;
	while((opt = getopt(argc, argv, "123ch")) != -1){
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
			case'c':
				c=1;
				break;
			case 'h':
				printf("Command lines:\n -1 for square grid\n -2 for triangular grid \n -3 for hexagonal grid\n -c for cold start(hot start default)\n");
				return 0;
		}
	}
	//Error checking
	if((sq+triangle+hex)==0){
		perror("No parse command given, use -h for help");
		exit(EXIT_FAILURE);
	}
	if((sq+triangle+hex)>1){
		perror("Please provide only one parse command at a time");
		exit(EXIT_FAILURE);
	}
	double T;
	const double J = 1.0; //Coupling constant
	//const double k = pow(1.38064852, -23); //Boltzmann constant
	
	double r1 = 50.0; //Number of times simulation is ran to ge the average
	int r2 = 100; //Number of temperatures
	int r3 = 1000;//1000; //Number of sweeps

	//int seed = 1997;
	//srand48(seed);
	int g[m][n]; //define grid
	
	/*Set up rng*/
	gsl_rng *gsl_mt = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long seed = 1999;
	gsl_rng_set(gsl_mt, seed);

	init_grid(g,c,gsl_mt);
	char title1[100];
	snprintf(title1, 100, "initial_configuration.txt");
	save_configuration(title1, g);

	//Allocate memroy for observables
	double *mag, *enrgy, *spec, *sus, *temp;
	mag = malloc(r2*sizeof(double));
	enrgy = malloc(r2*sizeof(double));
	temp = malloc(r2*sizeof(double));
	spec = malloc(r2*sizeof(double));
	sus = malloc(r2*sizeof(double));	

	//Arrays for autocorrelation
	double *mag2, *autocorr;
	mag2 = malloc(r3*sizeof(double));
	autocorr = malloc(r3*sizeof(double));

	//Initiate to zero
	for(int i=0; i<r2; i++){
		mag[i]=0.0;
		enrgy[i]=0.0;
		spec[i]=0.0;
		sus[i]=0.0;
	}

	//Triangular case
        if(triangle==1){
        for(int k=0; k<r1; k++){
                T=0.0001;
		double t; 
                for(int j=0; j<r2; j++){
			t=1/T; //Precalculate inverse
			//Perform sweeps for given temperature
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_triangular(g,J,t,gsl_mt);
                        }
			//Save configuration after a run
                        /*if(k==0){
                                char title2[100];
                                snprintf(title2, 100, "final_configuration.txt");
                                save_configuration(title2,g);
                        }*/
			//Calculate observables
                        double m1 = magnetisation(g);
                        double e1 = energy_triangular(g,J);
                        mag[j] += m1/r1;
                        enrgy[j] += e1/r1;
                        spec[j] += var_enrg_tri(g,J,e1)*t*t/r1;
                        sus[j] += var_mag(g,m1)*t/r1;
                        if(k==0) temp[j] = T;
                        T+=0.07;
                        init_grid(g,c,gsl_mt); //Reset grid
                }
        }
        }

	//Hexagonal case
        if(hex==1){
        for(int k=0; k<r1; k++){
                T=0.0001;
		double t;
                for(int j=0; j<r2; j++){
			t=1/T;
                        for(int i=0; i<r3; i++){
                                metropolis_sweep_hexagonal(g,J,t,gsl_mt);
                        }
                        /*if(k==0){
                                char title2[100];
                                snprintf(title2, 100, "final_configuration.txt");
                                save_configuration(title2,g);
                        }*/
                        double m1 = magnetisation(g);
                        double e1 = energy_hexagonal(g,J);
                        mag[j] += m1/r1;
                        enrgy[j] += e1/r1;
                        spec[j] += var_enrg_hex(g,J,e1)*t*t/r1;
                        sus[j] += var_mag(g,m1)*t/r1;
                        if(k==0) temp[j] = T;
                        T+=0.05;
                        init_grid(g,c,gsl_mt);
                }
        }
        }

	//Square lattice case
	if(sq==1){
	for(int k=0; k<r1; k++){
		T=0.0001;
		double t;
		for(int j=0; j<r2; j++){
			t=1/T;
			for(int i=0; i<r3; i++){
				metropolis_sweep(g,J,t,gsl_mt);
				//Measure autocorrelation 
				/*if(k==0 && (T==0.000100 || T==2.250100 || T==4.5001)){
					mag2[i] = magnetisation(g);
					autocorr[i] = mag2[i];//autocorrelation(i,r3,mag2);
					//snprintf(title, 100, "Autocorrelations.txt");
					//write_stats_autocorr(title, temp, autocorr, r2);
				}*/
				
			}
			//Autocorrelation code for different temperatures
			//if(k==0 && (T==2.2501 || T==4.2501)){
			if(k==0){
				//char title[100];
				//snprintf(title, 100, "Autocorrelations_T%f.txt", T);
				//write_stats_autocorr(title, autocorr, r3);
				printf("%f\n", T);
				char title2[100];
				snprintf(title2, 100, "final_configuration_T%f.txt", T);
				save_configuration(title2,g);
			}
                        double m1 = magnetisation(g);
                        double e1 = energy(g,J);
                        mag[j] += m1/r1;
                        enrgy[j] += e1/r1;
                        spec[j] += var_enrg(g,J,e1)*t*t/r1;
                        sus[j] += var_mag(g,m1)*t/r1;
                        if(k==0) temp[j] = T;
                        T+=0.05;
			//Code to save configurations for some metastable states
			/*int a=0;
			if((m1<0.9) && (T<2.0) && (a==0)){
				char title[100];
				snprintf(title, 100, "final_configuration_metastable_T%f_M%f.txt", T, m1);
				save_configuration(title,g);
				a++;
			}*/
                        init_grid(g,c,gsl_mt);
		}
	}
	}
	//Save observables and free memory
	char title[100];
	if(hex==1 && c==1) snprintf(title, 100, "stats_hexagonal_serial_cold.txt");
	if(triangle==1 && c==1) snprintf(title, 100, "stats_triangular_serial_cold.txt");
	if(sq==1 && c==1) snprintf(title, 100, "stats_square_serial_cold.txt");
        if(hex==1 && c==0) snprintf(title, 100, "stats_hexagonal_serial_hot.txt");
        if(triangle==1 && c==0) snprintf(title, 100, "stats_triangular_serial_hot.txt");
        if(sq==1 && c==0) snprintf(title, 100, "stats_square_serial_hot.txt");
	write_stats(title,mag,enrgy,sus,spec,temp,r2);
	free(mag);
	free(enrgy);
	free(temp);
	free(spec);
	free(sus);
	free(mag2);
	free(autocorr);
return 0;
}

//Initialise grid for both hot and cold starts
void init_grid(int g[m][n], int c, gsl_rng *gsl_mt){
	int i,j;
	//hot start
	if(c==1){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				g[i][j]=1;
			}
		}
	}
	//cold start
	if(c==0){
		for(i=0; i<m; i++){
			for(j=0; j<n; j++){
				if(gsl_rng_uniform(gsl_mt)<0.5) g[i][j]=-1;
				else{ 
					g[i][j]=1;
				}
			}
		}
	}
}


//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], const double J, double T, gsl_rng *gsl_mt){
	int i,j;
	double dE;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
		dE = 2.0*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]); //Energy difference if spin is flipped
		if(dE<=0) g[i][j] *=-1;
		else{ if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] *=-1; }
		}
	}
}

//Sweep for triangular lattice
void metropolis_sweep_triangular(int g[m][n], const double J, double T, gsl_rng *gsl_mt){
        int i,j;
	double dE;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
                        dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        if(dE<=0) g[i][j] *= -1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] *= -1;
                        }
                }
        }
}

//Sweep for hexagonal lattice
void metropolis_sweep_hexagonal(int g[m][n], const double J, double T, gsl_rng *gsl_mt){
        int i,j;
	double dE;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        //Calculate the energy difference and check if spin has to be flipped
	 		dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        if(dE<=0) g[i][j] *= -1;
                        else{
                                if(gsl_rng_uniform(gsl_mt)<(exp(-dE*T))) g[i][j] *= -1;
                        }
                }
        }
}


//Calculate the energy per site of the system for a square lattice
double energy(int g[m][n], const double J){
	int i,j;
	double E=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
		}
	}
return E/(2.0*m*n);
}

//Calculate the energy variance per site given the average energy for a square lattice
double var_enrg(int g[m][n], double J, double avg){
	int i,j;
	double dE;
	double var=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
			dE/=2.0;
			var += (dE - avg)*(dE - avg);
		}
	}
return var/(m*n);
}

//Calculate the energy per site of the system for a triangular lattce
double energy_triangular(int g[m][n], const double J){
        int i,j;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                }
        }
return E/(2.0*m*n);
}

//Calculate the average energy variance for triangular lattice
double var_enrg_tri(int g[m][n], double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
			dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)] + g[modulo(i+1,m)][modulo(j+((int)pow(-1,i)),n)] +g[modulo(i-1,m)][modulo(j+((int)pow(-1,i)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/(m*n);
}

//Calculate the average magnetisation varaince
double var_mag(int g[m][n], double avg){
	int i,j;
	double var=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			var += (g[i][j] - avg)*(g[i][j] - avg);
		}
	}
return var/(m*n);
}


//Calculate the energy per site of the system for a hexagonal lattice
double energy_hexagonal(int g[m][n], const double J){
        int i,j;
        double E=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
                        E+=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                }
        }
return E/(2.0*m*n);
}

//Calculate the average enrgy variance for a hexagonal lattice
double var_enrg_hex(int g[m][n], double J, double avg){
        int i,j;
        double dE;
        double var=0;
        for(i=0; i<m; i++){
                for(j=0; j<n; j++){
			dE=-J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+((int)pow(-1,i+j)),n)]);
                        dE/=2.0;
                        var += (dE - avg)*(dE - avg);
                }
        }
return var/(m*n);
}

//Calculate magnetisation per site
double magnetisation(int g[m][n]){
	int i,j;
	double M=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			M+=g[i][j];
		}
	}
return fabs(M/(m*n));
}

//Print out the matrix
void print_matrix(int g[m][n]){
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%d ", g[i][j]);
		}
	printf("\n");
	}
}

//Reset grid after every temp increase
void reset_grid(int G[m][n], int g[m][n]){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			g[i][j]=G[i][j];
		}
	}
}

//Write stats into file
void write_stats(char *title, double *mag, double *energy, double *sus, double *spec, double *T, int r2){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r2; i++){
		fprintf(fp,"%f %f %f %f %f\n", T[i], mag[i], energy[i], sus[i], spec[i]);
	}
	fclose(fp);
}

//Write autocorrelation data into file
void write_stats_autocorr(char *title, double *autocorr, int r3){
	FILE *fp = fopen(title, "w");
	for(int i=0; i<r3; i++){
			fprintf(fp, "%f \n", autocorr[i]);
	}
	fclose(fp);
}

//Save system state into file
void save_configuration(char *title, int g[m][n]){
	FILE *fp = fopen(title, "w");
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			fprintf(fp, "%d\t", g[i][j]);
		}
		fprintf(fp, "\n");
	}
}

//Modulo function for handling -1(mod)n
int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}

//Calculate autocorrelation
double autocorrelation(int t, int N, double *x){
	double x1=0;
	double x2=0;
	for(int i=0; i<(N-t-1); i++){
		x2+=x[i];
	}
	x2/=(N-t-1);
	for(int j=0; j<(N-t-1); j++){
		x1 += x[j+t]*(x[j] - x2);
	}
	x1/=(N-t);
return x1;
}

