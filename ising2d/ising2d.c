#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<mpi.h>

#define m 50
#define n 50

void init_grid(int g[m][n]);
void metropolis_sweep(int g[m][n], const double , double, const double);
double energy(int g[m][n], const double);
double magnetisation(int g[m][n]);
void print_matrix(int g[m][n]);
void write_stats(double *, double *, double *, int);
void reset_grid(int G[m][n], int g[m][n]);
int modulo(int, int);

int main(){
	double T = 0; //Start at 0 degrees Celsius
	const double J = 1.0; //Coupling constant
	const double k = pow(1.38064852, -23); //Boltzmann constant
	int r1 = 100;
	int r2 = 100;

	int seed = 1997;
	srand48(seed);
	int g[m][n]; //define grid
	
	init_grid(g);

	//print_matrix(g);
	printf("%f\n", magnetisation(g));

	double *mag, *enrgy, *temp;
	mag = malloc(100*sizeof(double));
	enrgy = malloc(100*sizeof(double));
	temp = malloc(100*sizeof(double));
	
	for(int j=0; j<r1; j++){
		for(int i=0; i<r2; i++){
			metropolis_sweep(g,J,T,k);
			//energy(g,m,n);
			//printf("%f\n", magnetisation(g));
		}
		mag[j] = magnetisation(g);
		enrgy[j] = energy(g,J);
		temp[j] = T;
		T+=0.05;
		//print_matrix(g);
		init_grid(g);
	}
	write_stats(mag,enrgy,temp,r1);
	//printf("%f\n", magnetisation(g));
	//g[3][5]=9;
	//print_matrix(g);
	free(mag);
	free(enrgy);
	free(temp);
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
	//print_matrix(g);
}

//Sweep through the lattice using the metropolis algorithm
void metropolis_sweep(int g[m][n], const double J, double T, const double k){
	int i,j;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			//Calculate the energy difference and check if spin has to be flipped
			double dE = 2*J*g[i][j]*(g[modulo(i+1,m)][j] + g[modulo(i-1,m)][j] + g[i][modulo(j+1,n)] + g[i][modulo(j-1,n)]);
			//double dE = J*g[i][j]*(g[(i+1)%m][j]+g[(i-1)&m][j]+g[i][(j+1)%n]+g[i][(j-1)%n]);
			if(dE<0) g[i][j] *= -1;
			else{
				if(drand48()<(exp(-dE/T))){
					//printf("%f %f %f\n", dE, T, exp(-dE/T));
					g[i][j] *= -1;
				}
			}
		}
	}
}

//Calculate the energy of the system
double energy(int g[m][n], const double J){
	int i,j;
	double E=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			E+=-J*g[i][j]*(g[modulo(i+1,m)][j]+g[modulo(i-1,m)][j]+g[i][modulo(j+1,n)]+g[i][modulo(j-1,n)]);
		}
	}
return E;
}

//Calculate magnetisation persite
double magnetisation(int g[m][n]){
	int i,j;
	double M=0;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			M+=g[i][j];
		}
	}
return fabs(M)/(m*n);
}

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

void write_stats(double *mag, double *energy, double *T, int r1){
	FILE *fp = fopen("stats_serial.txt", "w");
	for(int i=0; i<r1; i++){
		fprintf(fp,"%f %f %f\n", T[i], mag[i], energy[i]);
	}
	fclose(fp);
}

int modulo(int a, int b){
	int r = a%b;
	if(r<0) r+=b;
return r;
}
