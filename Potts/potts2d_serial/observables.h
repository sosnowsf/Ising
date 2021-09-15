#include"grid_size.h"
#include"hamiltonian.h"
#include"modulo.h"

double energy(int g[m][n], const double);
double energy2(int g[m][n], const double);
double energy_triangular(int g[m][n], const double);
double energy_triangular2(int g[m][n], const double);
double energy_hexagonal(int g[m][n], const double);
double energy_hexagonal2(int g[m][n], const double);
double magnetisation(int g[m][n]);
double magnetisation2(int g[m][n]);
double var_mag(int g[m][n], double);
double var_enrg(int g[m][n], double, double);
double var_enrg_tri(int g[m][n], double, double);
double var_enrg_hex(int g[m][n], double, double);
