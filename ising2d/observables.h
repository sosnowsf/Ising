#include"grid_size.h"
#include"modulo.h"

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
double var_mag(int g[m][n], int, int, double);
double var_mag2d(int g[m][n], int, int, int, int, double);
double var_enrg(int g[m][n], int, int, double, double);
double var_enrg2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_tri(int g[m][n], int, int, double, double);
double var_enrg_tri2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_hex(int g[m][n], int, int, double, double);
double var_enrg_hex2d(int g[m][n], int, int, int, int, double, double);
