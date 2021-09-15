#include"grid_size.h"

double energy_2d(int g[m][n], const double, int, int, int, int);
double energy_triangular_2d(int g[m][n], const double, int, int, int, int);
double energy_hexagonal_2d(int g[m][n], const double, int, int, int, int);
double magnetisation_2d(int g[m][n], int, int, int, int);
double var_mag_2d(int g[m][n], int, int, int, int, double);
double var_enrg_2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_tri_2d(int g[m][n], int, int, int, int, double, double);
double var_enrg_hex_2d(int g[m][n], int, int, int, int, double, double);
double energy(int g[m][n], const double, int, int);
double energy_triangular(int g[m][n], const double, int, int);
double energy_hexagonal(int g[m][n], const double, int, int);
double magnetisation(int g[m][n], int, int);
double var_mag(int g[m][n], int, int, double);
double var_enrg(int g[m][n], int, int, double, double);
double var_enrg_tri(int g[m][n], int, int, double, double);
double var_enrg_hex(int g[m][n], int, int, double, double);

