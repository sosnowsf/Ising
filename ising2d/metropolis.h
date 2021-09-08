#include"grid_size.h"
#include"modulo.h"

void metropolis_sweep(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep_triangular2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], int, int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal2d(int g[m][n], int, int, int, int, const double , double, gsl_rng *);
