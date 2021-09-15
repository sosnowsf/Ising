#include"grid_size.h"
#include"hamiltonian.h"
#include"modulo.h"

void metropolis_sweep(int g[m][n], int, const double , double, gsl_rng *);
void metropolis_sweep_triangular(int g[m][n], int, const double , double, gsl_rng *);
void metropolis_sweep_hexagonal(int g[m][n], int, const double , double, gsl_rng *);
