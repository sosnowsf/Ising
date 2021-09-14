#include"grid_size.h"

void metropolis_sweep(int g[x][y][z], int, int, const double , double, gsl_rng *);
void metropolis_sweep2d(int g[x][y][z], int, int, int, int, const double , double, gsl_rng *);
void metropolis_sweep3d(int g[x][y][z], int, int, int, int, int, int, const double , double, gsl_rng *);
