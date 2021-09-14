#include"grid_size.h"

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
