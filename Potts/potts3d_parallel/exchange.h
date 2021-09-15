#include"grid_size.h"
#include"modulo.h"
#include"hamiltonian.h"

void exchange(int g[x][y][z], int, int, int, int, MPI_Comm comm);
void exchange2d(int g[x][y][z], int, int, int, int, int, int, int, int, MPI_Comm);
void exchange3d(int g[x][y][z], int, int, int, int, int, int, int, int, int, int, int, int, MPI_Comm);
