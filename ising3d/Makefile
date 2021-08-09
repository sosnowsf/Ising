CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

ising3d_parallel: ising3d_parallel.c decomp1d.o
	$(MPICC) $(CCFLAGS)  -o $@ $^
ising3d: ising3d.c
	$(CC) $(CCFLAGS) -o $@ ising3d.c
