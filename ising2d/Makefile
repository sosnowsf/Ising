CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall 

ising2d_parallel: ising2d_parallel.c decomp1d.o
	$(MPICC) $(CCFLAGS) -o $@ $^
ising2d: ising2d.c
	$(CC) $(CCFLAGS) -o $@ ising2d.c
