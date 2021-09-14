CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean
all: ising3d_parallel ising3d

ising3d_parallel: ising3d_parallel.c decomp1d.o modulo.o exchange.o observables.o mpi_io.o metropolis.o
	$(MPICC) $(CCFLAGS)  -o $@ $^
ising3d: ising3d.c
	$(CC) $(CCFLAGS) -o $@ ising3d.c
clean:
	$(RM) *.o *.txt ising3d ising3d_parallel
