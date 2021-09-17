CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -W -lgsl -lgslcblas

.PHONY: clean
all:  ising2d_parallel_1d ising2d_parallel_2d ising2d

ising2d_parallel_1d: ising2d_parallel_1d.c decomp1d.o modulo.o mpi_io.o exchange.o observables.o metropolis.o
	$(MPICC) $(CCFLAGS) -o $@ $^
ising2d_parallel_2d: ising2d_parallel_2d.c decomp1d.o modulo.o mpi_io.o exchange.o observables.o metropolis.o
	$(MPICC) $(CCFLAGS) -o $@ $^
ising2d: ising2d.c
	$(CC) $(CCFLAGS) -o $@ ising2d.c
clean:
	$(RM) *.o *.txt ising2d ising2d_parallel_1d ising2d_parallel_2d
