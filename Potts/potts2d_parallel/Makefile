CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean
all: potts2d_parallel potts2d_parallel_2d

potts2d_parallel: potts2d_parallel.c decomp1d.o modulo.o hamiltonian.o exchange.o observables.o metropolis.o mpi_io.o
	$(MPICC) $(CCFLAGS) -o $@ $^
potts2d_parallel_2d: potts2d_parallel_2d.c decomp1d.o modulo.o hamiltonian.o exchange.o observables.o metropolis.o mpi_io.o
	$(MPICC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts2d_parallel potts2d_parallel_2d
