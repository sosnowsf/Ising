CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

all: potts3d_parallel potts3d_parallel_2d potts3d_parallel_3d
.PHONY: clean

potts3d_parallel: potts3d_parallel.c decomp1d.o hamiltonian.o modulo.o exchange.o metropolis.o observables.o mpi_io.o 
	$(MPICC) $(CCFLAGS) -o $@ $^
potts3d_parallel_2d: potts3d_parallel_2d.c decomp1d.o hamiltonian.o modulo.o exchange.o metropolis.o observables.o mpi_io.o
	$(MPICC) $(CCFLAGS) -o $@ $^
potts3d_parallel_3d: potts3d_parallel_3d.c decomp1d.o hamiltonian.o modulo.o exchange.o metropolis.o observables.o mpi_io.o
	$(MPICC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts3d_parallel potts3d_parallel_2d potts3d_parallel_3d
