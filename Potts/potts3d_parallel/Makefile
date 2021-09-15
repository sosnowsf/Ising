CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean

potts3d_parallel: potts3d_parallel.c decomp1d.o hamiltonian.o modulo.o exchange.o metropolis.o observables.o 
	$(MPICC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts3d_parallel
