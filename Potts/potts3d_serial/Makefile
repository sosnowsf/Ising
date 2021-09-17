CC = gcc
#MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean

potts3d_serial: potts3d_serial.c hamiltonian.o modulo.o observables.o metropolis.o
	$(CC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts3d_serial
