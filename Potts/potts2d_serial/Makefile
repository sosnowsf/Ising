CC = gcc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean

potts2d_serial: potts2d_serial.c hamiltonian.o modulo.o metropolis.o observables.o
	$(CC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts2d_serial 
