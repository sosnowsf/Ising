CC = gcc
#MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

.PHONY: clean
#all: potts2d_serial potts2d_parallel potts3d_serial potts3d_parallel

#potts2d_parallel: potts2d_parallel.c decomp1d.o
#	$(MPICC) $(CCFLAGS) -o $@ $^
#potts2d_serial: potts2d_serial.c
#	$(CC) $(CCFLAGS) -o $@ $^
#potts3d_parallel: potts3d_parallel.c decomp1d.o
#	$(MPICC) $(CCFLAGS) -o $@ $^
potts3d_serial: potts3d_serial.c hamiltonian.o modulo.o observables.o metropolis.o
	$(CC) $(CCFLAGS) -o $@ $^
clean:
	$(RM) *.o *.txt potts3d_serial
