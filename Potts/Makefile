CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -lgsl -lgslcblas

potts2d_serial: potts2d_serial.c
	$(CC) $(CCFLAGS) -o $@ $^
