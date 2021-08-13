CC = gcc
MPICC = mpicc
CCFLAGS = -lm -Wall -W -lgsl -lgslcblas

.PHONY: clean
all: ising2d_parallel ising2d

#ising2d_parallel_subarray: ising2d_parallel_subarray.c decomp1d.o
#	$(MPICC) $(CCFLAGS) -o $@ $^
ising2d_parallel: ising2d_parallel.c decomp1d.o
	$(MPICC) $(CCFLAGS) -o $@ $^
ising2d: ising2d.c
	$(CC) $(CCFLAGS) -o $@ ising2d.c
clean:
	$(RM) *.o *.txt ising2d_parallel ising2d
