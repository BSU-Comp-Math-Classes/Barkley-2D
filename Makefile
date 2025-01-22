CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

barkley: barkley.c
	$(CC) $(CFLAGS) -o barkley barkley.c

clean:
	rm -f barkley

