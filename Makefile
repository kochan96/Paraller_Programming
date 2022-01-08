CC	= gcc   
# uncomment for paraller
CFLAGS    = -fopenmp
LIBS	  = -lm

all: main

main: main.c
	$(CC) $(CFLAGS) -o qr_factorization main.c common.c $(LIBS)

clean: 
	rm -f *.o qr_factorization