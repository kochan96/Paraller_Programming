CC	= gcc   
# uncomment for paraller
# CFLAGS    = -fopenmp
LIBS	  = -lm

all: main

main: main.c
	$(CC) $(CFLAGS) -o main main.c $(LIBS)

clean: 
	rm -f *.o main