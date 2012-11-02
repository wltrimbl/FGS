CC=gcc
CFLAGS=-O3 -Wall
SRCS=util_lib.c hmm_lib.c run_hmm.c
OBJ=$(SRCS:.c=.o)

all: fgs

#hmm.obj: $(SRCS)
#	$(CC) -c $(CFLAGS) $(SRCS)

fgs:	$(OBJ)
	$(CC) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o  -lm 

.o:
	$(CC) -c $(CFLAGS) $(SRCS)

clean:
	rm -rf *.o FragGeneScan* *~
