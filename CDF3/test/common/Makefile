srcdir = .


CC = /home/kgao/MPICH/bin/mpicc
CFLAGS = -g -O2

OBJS = testutils.o

all: libtestutils.a

libtestutils.a: $(OBJS)
	ar cru libtestutils.a $(OBJS)
testutils.o: testutils.c testutils.h

clean:
	rm -f *.o *.a
