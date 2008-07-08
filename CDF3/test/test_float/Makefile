srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common

CC		= /home/kgao/MPICH/bin/mpicc
CFLAGS          = -g -O2 $(INCLUDES) 

SRCS		= test_read.c test_read_indep.c test_write.c test_write_indep.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_read test_read_indep test_write test_write_indep
LDFLAGS		= -L$(LIBDIR) -L..//common 
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_read: test_read.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_read_indep: test_read_indep.o 
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_write: test_write.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_write_indep: test_write_indep.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
		rm -f $(OBJS) $(PROGS) *.nc
