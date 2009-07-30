srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common

CC		= /usr/local//bin/mpicc
CFLAGS          = -g -O2 $(INCLUDES) 

SRCS		= test_check_header.c test_check_header1.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_check_header test_check_header1
LDFLAGS		= -L$(LIBDIR) -L..//common 
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_check_header: test_check_header.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

test_check_header1: test_check_header1.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
		rm -f $(OBJS) $(PROGS) *.nc
