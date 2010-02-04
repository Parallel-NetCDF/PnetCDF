srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common -I../../src/lib

CC		= /usr/local/bin/mpicc
CFLAGS          = -g -O2 $(INCLUDES) 

SRCS		= test_inq_format.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_inq_format
LDFLAGS		= -L$(LIBDIR) -L..//common 
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_inq_format: test_inq_format.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
clean:
	rm -f $(OBJS) $(PROGS) *.nc

distclean: clean
	rm -f Makefile
