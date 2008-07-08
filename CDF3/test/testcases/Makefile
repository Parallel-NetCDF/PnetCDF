srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) 

CC		= /home/kgao/MPICH/bin/mpicc
CFLAGS          = -g -O2 $(INCLUDES) 

SRCS		= $(srcdir)/ncmpi_vars_null_stride.c
OBJS		= $(SRCS:.c=.o)
PROGS		= ncmpi_vars_null_stride
LDFLAGS		= -L$(LIBDIR) 
LIBS		= -lpnetcdf 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

ncmpi_vars_null_stride: ncmpi_vars_null_stride.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
clean:
		rm -f $(OBJS) $(PROGS) *.nc

test: ncmpi_vars_null_stride
	./ncmpi_vars_null_stride test.nc
	@echo '*** Success ***'

testing: test
