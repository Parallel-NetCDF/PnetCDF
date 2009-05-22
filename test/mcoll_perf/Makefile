srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR)

CC		= /home/kgao/mpich2-1.0.8/bin/mpicc
CC		= /bgsys/drivers/ppcfloor/comm/bin/mpicc 
CFLAGS          = -g $(INCLUDES) 
SRCS		= mcoll_perf.c single_mcoll_perf_record.c
OBJS		= $(SRCS:.c=.o)
PROGS		= mcoll_perf single_mcoll_perf_record
LDFLAGS		= -L$(LIBDIR) 
LIBS		= -lpnetcdf  


.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

mcoll_perf: mcoll_perf.o
	$(CC) -g -o $@ $@.o $(LDFLAGS) $(LIBS)
single_mcoll_perf_record: single_mcoll_perf_record.o
	$(CC) -g -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
		rm -f $(OBJS) $(PROGS) *.nc
