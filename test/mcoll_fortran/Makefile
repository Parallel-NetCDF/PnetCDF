srcdir = .


include ../../macros.make

ALL: mcoll_testf 

INCDIR = $(srcdir)/../../src/lib
INCDIRF = ../../src/libf/
LNKDIR = ../../src/lib

MPICC = /usr/local/bin/mpicc
MPIF77 = /usr/local/bin/mpif77

EXECS = mcoll_testf 

mcoll_testf: mcoll_testf.F
	$(LINK.F) -o mcoll_testf $(srcdir)/mcoll_testf.F -I$(INCDIRF) -L$(LNKDIR) -lpnetcdf -lm $(LIBS)

clean:
	rm -f $(EXECS) *.o *.nc

