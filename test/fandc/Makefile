ALL: pnctestf pnctest csnap pnf_test

INCDIR = /home/gropp/projects/parallel-netcdf/src/lib
INCDIRF = /home/gropp/projects/parallel-netcdf/src/libf
LNKDIR = /home/gropp/projects/parallel-netcdf/src/lib

MPICC = /home/gropp/tmp/mpi2-inst/bin/mpicc
MPIF77 = /home/gropp/tmp/mpi2-inst/bin/mpif77

EXECS = pnctestf pnctest csnap

pnctestf: pnctestf.F
	$(MPIF77) -o pnctestf pnctestf.F -I$(INCDIRF) -L$(LNKDIR) -lpnetcdf -lm

pnctest: pnctest.c
	$(MPICC) -o pnctest pnctest.c -I$(INCDIR) -L$(LNKDIR) -lpnetcdf -lm

csnap:  csnap.o
	$(MPICC) -L$(LNKDIR) -o csnap csnap.o -lpnetcdf -lm

csnap.o:  csnap.c
	$(MPICC) -I$(INCDIR) -c csnap.c

pnf_test:  pnf_test.o
	$(MPIF77) -L$(LNKDIR) -o pnf_test pnf_test.o -lpnetcdf -lm

pnf_test.o:  pnf_test.F
	$(MPIF77) -I$(INCDIRF) -c pnf_test.F

clean:
	rm -f $(EXECS) *.o *.nc

