ALL: pnctestf pnctest

INCDIR = /home/gropp/projects/parallel-netcdf/src/lib
INCDIRF = /home/gropp/projects/parallel-netcdf/src/libf
LNKDIR = /home/gropp/projects/parallel-netcdf/src/lib

MPICC = /home/gropp/tmp/mpi2-inst/bin/mpicc
MPIF77 = /home/gropp/tmp/mpi2-inst/bin/mpif77

pnctestf: pnctestf.F
	$(MPIF77) -o pnctestf pnctestf.F -I$(INCDIRF) -L$(LNKDIR) -lpnetcdf -lm

pnctest: pnctest.c
	$(MPICC) -o pnctest pnctest.c -I$(INCDIR) -L$(LNKDIR) -lpnetcdf -lm

clean:
	rm -f pnctestf pnctest *.o *.nc

