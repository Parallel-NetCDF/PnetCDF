all: pres_temp_4D_wr

pres_temp_4D_wr: pres_temp_4D_wr.o 
	mpicc -o pres_temp_4D_wr pres_temp_4D_wr.o -L../../src/lib -L../common -lpnetcdf -ltestutils   
pres_temp_4D_wr.o: pres_temp_4D_wr.c 
	mpicc -c -g -O2 -I./../../src/lib -I. pres_temp_4D_wr.c


clean:
	rm *.o pres_temp_4D_wr
