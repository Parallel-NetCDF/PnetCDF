## Test Installed PnetCDF

This folder contains a Makefile and run scripts to test an already installed
PnetCDF library, using a subset of test programs from this software
distribution. All test programs are designed to run on 4 MPI processes.

### Usage
* First, set the following four environment variables.
  + `CC` - MPI C compiler
  + `FC` - MPI Fortran compiler
  + `RUN_CMD` - Run command, e.g. mpiexec or srun
  + `OUT_DIR` - Folder to store the output files
  + `PNETCDF_DIR` - Installation of PnetCDF library

  Below shows an example of commands used on Perlmutter at NERSC.
  ```console
  module load cray-parallel-netcdf
  export CC=cc
  export FC=ftn
  export RUN_CMD=srun
  export OUT_DIR=$SCRATCH
  export PNETCDF_DIR=$PE_PARALLEL_NETCDF_DIR
  ```
* Run command `make` to compile the test programs and generate run script files.
  For instance,
  ```console
  make -j 8
  ```
  <details>
  <summary>Example output of make command shown on screen (click to expand)</summary>

  ```console
  login02::test_installed(1:00pm) #1096 make -j 8
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -I../common -c ../common/testutils.c -o testutils.o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/utils.F90 -o utils.o
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../C/pres_temp_4D_wr.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  pres_temp_4D_wr.o  testutils.o -lpnetcdf -o pres_temp_4D_wr
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../C/pres_temp_4D_rd.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  pres_temp_4D_rd.o  testutils.o -lpnetcdf -o pres_temp_4D_rd
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../header/header_consistency.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  header_consistency.o  testutils.o -lpnetcdf -o header_consistency
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/flexible_bput.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  flexible_bput.o  testutils.o -lpnetcdf -o flexible_bput
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/interleaved.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  interleaved.o  testutils.o -lpnetcdf -o interleaved
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/req_all.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  req_all.o  testutils.o -lpnetcdf -o req_all
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/i_varn_indef.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  i_varn_indef.o  testutils.o -lpnetcdf -o i_varn_indef
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/large_num_reqs.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  large_num_reqs.o  testutils.o -lpnetcdf -o large_num_reqs
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/test_bput.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_bput.o  testutils.o -lpnetcdf -o test_bput
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/i_varn_int64.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  i_varn_int64.o  testutils.o -lpnetcdf -o i_varn_int64
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/mcoll_perf.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  mcoll_perf.o  testutils.o -lpnetcdf -o mcoll_perf
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../nonblocking/wait_after_indep.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  wait_after_indep.o  testutils.o -lpnetcdf -o wait_after_indep
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/add_var.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  add_var.o  testutils.o -lpnetcdf -o add_var
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/alignment_test.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  alignment_test.o  testutils.o -lpnetcdf -o alignment_test
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/buftype_free.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  buftype_free.o  testutils.o -lpnetcdf -o buftype_free
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/check_striping.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  check_striping.o  testutils.o -lpnetcdf -o check_striping
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/check_type.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  check_type.o  testutils.o -lpnetcdf -o check_type
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/collective_error.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  collective_error.o  testutils.o -lpnetcdf -o collective_error
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/flexible.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  flexible.o  testutils.o -lpnetcdf -o flexible
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/flexible2.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  flexible2.o  testutils.o -lpnetcdf -o flexible2
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/flexible_varm.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  flexible_varm.o  testutils.o -lpnetcdf -o flexible_varm
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/inq_num_vars.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  inq_num_vars.o  testutils.o -lpnetcdf -o inq_num_vars
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/inq_recsize.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  inq_recsize.o  testutils.o -lpnetcdf -o inq_recsize
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/ivarn.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  ivarn.o  testutils.o -lpnetcdf -o ivarn
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/large_var_cdf5.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  large_var_cdf5.o  testutils.o -lpnetcdf -o large_var_cdf5
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/last_large_var.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  last_large_var.o  testutils.o -lpnetcdf -o last_large_var
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/mix_collectives.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  mix_collectives.o  testutils.o -lpnetcdf -o mix_collectives
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/modes.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  modes.o  testutils.o -lpnetcdf -o modes
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/ncmpi_vars_null_stride.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  ncmpi_vars_null_stride.o  testutils.o -lpnetcdf -o ncmpi_vars_null_stride
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/noclobber.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  noclobber.o  testutils.o -lpnetcdf -o noclobber
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/nonblocking.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  nonblocking.o  testutils.o -lpnetcdf -o nonblocking
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/one_record.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  one_record.o  testutils.o -lpnetcdf -o one_record
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/record.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  record.o  testutils.o -lpnetcdf -o record
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/redef1.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  redef1.o  testutils.o -lpnetcdf -o redef1
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/scalar.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  scalar.o  testutils.o -lpnetcdf -o scalar
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_erange.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_erange.o  testutils.o -lpnetcdf -o test_erange
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_fillvalue.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_fillvalue.o  testutils.o -lpnetcdf -o test_fillvalue
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_get_varn.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_get_varn.o  testutils.o -lpnetcdf -o test_get_varn
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_vard.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_vard.o  testutils.o -lpnetcdf -o test_vard
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_vard_multiple.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_vard_multiple.o  testutils.o -lpnetcdf -o test_vard_multiple
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_vard_rec.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_vard_rec.o  testutils.o -lpnetcdf -o test_vard_rec
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/test_varm.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  test_varm.o  testutils.o -lpnetcdf -o test_varm
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/tst_def_var_fill.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  tst_def_var_fill.o  testutils.o -lpnetcdf -o tst_def_var_fill
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/tst_dimsizes.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  tst_dimsizes.o  testutils.o -lpnetcdf -o tst_dimsizes
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/tst_free_comm.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  tst_free_comm.o  testutils.o -lpnetcdf -o tst_free_comm
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/tst_max_var_dims.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  tst_max_var_dims.o  testutils.o -lpnetcdf -o tst_max_var_dims
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/tst_version.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  tst_version.o  testutils.o -lpnetcdf -o tst_version
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/varn_contig.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  varn_contig.o  testutils.o -lpnetcdf -o varn_contig
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/varn_int.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  varn_int.o  testutils.o -lpnetcdf -o varn_int
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../testcases/vectors.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  vectors.o  testutils.o -lpnetcdf -o vectors
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../common/testutilsf.F90 -o testutilsf.o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o f90tst_parallel.exe90 ../F90/f90tst_parallel.f90 testutilsf.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o f90tst_parallel2.exe90 ../F90/f90tst_parallel2.f90 testutilsf.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o f90tst_parallel3.exe90 ../F90/f90tst_parallel3.f90 testutilsf.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o f90tst_parallel4.exe90 ../F90/f90tst_parallel4.f90 testutilsf.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/block_cyclic.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  block_cyclic.o  testutils.o -lpnetcdf -o block_cyclic
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/bput_varn_int64.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  bput_varn_int64.o  testutils.o -lpnetcdf -o bput_varn_int64
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/bput_varn_uint.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  bput_varn_uint.o  testutils.o -lpnetcdf -o bput_varn_uint
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/collective_write.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  collective_write.o  testutils.o -lpnetcdf -o collective_write
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/column_wise.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  column_wise.o  testutils.o -lpnetcdf -o column_wise
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/create_open.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  create_open.o  testutils.o -lpnetcdf -o create_open
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/fill_mode.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  fill_mode.o  testutils.o -lpnetcdf -o fill_mode
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/flexible_api.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  flexible_api.o  testutils.o -lpnetcdf -o flexible_api
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/get_info.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  get_info.o  testutils.o -lpnetcdf -o get_info
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/ghost_cell.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  ghost_cell.o  testutils.o -lpnetcdf -o ghost_cell
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/global_attributes.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  global_attributes.o  testutils.o -lpnetcdf -o global_attributes
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/hints.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  hints.o  testutils.o -lpnetcdf -o hints
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/mput.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  mput.o  testutils.o -lpnetcdf -o mput
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/nonblocking_write.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  nonblocking_write.o  testutils.o -lpnetcdf -o nonblocking_write
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/nonblocking_write_in_def.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  nonblocking_write_in_def.o  testutils.o -lpnetcdf -o nonblocking_write_in_def
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/put_vara.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  put_vara.o  testutils.o -lpnetcdf -o put_vara
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/get_vara.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  get_vara.o  testutils.o -lpnetcdf -o get_vara
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/put_varn_float.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  put_varn_float.o  testutils.o -lpnetcdf -o put_varn_float
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/put_varn_int.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  put_varn_int.o  testutils.o -lpnetcdf -o put_varn_int
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/time_var.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  time_var.o  testutils.o -lpnetcdf -o time_var
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/transpose2D.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  transpose2D.o  testutils.o -lpnetcdf -o transpose2D
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/transpose.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  transpose.o  testutils.o -lpnetcdf -o transpose
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/vard_int.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  vard_int.o  testutils.o -lpnetcdf -o vard_int
  cc -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -c ../../examples/C/vard_mvars.c
  cc -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib  vard_mvars.o  testutils.o -lpnetcdf -o vard_mvars
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/block_cyclic.f -o block_cyclic.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o block_cyclic.exe77 block_cyclic.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/bput_varn_int8.f -o bput_varn_int8.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o bput_varn_int8.exe77 bput_varn_int8.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/column_wise.f -o column_wise.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o column_wise.exe77 column_wise.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/fill_mode.f -o fill_mode.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o fill_mode.exe77 fill_mode.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/flexible_api.f -o flexible_api.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o flexible_api.exe77 flexible_api.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/get_info.f -o get_info.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o get_info.exe77 get_info.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/hints.f -o hints.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o hints.exe77 hints.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/i_varn_real.f -o i_varn_real.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o i_varn_real.exe77 i_varn_real.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/nonblocking_write.f -o nonblocking_write.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o nonblocking_write.exe77 nonblocking_write.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/put_vara.f -o put_vara.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_vara.exe77 put_vara.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/put_varn_int.f -o put_varn_int.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_varn_int.exe77 put_varn_int.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/put_varn_real.f -o put_varn_real.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_varn_real.exe77 put_varn_real.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/time_var.f -o time_var.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o time_var.exe77 time_var.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/transpose.f -o transpose.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o transpose.exe77 transpose.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F77/vard_int.f -o vard_int.77o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o vard_int.exe77 vard_int.77o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/block_cyclic.f90 -o block_cyclic.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o block_cyclic.exe90 block_cyclic.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/column_wise.f90 -o column_wise.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o column_wise.exe90 column_wise.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/fill_mode.f90 -o fill_mode.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o fill_mode.exe90 fill_mode.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/flexible_api.f90 -o flexible_api.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o flexible_api.exe90 flexible_api.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/get_info.f90 -o get_info.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o get_info.exe90 get_info.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/hints.f90 -o hints.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o hints.exe90 hints.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/nonblocking_write.f90 -o nonblocking_write.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o nonblocking_write.exe90 nonblocking_write.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/put_var.f90 -o put_var.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_var.exe90 put_var.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/put_varn_int.f90 -o put_varn_int.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_varn_int.exe90 put_varn_int.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/put_varn_real.f90 -o put_varn_real.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o put_varn_real.exe90 put_varn_real.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/transpose.f90 -o transpose.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o transpose.exe90 transpose.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -c ../../examples/F90/vard_int.f90 -o vard_int.90o
  ftn -I/opt/cray/pe/parallel-netcdf/1.12.3.1/include -I../common -w -fallow-argument-mismatch -o vard_int.exe90 vard_int.90o utils.o -L/opt/cray/pe/parallel-netcdf/1.12.3.1/lib -lpnetcdf
  ```
  </details>

* Two shell script files, `batch.sh` and `interactive.sh` will be created after
  make command. They can be used to run in batch and interactive modes,
  respectively. Please edit them to adjust batch environment settings, such as
  the account number, job queue, number of nodes, and number of MPI tasks per
  node.

  <details>
  <summary>Example run output shown on screen (click to expand)</summary>

  ```console
  login02::test_installed(1:00pm) #1097 sbatch batch.sh
  Submitted batch job 5291610
  login02::test_installed(1:01pm) #1098 sqs
  JOBID            ST USER      NAME          NODES TIME_LIMIT       TIME  SUBMIT_TIME          QOS             START_TIME           FEATURES       NODELIST(REASON
  5291610          R  wkliao    pnetcdf_test  1           5:00       0:44  2023-02-05T13:00:17  debug           2023-02-05T13:00:33  cpu            nid005504      

  login02::test_installed(1:04pm) #1099 cat out.5291610
  *** TESTING C   pres_temp_4D_wr for writing classic file           ------ pass
  *** TESTING C   pres_temp_4D_rd for reading classic file           ------ pass
  *** TESTING C   header_consistency for header consistency          ------ pass
  *** TESTING C   flexible_bput for flexible bput_varm               ------ pass
  *** TESTING C   interleaved for writing interleaved fileviews      ------ pass
  *** TESTING C   req_all for NC_REQ_ALL                             ------ pass
  *** TESTING C   i_varn_indef for iput/iget varn in define mode     ------ pass
  *** TESTING C   large_num_reqs for large number of iput/iget       ------ pass
  *** TESTING C   test_bput for bput API                             ------ pass
  *** TESTING C   i_varn_int64 for iput/iget varn                    ------ pass
  *** TESTING C   mcoll_perf for mput/iput APIs                      ------ pass
  *** TESTING C   wait_after_indep for ncmpi_end_indep_data          ------ pass
  *** TESTING C   add_var for checking offsets of new variables      ------ pass
  *** TESTING C   alignment_test for alignment                       ------ pass
  *** TESTING C   buftype_free for free buftype in flexible API      ------ pass
  *** TESTING C   check_striping for striping info                   ------ pass
  *** TESTING C   check_type for checking for type conflict          ------ pass
  *** TESTING C   collective_error for collective abort              ------ pass
  *** TESTING C   flexible for flexible put and get                  ------ pass
  *** TESTING C   flexible2 for flexible APIs                        ------ pass
  *** TESTING C   flexible_varm for flexible varm APIs               ------ pass
  *** TESTING C   inq_num_vars for no. record/fixed variables        ------ pass
  *** TESTING C   inq_recsize for inquiring record size              ------ pass
  *** TESTING C   ivarn for ncmpi_iput_varn_<type>()                 ------ pass
  *** TESTING C   large_var_cdf5 for large var in CDF-5              ------ pass
  *** TESTING C   last_large_var for last large var in CDF-1/2       ------ pass
  *** TESTING C   mix_collectives for get/put varm                   ------ pass
  *** TESTING C   modes for file create/open modes                   ------ pass
  *** TESTING C   ncmpi_vars_null_stride for NULL stride             ------ pass
  *** TESTING C   noclobber for NC_NOCLOBBER and NC_EEXIST           ------ pass
  *** TESTING C   nonblocking for using ncmpi_iput_vara_int()        ------ pass
  *** TESTING C   one_record for only one record variable            ------ pass
  *** TESTING C   record for write records in reversed order         ------ pass
  *** TESTING C   redef1 for entering re-define mode                 ------ pass
  *** TESTING C   scalar for get/put scalar variables                ------ pass
  *** TESTING C   test_erange for checking for NC_ERANGE             ------ pass
  *** TESTING C   test_fillvalue for _FillValue for NC_GLOBAL        ------ pass
  *** TESTING C   test_get_varn for get_varn                         ------ pass
  *** TESTING C   test_vard for vard put and get                     ------ pass
  *** TESTING C   test_vard_multiple for vard to 2 variables         ------ pass
  *** TESTING C   test_vard_rec for vard put on record var           ------ pass
  *** TESTING C   test_varm for get/put varm                         ------ pass
  *** TESTING C   tst_def_var_fill for def_var_fill                  ------ pass
  *** TESTING C   tst_dimsizes for defining max dimension sizes      ------ pass
  *** TESTING C   tst_free_comm for freeing MPI communicator         ------ pass
  *** TESTING C   tst_max_var_dims for checking NC_MAX_VAR_DIMS      ------ skip
  *** TESTING C   tst_version for PnetCDF library version            ------ pass
  *** TESTING C   varn_contig for put_varn with contig fileview      ------ pass
  *** TESTING C   varn_int for ncmpi_put_varn_int_all()              ------ pass
  *** TESTING C   vectors for put_vara/get_vara                      ------ pass
  *** TESTING F90 f90tst_parallel.exe90                              ------ pass
  *** TESTING F90 f90tst_parallel2.exe90 for strided access          ------ pass
  *** TESTING F90 f90tst_parallel3.exe90                             ------ pass
  *** TESTING F90 f90tst_parallel4.exe90                             ------ pass
  *** TESTING C   examples/C/block_cyclic                            ------ pass
  *** TESTING C   examples/C/bput_varn_int64                         ------ pass
  *** TESTING C   examples/C/bput_varn_uint                          ------ pass
  *** TESTING C   examples/C/collective_write                        ------ pass
  *** TESTING C   examples/C/column_wise                             ------ pass
  *** TESTING C   examples/C/create_open                             ------ pass
  *** TESTING C   examples/C/fill_mode                               ------ pass
  *** TESTING C   examples/C/flexible_api                            ------ pass
  *** TESTING C   examples/C/get_info                                ------ pass
  *** TESTING C   examples/C/ghost_cell                              ------ pass
  *** TESTING C   examples/C/global_attributes                       ------ pass
  *** TESTING C   examples/C/hints                                   ------ pass
  *** TESTING C   examples/C/mput                                    ------ pass
  *** TESTING C   examples/C/nonblocking_write                       ------ pass
  *** TESTING C   examples/C/nonblocking_write_in_def                ------ pass
  *** TESTING C   examples/C/put_vara                                ------ pass
  *** TESTING C   examples/C/get_vara                                ------ pass
  *** TESTING C   examples/C/put_varn_float                          ------ pass
  *** TESTING C   examples/C/put_varn_int                            ------ pass
  *** TESTING C   examples/C/time_var                                ------ pass
  *** TESTING C   examples/C/transpose2D                             ------ pass
  *** TESTING C   examples/C/transpose                               ------ pass
  *** TESTING C   examples/C/vard_int                                ------ pass
  *** TESTING C   examples/C/vard_mvars                              ------ pass
  *** TESTING F77 examples/F77/block_cyclic.exe77                    ------ pass
  *** TESTING F77 examples/F77/bput_varn_int8.exe77                  ------ pass
  *** TESTING F77 examples/F77/column_wise.exe77                     ------ pass
  *** TESTING F77 examples/F77/fill_mode.exe77                       ------ pass
  *** TESTING F77 examples/F77/flexible_api.exe77                    ------ pass
  *** TESTING F77 examples/F77/get_info.exe77                        ------ pass
  *** TESTING F77 examples/F77/hints.exe77                           ------ pass
  *** TESTING F77 examples/F77/i_varn_real.exe77                     ------ pass
  *** TESTING F77 examples/F77/nonblocking_write.exe77               ------ pass
  *** TESTING F77 examples/F77/put_vara.exe77                        ------ pass
  *** TESTING F77 examples/F77/put_varn_int.exe77                    ------ pass
  *** TESTING F77 examples/F77/put_varn_real.exe77                   ------ pass
  *** TESTING F77 examples/F77/time_var.exe77                        ------ pass
  *** TESTING F77 examples/F77/transpose.exe77                       ------ pass
  *** TESTING F77 examples/F77/vard_int.exe77                        ------ pass
  *** TESTING F90 examples/F90/block_cyclic.exe90                    ------ pass
  *** TESTING F90 examples/F90/column_wise.exe90                     ------ pass
  *** TESTING F90 examples/F90/fill_mode.exe90                       ------ pass
  *** TESTING F90 examples/F90/flexible_api.exe90                    ------ pass
  *** TESTING F90 examples/F90/get_info.exe90                        ------ pass
  *** TESTING F90 examples/F90/hints.exe90                           ------ pass
  *** TESTING F90 examples/F90/nonblocking_write.exe90               ------ pass
  *** TESTING F90 examples/F90/put_var.exe90                         ------ pass
  *** TESTING F90 examples/F90/put_varn_int.exe90                    ------ pass
  *** TESTING F90 examples/F90/put_varn_real.exe90                   ------ pass
  *** TESTING F90 examples/F90/transpose.exe90                       ------ pass
  *** TESTING F90 examples/F90/vard_int.exe90                        ------ pass

  Total number of tested programs: 105

  ```
  </details>
