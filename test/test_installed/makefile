.SUFFIXES: .o .exe77 .exe90 .77o .90o
.PRECIOUS: %.o %.77o %.90o

CPPFLAGS = -I$(PNETCDF_DIR)/include -I../common
CFLAGS = $(CPPFLAGS)
FCFLAGS = -w -fallow-argument-mismatch
LDFLAGS = -L$(PNETCDF_DIR)/lib
LDLIBS = testutils.o -lpnetcdf

TEST_PROGS =

C_src = ../C/pres_temp_4D_wr.c \
        ../C/pres_temp_4D_rd.c
TEST_PROGS += $(C_src:../C/%.c=%)
%.o: ../C/%.c
	$(CC) $(CFLAGS) -c $<

header_src = ../header/header_consistency.c
TEST_PROGS += $(header_src:../header/%.c=%)
%.o: ../header/%.c
	$(CC) $(CFLAGS) -c $<

nonblocking_src = ../nonblocking/flexible_bput.c \
                  ../nonblocking/interleaved.c \
                  ../nonblocking/req_all.c \
                  ../nonblocking/i_varn_indef.c \
                  ../nonblocking/large_num_reqs.c \
                  ../nonblocking/test_bput.c \
                  ../nonblocking/i_varn_int64.c \
                  ../nonblocking/mcoll_perf.c \
                  ../nonblocking/wait_after_indep.c
TEST_PROGS += $(nonblocking_src:../nonblocking/%.c=%)
%.o: ../nonblocking/%.c
	$(CC) $(CFLAGS) -c $<

testcases_src = ../testcases/add_var.c \
                ../testcases/alignment_test.c \
                ../testcases/buftype_free.c \
                ../testcases/check_striping.c \
                ../testcases/check_type.c \
                ../testcases/collective_error.c \
                ../testcases/flexible.c \
                ../testcases/flexible2.c \
                ../testcases/flexible_varm.c \
                ../testcases/inq_num_vars.c \
                ../testcases/inq_recsize.c \
                ../testcases/ivarn.c \
                ../testcases/large_var_cdf5.c \
                ../testcases/last_large_var.c \
                ../testcases/mix_collectives.c \
                ../testcases/modes.c \
                ../testcases/ncmpi_vars_null_stride.c \
                ../testcases/noclobber.c \
                ../testcases/nonblocking.c \
                ../testcases/one_record.c \
                ../testcases/record.c \
                ../testcases/redef1.c \
                ../testcases/scalar.c \
                ../testcases/test_erange.c \
                ../testcases/test_fillvalue.c \
                ../testcases/test_get_varn.c \
                ../testcases/test_vard.c \
                ../testcases/test_vard_multiple.c \
                ../testcases/test_vard_rec.c \
                ../testcases/test_varm.c \
                ../testcases/tst_def_var_fill.c \
                ../testcases/tst_dimsizes.c \
                ../testcases/tst_free_comm.c \
                ../testcases/tst_max_var_dims.c \
                ../testcases/tst_version.c \
                ../testcases/varn_contig.c \
                ../testcases/varn_int.c \
                ../testcases/vectors.c
TEST_PROGS += $(testcases_src:../testcases/%.c=%)
%.o: ../testcases/%.c
	$(CC) $(CFLAGS) -c $<

EXAMPLE_PROGS =
examples_C_src = ../../examples/C/block_cyclic.c \
                 ../../examples/C/bput_varn_int64.c \
                 ../../examples/C/bput_varn_uint.c \
                 ../../examples/C/collective_write.c \
                 ../../examples/C/column_wise.c \
                 ../../examples/C/create_open.c \
                 ../../examples/C/fill_mode.c \
                 ../../examples/C/flexible_api.c \
                 ../../examples/C/get_info.c \
                 ../../examples/C/ghost_cell.c \
                 ../../examples/C/global_attributes.c \
                 ../../examples/C/hints.c \
                 ../../examples/C/mput.c \
                 ../../examples/C/nonblocking_write.c \
                 ../../examples/C/nonblocking_write_in_def.c \
                 ../../examples/C/put_vara.c \
                 ../../examples/C/get_vara.c \
                 ../../examples/C/put_varn_float.c \
                 ../../examples/C/put_varn_int.c \
                 ../../examples/C/time_var.c \
                 ../../examples/C/transpose2D.c \
                 ../../examples/C/transpose.c \
                 ../../examples/C/vard_int.c \
                 ../../examples/C/vard_mvars.c
EXAMPLE_PROGS += $(examples_C_src:../../examples/C/%.c=%)
%.o: ../../examples/C/%.c
	$(CC) $(CFLAGS) -c $<

F90_src = ../F90/f90tst_parallel.f90 \
          ../F90/f90tst_parallel2.f90 \
          ../F90/f90tst_parallel3.f90 \
          ../F90/f90tst_parallel4.f90
TEST_PROGS += $(F90_src:../F90/%.f90=%.exe90)

examples_F77_src = ../../examples/F77/block_cyclic.f \
                   ../../examples/F77/bput_varn_int8.f \
                   ../../examples/F77/column_wise.f \
                   ../../examples/F77/fill_mode.f \
                   ../../examples/F77/flexible_api.f \
                   ../../examples/F77/get_info.f \
                   ../../examples/F77/hints.f \
                   ../../examples/F77/i_varn_real.f \
                   ../../examples/F77/nonblocking_write.f \
                   ../../examples/F77/put_vara.f \
                   ../../examples/F77/put_varn_int.f \
                   ../../examples/F77/put_varn_real.f \
                   ../../examples/F77/time_var.f \
                   ../../examples/F77/transpose.f \
                   ../../examples/F77/vard_int.f
EXAMPLE_PROGS += $(examples_F77_src:../../examples/F77/%.f=%.exe77)

%.77o: ../../examples/F77/%.f
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $< -o $@

examples_F90_src = ../../examples/F90/block_cyclic.f90 \
                   ../../examples/F90/column_wise.f90 \
                   ../../examples/F90/fill_mode.f90 \
                   ../../examples/F90/flexible_api.f90 \
                   ../../examples/F90/get_info.f90 \
                   ../../examples/F90/hints.f90 \
                   ../../examples/F90/nonblocking_write.f90 \
                   ../../examples/F90/put_var.f90 \
                   ../../examples/F90/put_varn_int.f90 \
                   ../../examples/F90/put_varn_real.f90 \
                   ../../examples/F90/transpose.f90 \
                   ../../examples/F90/vard_int.f90
EXAMPLE_PROGS += $(examples_F90_src:../../examples/F90/%.f90=%.exe90)

all: env_check testutils.o utils.o $(TEST_PROGS) $(EXAMPLE_PROGS) batch.sh interactive.sh

f90tst_parallel.exe90: ../F90/f90tst_parallel.f90 testutilsf.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< testutilsf.o $(LDFLAGS) -lpnetcdf
f90tst_parallel2.exe90: ../F90/f90tst_parallel2.f90 testutilsf.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< testutilsf.o $(LDFLAGS) -lpnetcdf
f90tst_parallel3.exe90: ../F90/f90tst_parallel3.f90 testutilsf.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< testutilsf.o $(LDFLAGS) -lpnetcdf
f90tst_parallel4.exe90: ../F90/f90tst_parallel4.f90 testutilsf.o
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< testutilsf.o $(LDFLAGS) -lpnetcdf

%.90o: ../../examples/F90/%.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $< -o $@

.77o.exe77:
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< utils.o $(LDFLAGS) -lpnetcdf

.90o.exe90:
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $@ $< utils.o $(LDFLAGS) -lpnetcdf

env_check:
	@if test "x$(CC)" = x ; then \
	    echo "Please set environment CC to the MPI C compiler" ; \
	    exit 1; \
	fi
	@if test "x$(FC)" = x ; then \
	    echo "Please set environment FC to the MPI Fortran compiler" ; \
	    exit 1; \
	fi
	@if test "x$(RUN_CMD)" = x ; then \
	    echo "Please set environment RUN_CMD to the MPI run command (e.g. mpiexec or srun)" ; \
	    exit 1; \
	fi
	@if test "x$(OUT_DIR)" = x ; then \
	    echo "Please set environment OUT_DIR to the folder for storing output files" ; \
	    exit 1; \
	fi
	@if test "x$(PNETCDF_DIR)" = x ; then \
	    echo "Please set environment PNETCDF_DIR to the PnetCDF install path" ; \
	    exit 1; \
	fi

testutils.o: ../common/testutils.c
	$(CC) $(CFLAGS) -I../common -c $< -o $@

testutilsf.o: ../common/testutilsf.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $< -o $@

utils.o: ../../examples/F90/utils.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $< -o $@

interactive.sh: interactive.sh.in makefile
	@sed -e "s|__TEST_PROGS__|$(TEST_PROGS)|" $< > $@
	@sed -i -e "s|__EXAMPLE_PROGS__|$(EXAMPLE_PROGS)|" $@
	@sed -i -e "s|__RUN_CMD__|$(RUN_CMD)|" $@
	@sed -i -e "s|__OUT_DIR__|$(OUT_DIR)|" $@

batch.sh: batch.sh.in makefile
	@sed -e "s|__TEST_PROGS__|$(TEST_PROGS)|" $< > $@
	@sed -i -e "s|__EXAMPLE_PROGS__|$(EXAMPLE_PROGS)|" $@
	@sed -i -e "s|__RUN_CMD__|$(RUN_CMD)|" $@
	@sed -i -e "s|__OUT_DIR__|$(OUT_DIR)|" $@

clean:
	rm -f *.o *.77o *.90o $(TEST_PROGS) $(EXAMPLE_PROGS) batch.sh interactive.sh

