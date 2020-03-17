!  Copyright (C) 2003, Northwestern University and Argonne National
!  Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

C Copyright 1996-2006, UCAR/Unidata
C See COPYRIGHT file for copying and redistribution conditions.

C Test driver for netCDF-3 interface.  This program performs tests
C against
C the netCDF-3 specification for all user-level functions in an
C implementation of the netCDF library.


        subroutine usage()
        implicit        none
        include "pnetcdf.inc"
#include "tests.inc"

        integer MY_LEN_TRIM
        call error('Usage: '//TRIM(progname)//
     +             ' [-hrv -n <MAX_NMPT>]')
        call error('       nf_test [-h] Print help' )
        call error('   [-1] test CDF-1 format' )
        call error('   [-2] test CDF-2 format' )
        call error('   [-5] test CDF-5 format' )
#ifdef ENABLE_NETCDF4
        call error('   [-4] test NetCDF-4 classic-model format' )
#endif
        call error('   [-r] Just do read-only tests' )
        call error(
     +  '   [-d directory] directory for storing input/output files' )
        call error('   [-v] Verbose mode' )
        call error(
     +  '   [-n <max>] max. number of messages per test (Default: 20)')
        end


        subroutine report_test
        implicit        none
        character*1024   msg
#include "tests.inc"

        integer MY_LEN_TRIM
        if (cdf_format .EQ. 4) then
          write(msg,"(A)") '*** TESTING F77 '//
     +                     TRIM(PROGNAME)//
     +                     ' for NetCDF-4'
        else
          write(msg,"(A,I1)") '*** TESTING F77 '//
     +                        TRIM(PROGNAME)//
     +                        ' for CDF-', cdf_format
        endif
        if (nfailsTotal .ne. 0) then
          write(*,*) TRIM(PROGNAME)//
     +               ' expects to see 0 failure ... '//
     +               'Total number of failures: ', nfailsTotal
        endif
        call pass_fail(nfailsTotal, msg)
        end

        subroutine test(name, func)
        implicit        none
        character*(*)   name
        character*25    name_str
        integer         name_len
        external        func
        include "pnetcdf.inc"
#include "tests.inc"

        integer MY_LEN_TRIM
        name_len = MY_LEN_TRIM(name)
        name_str(1:name_len) = name(:)
        name_str(name_len+1:25) = ' '
        if (verbose) write(*, 1) name_str
1       format('*** Testing ', a, ' ... ')

        nfails = 0
        call func()
        nfailsTotal = nfailsTotal + nfails
        if ( nfails .ne. 0) then
            print *, ' '
            print *, '  ### ', nfails, ' FAILURES TESTING ', name,
     +               '! Stop ... ###'
            call report_test
            stop 2
        end if
        end


#if _CRAYIEEE
! which machines need this?
        subroutine getarg(iarg, carg)
        implicit        none
        integer iarg
        character*(*)   carg
        integer ilen
        integer ierror
        call PXFGETARG(iarg, carg, ilen, ierror)
        end
#endif

        program nf_test
#if defined(VISUAL_CPLUSPLUS)
!       DIGITAL Visual Fortran needs DFLIB for getarg
        USE DFLIB
!       DIGITAL Visual Fortran needs DFPORT for iargc
        USE DFPORT
        implicit        none
#elif defined(NAGFOR)
        USE F90_UNIX_ENV, only : iargc, getarg
        implicit none
#else
        implicit        none
        integer         iargc
#endif
#if defined(__crayx1)
        integer         ipxfargc
#endif
        include "pnetcdf.inc"
#include "tests.inc"

        integer         argc
        character*80    arg
        integer         iarg
        integer         iopt
        character*1     opt
        integer         lastopt
        logical         skiparg
        integer         err
        integer MY_LEN_TRIM

        external        test_nfmpi_strerror
        external        test_nfmpi_open
        external        test_nfmpi_close
        external        test_nfmpi_inq
        external        test_nfmpi_inq_dimid
        external        test_nfmpi_inq_dim
        external        test_nfmpi_inq_dimlen
        external        test_nfmpi_inq_dimname
        external        test_nfmpi_inq_varid
        external        test_nfmpi_inq_var
        external        test_nfmpi_inq_natts
        external        test_nfmpi_inq_ndims
        external        test_nfmpi_inq_nvars
        external        test_nfmpi_inq_unlimdim
        external        test_nfmpi_inq_vardimid
        external        test_nfmpi_inq_varname
        external        test_nfmpi_inq_varnatts
        external        test_nfmpi_inq_varndims
        external        test_nfmpi_inq_vartype
        external        test_nfmpi_get_var1_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_var1_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_var1_int2
#endif
        external        test_nfmpi_get_var1_int
        external        test_nfmpi_get_var1_real
        external        test_nfmpi_get_var1_double
        external        test_nfmpi_get_var1_int8
        external        test_nfmpi_get_var_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_var_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_var_int2
#endif
        external        test_nfmpi_get_var_int
        external        test_nfmpi_get_var_real
        external        test_nfmpi_get_var_double
        external        test_nfmpi_get_var_int8
        external        test_nfmpi_get_vara_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_vara_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_vara_int2
#endif
        external        test_nfmpi_get_vara_int
        external        test_nfmpi_get_vara_real
        external        test_nfmpi_get_vara_double
        external        test_nfmpi_get_vara_int8
        external        test_nfmpi_get_vars_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_vars_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_vars_int2
#endif
        external        test_nfmpi_get_vars_int
        external        test_nfmpi_get_vars_real
        external        test_nfmpi_get_vars_double
        external        test_nfmpi_get_vars_int8

        external        test_nfmpi_get_varm_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_varm_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_varm_int2
#endif
        external        test_nfmpi_get_varm_int
        external        test_nfmpi_get_varm_real
        external        test_nfmpi_get_varm_double
        external        test_nfmpi_get_varm_int8

        external        test_nfmpi_iget_var1_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iget_var1_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iget_var1_int2
#endif
        external        test_nfmpi_iget_var1_int
        external        test_nfmpi_iget_var1_real
        external        test_nfmpi_iget_var1_double
        external        test_nfmpi_iget_var1_int8
        external        test_nfmpi_iget_var_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iget_var_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iget_var_int2
#endif
        external        test_nfmpi_iget_var_int
        external        test_nfmpi_iget_var_real
        external        test_nfmpi_iget_var_double
        external        test_nfmpi_iget_var_int8
        external        test_nfmpi_iget_vara_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iget_vara_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iget_vara_int2
#endif
        external        test_nfmpi_iget_vara_int
        external        test_nfmpi_iget_vara_real
        external        test_nfmpi_iget_vara_double
        external        test_nfmpi_iget_vara_int8
        external        test_nfmpi_iget_vars_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iget_vars_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iget_vars_int2
#endif
        external        test_nfmpi_iget_vars_int
        external        test_nfmpi_iget_vars_real
        external        test_nfmpi_iget_vars_double
        external        test_nfmpi_iget_vars_int8

        external        test_nfmpi_iget_varm_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iget_varm_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iget_varm_int2
#endif
        external        test_nfmpi_iget_varm_int
        external        test_nfmpi_iget_varm_real
        external        test_nfmpi_iget_varm_double
        external        test_nfmpi_iget_varm_int8

        external        test_nfmpi_get_att_text
#if defined(NF_INT1_T)
        external        test_nfmpi_get_att_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_get_att_int2
#endif
        external        test_nfmpi_get_att_int
        external        test_nfmpi_get_att_real
        external        test_nfmpi_get_att_double
        external        test_nfmpi_get_att_int8
        external        test_nfmpi_inq_att
        external        test_nfmpi_inq_attname
        external        test_nfmpi_inq_attid
        external        test_nfmpi_inq_attlen
        external        test_nfmpi_inq_atttype
        external        test_nfmpi_create
        external        test_nfmpi_redef
        external        test_nfmpi_enddef
        external        test_nfmpi_sync
        external        test_nfmpi_flush
        external        test_nfmpi_abort
        external        test_nfmpi_def_dim
        external        test_nfmpi_rename_dim
        external        test_nfmpi_def_var
        external        test_nfmpi_put_var1_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_var1_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_var1_int2
#endif
        external        test_nfmpi_put_var1_int
        external        test_nfmpi_put_var1_real
        external        test_nfmpi_put_var1_double
        external        test_nfmpi_put_var1_int8
        external        test_nfmpi_put_var_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_var_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_var_int2
#endif
        external        test_nfmpi_put_var_int
        external        test_nfmpi_put_var_real
        external        test_nfmpi_put_var_double
        external        test_nfmpi_put_var_int8
        external        test_nfmpi_put_vara_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_vara_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_vara_int2
#endif
        external        test_nfmpi_put_vara_int
        external        test_nfmpi_put_vara_real
        external        test_nfmpi_put_vara_double
        external        test_nfmpi_put_vara_int8
        external        test_nfmpi_put_vars_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_vars_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_vars_int2
#endif
        external        test_nfmpi_put_vars_int
        external        test_nfmpi_put_vars_real
        external        test_nfmpi_put_vars_double
        external        test_nfmpi_put_vars_int8

        external        test_nfmpi_put_varm_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_varm_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_varm_int2
#endif
        external        test_nfmpi_put_varm_int
        external        test_nfmpi_put_varm_real
        external        test_nfmpi_put_varm_double
        external        test_nfmpi_put_varm_int8

        external        test_nfmpi_iput_var1_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iput_var1_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iput_var1_int2
#endif
        external        test_nfmpi_iput_var1_int
        external        test_nfmpi_iput_var1_real
        external        test_nfmpi_iput_var1_double
        external        test_nfmpi_iput_var1_int8
        external        test_nfmpi_iput_var_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iput_var_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iput_var_int2
#endif
        external        test_nfmpi_iput_var_int
        external        test_nfmpi_iput_var_real
        external        test_nfmpi_iput_var_double
        external        test_nfmpi_iput_var_int8
        external        test_nfmpi_iput_vara_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iput_vara_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iput_vara_int2
#endif
        external        test_nfmpi_iput_vara_int
        external        test_nfmpi_iput_vara_real
        external        test_nfmpi_iput_vara_double
        external        test_nfmpi_iput_vara_int8
        external        test_nfmpi_iput_vars_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iput_vars_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iput_vars_int2
#endif
        external        test_nfmpi_iput_vars_int
        external        test_nfmpi_iput_vars_real
        external        test_nfmpi_iput_vars_double
        external        test_nfmpi_iput_vars_int8

        external        test_nfmpi_iput_varm_text
#if defined(NF_INT1_T)
        external        test_nfmpi_iput_varm_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_iput_varm_int2
#endif
        external        test_nfmpi_iput_varm_int
        external        test_nfmpi_iput_varm_real
        external        test_nfmpi_iput_varm_double
        external        test_nfmpi_iput_varm_int8

        external        test_nfmpi_rename_var
        external        test_nfmpi_put_att_text
#if defined(NF_INT1_T)
        external        test_nfmpi_put_att_int1
#endif
#if defined(NF_INT2_T)
        external        test_nfmpi_put_att_int2
#endif
        external        test_nfmpi_put_att_int
        external        test_nfmpi_put_att_real
        external        test_nfmpi_put_att_double
        external        test_nfmpi_put_att_int8
        external        test_nfmpi_copy_att
        external        test_nfmpi_rename_att
        external        test_nfmpi_del_att
        external        test_nfmpi_set_fill
        external        test_nfmpi_set_default_format
        external        nc_ignorefpe

        call MPI_INIT(err)
        comm = MPI_COMM_WORLD

        call nc_ignorefpe(1)

        testfile = 'test.nc'
        scratch = 'scratch.nc'

        nfailsTotal = 0
        call getarg(0, progname)
        readonly = .false.      !/* assume may write in test dir as default */
        verbose = .false.
        max_nmpt = 20
        skiparg = .false.
        cdf_format = 1
        extra_flags = 0

#if defined(__crayx1)
        argc = ipxfargc()
#else
        argc = iargc()
#endif
        call getarg(0, PROGNAME)

        do 1, iarg = 1, argc
            if (skiparg) then
                skiparg = .false.
            else
                call getarg(iarg, arg)
                if (arg(1:1) .eq. '-') then
                    lastopt = index(arg, ' ') - 1
                    do 2, iopt = 2, lastopt
                        opt = arg(iopt:iopt)
                        if (opt .eq. 'r') then
                            readonly = .true.
                        else if (opt .eq. 'v') then
                            verbose = .true.
                        else if (opt .eq. 'n') then
                            call getarg(iarg+1, arg)
                            ! NOTE: The UNICOS 8 fort77(1) compiler does
                            ! not support list-directed I/O from an internal
                            ! file -- so we use a format specification.
                            read (arg, '(i6)') max_nmpt
                            skiparg = .true.
                            go to 1
                        else if (opt .eq. '1') then
                            cdf_format = 1
                        else if (opt .eq. '2') then
                            cdf_format = 2
                            extra_flags = NF_64BIT_OFFSET
                        else if (opt .eq. '4') then
                            cdf_format = 4
                            extra_flags=IOR(NF_NETCDF4,NF_CLASSIC_MODEL)
                        else if (opt .eq. '5') then
                            cdf_format = 5
                            extra_flags = NF_64BIT_DATA
                        else if (opt .eq. 'd') then
                            call getarg(iarg+1, arg)
                            testfile = TRIM(arg)//
     +                                 "/test.nc"
                            scratch = TRIM(arg) //
     +                                 "/scratch.nc"
                            skiparg = .true.
                            go to 1
                        else
                            call usage
                            call ud_exit(1)
                        end if
    2           continue
                else
                    call usage
                    call ud_exit(1)
                end if
            end if
1       continue

#ifndef ENABLE_NETCDF4
        if (cdf_format .EQ. 4) then
            call error(
     +      "Error: NetCDF-4 support is not enabled at configure time")
            call MPI_Finalize(err)
            stop 1
        endif
#endif

        numGatts = 6
        numVars  = 136
        numTypes = 6
        if (cdf_format .EQ. 5) then
            numGatts = NGATTS
            numVars  = NVARS
            numTypes = NTYPES
        endif

C       PVFS2 driver has a problem of ADIOI_Set_lock when data sieving is
C       enabled
        call MPI_Info_create(info, err)
C       call MPI_Info_set(info, "romio_pvfs2_posix_write", "enable",err)
C       disable MPI-IO data sieving
        call MPI_Info_set(info, "romio_ds_write", "disable", err)
        call MPI_Info_set(info, "romio_lustre_ds_in_coll","disable",err)

C       /* Initialize global variables defining test file */
        call init_gvars

C       create file test.nc for testing read operations */
        call write_file(testfile)
        if (nfailsTotal .GT. 0) then
            call MPI_Info_free(info, err)
            call ud_exit(1)
        end if

C       /* delete any existing scratch netCDF file */
        if ( .not. readonly )
     +      err = nfmpi_delete(scratch, MPI_INFO_NULL)

C       /* Test read-only functions, using pregenerated test-file */
        call test('nfmpi_strerror', test_nfmpi_strerror)
        call test('nfmpi_open', test_nfmpi_open)
        call test('nfmpi_close', test_nfmpi_close)
        call test('nfmpi_inq', test_nfmpi_inq)
        call test('nfmpi_inq_dimid', test_nfmpi_inq_dimid)
        call test('nfmpi_inq_dim', test_nfmpi_inq_dim)
        call test('nfmpi_inq_dimlen', test_nfmpi_inq_dimlen)
        call test('nfmpi_inq_dimname', test_nfmpi_inq_dimname)
        call test('nfmpi_inq_varid', test_nfmpi_inq_varid)
        call test('nfmpi_inq_var', test_nfmpi_inq_var)
        call test('nfmpi_inq_natts', test_nfmpi_inq_natts)
        call test('nfmpi_inq_ndims', test_nfmpi_inq_ndims)
        call test('nfmpi_inq_nvars', test_nfmpi_inq_nvars)
        call test('nfmpi_inq_unlimdim', test_nfmpi_inq_unlimdim)
        call test('nfmpi_inq_vardimid', test_nfmpi_inq_vardimid)
        call test('nfmpi_inq_varname', test_nfmpi_inq_varname)
        call test('nfmpi_inq_varnatts', test_nfmpi_inq_varnatts)
        call test('nfmpi_inq_varndims', test_nfmpi_inq_varndims)
        call test('nfmpi_inq_vartype', test_nfmpi_inq_vartype)

        call test('nfmpi_get_var1_text', test_nfmpi_get_var1_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_var1_int1', test_nfmpi_get_var1_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_var1_int2', test_nfmpi_get_var1_int2)
#endif
        call test('nfmpi_get_var1_int', test_nfmpi_get_var1_int)
        call test('nfmpi_get_var1_real', test_nfmpi_get_var1_real)
        call test('nfmpi_get_var1_double', test_nfmpi_get_var1_double)
        call test('nfmpi_get_var1_int8', test_nfmpi_get_var1_int8)

        call test('nfmpi_get_var_text', test_nfmpi_get_var_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_var_int1', test_nfmpi_get_var_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_var_int2', test_nfmpi_get_var_int2)
#endif
        call test('nfmpi_get_var_int', test_nfmpi_get_var_int)
        call test('nfmpi_get_var_real', test_nfmpi_get_var_real)
        call test('nfmpi_get_var_double', test_nfmpi_get_var_double)
        call test('nfmpi_get_var_int8', test_nfmpi_get_var_int8)

        call test('nfmpi_get_vara_text', test_nfmpi_get_vara_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_vara_int1', test_nfmpi_get_vara_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_vara_int2', test_nfmpi_get_vara_int2)
#endif
        call test('nfmpi_get_vara_int', test_nfmpi_get_vara_int)
        call test('nfmpi_get_vara_real', test_nfmpi_get_vara_real)
        call test('nfmpi_get_vara_double', test_nfmpi_get_vara_double)
        call test('nfmpi_get_vara_int8', test_nfmpi_get_vara_int8)

        call test('nfmpi_get_vars_text', test_nfmpi_get_vars_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_vars_int1', test_nfmpi_get_vars_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_vars_int2', test_nfmpi_get_vars_int2)
#endif
        call test('nfmpi_get_vars_int', test_nfmpi_get_vars_int)
        call test('nfmpi_get_vars_real', test_nfmpi_get_vars_real)
        call test('nfmpi_get_vars_double', test_nfmpi_get_vars_double)
        call test('nfmpi_get_vars_int8', test_nfmpi_get_vars_int8)

        call test('nfmpi_get_varm_text', test_nfmpi_get_varm_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_varm_int1', test_nfmpi_get_varm_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_varm_int2', test_nfmpi_get_varm_int2)
#endif
        call test('nfmpi_get_varm_int', test_nfmpi_get_varm_int)
        call test('nfmpi_get_varm_real', test_nfmpi_get_varm_real)
        call test('nfmpi_get_varm_double', test_nfmpi_get_varm_double)
        call test('nfmpi_get_varm_int8', test_nfmpi_get_varm_int8)

        if (cdf_format .ne. 4) then
        call test('nfmpi_iget_var1_text', test_nfmpi_iget_var1_text)
#if defined(NF_INT1_T)
        call test('nfmpi_iget_var1_int1', test_nfmpi_iget_var1_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_iget_var1_int2', test_nfmpi_iget_var1_int2)
#endif
        call test('nfmpi_iget_var1_int', test_nfmpi_iget_var1_int)
        call test('nfmpi_iget_var1_real', test_nfmpi_iget_var1_real)
        call test('nfmpi_iget_var1_double', test_nfmpi_iget_var1_double)
        call test('nfmpi_iget_var1_int8', test_nfmpi_iget_var1_int8)

        call test('nfmpi_iget_var_text', test_nfmpi_iget_var_text)
#if defined(NF_INT1_T)
        call test('nfmpi_iget_var_int1', test_nfmpi_iget_var_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_iget_var_int2', test_nfmpi_iget_var_int2)
#endif
        call test('nfmpi_iget_var_int', test_nfmpi_iget_var_int)
        call test('nfmpi_iget_var_real', test_nfmpi_iget_var_real)
        call test('nfmpi_iget_var_double', test_nfmpi_iget_var_double)
        call test('nfmpi_iget_var_int8', test_nfmpi_iget_var_int8)

        call test('nfmpi_iget_vara_text', test_nfmpi_iget_vara_text)
#if defined(NF_INT1_T)
        call test('nfmpi_iget_vara_int1', test_nfmpi_iget_vara_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_iget_vara_int2', test_nfmpi_iget_vara_int2)
#endif
        call test('nfmpi_iget_vara_int', test_nfmpi_iget_vara_int)
        call test('nfmpi_iget_vara_real', test_nfmpi_iget_vara_real)
        call test('nfmpi_iget_vara_double', test_nfmpi_iget_vara_double)
        call test('nfmpi_iget_vara_int8', test_nfmpi_iget_vara_int8)

        call test('nfmpi_iget_vars_text', test_nfmpi_iget_vars_text)
#if defined(NF_INT1_T)
        call test('nfmpi_iget_vars_int1', test_nfmpi_iget_vars_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_iget_vars_int2', test_nfmpi_iget_vars_int2)
#endif
        call test('nfmpi_iget_vars_int', test_nfmpi_iget_vars_int)
        call test('nfmpi_iget_vars_real', test_nfmpi_iget_vars_real)
        call test('nfmpi_iget_vars_double', test_nfmpi_iget_vars_double)
        call test('nfmpi_iget_vars_int8', test_nfmpi_iget_vars_int8)

        call test('nfmpi_iget_varm_text', test_nfmpi_iget_varm_text)
#if defined(NF_INT1_T)
        call test('nfmpi_iget_varm_int1', test_nfmpi_iget_varm_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_iget_varm_int2', test_nfmpi_iget_varm_int2)
#endif
        call test('nfmpi_iget_varm_int', test_nfmpi_iget_varm_int)
        call test('nfmpi_iget_varm_real', test_nfmpi_iget_varm_real)
        call test('nfmpi_iget_varm_double', test_nfmpi_iget_varm_double)
        call test('nfmpi_iget_varm_int8', test_nfmpi_iget_varm_int8)
        end if

        call test('nfmpi_get_att_text', test_nfmpi_get_att_text)
#if defined(NF_INT1_T)
        call test('nfmpi_get_att_int1', test_nfmpi_get_att_int1)
#endif
#if defined(NF_INT2_T)
        call test('nfmpi_get_att_int2', test_nfmpi_get_att_int2)
#endif
        call test('nfmpi_get_att_int', test_nfmpi_get_att_int)
        call test('nfmpi_get_att_real', test_nfmpi_get_att_real)
        call test('nfmpi_get_att_double', test_nfmpi_get_att_double)
        call test('nfmpi_get_att_int8', test_nfmpi_get_att_int8)
        call test('nfmpi_inq_att', test_nfmpi_inq_att)
        call test('nfmpi_inq_attname', test_nfmpi_inq_attname)
        call test('nfmpi_inq_attid', test_nfmpi_inq_attid)
        call test('nfmpi_inq_attlen', test_nfmpi_inq_attlen)
        call test('nfmpi_inq_atttype', test_nfmpi_inq_atttype)

C           /* Test write functions */
        if (.not. readonly) then
            call test('nfmpi_create', test_nfmpi_create)
            call test('nfmpi_redef', test_nfmpi_redef)
            call test('nfmpi_enddef', test_nfmpi_enddef)
            call test('nfmpi_sync', test_nfmpi_sync)
            if (cdf_format .ne. 4) then
                call test('nfmpi_flush', test_nfmpi_flush)
            endif
            call test('nfmpi_abort', test_nfmpi_abort)
            call test('nfmpi_def_dim', test_nfmpi_def_dim)
            call test('nfmpi_rename_dim', test_nfmpi_rename_dim)
            call test('nfmpi_def_var', test_nfmpi_def_var)
            call test('nfmpi_put_var1_text', test_nfmpi_put_var1_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_var1_int1', test_nfmpi_put_var1_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_var1_int2', test_nfmpi_put_var1_int2)
#endif
            call test('nfmpi_put_var1_int', test_nfmpi_put_var1_int)
            call test('nfmpi_put_var1_real', test_nfmpi_put_var1_real)
            call test('nfmpi_put_var1_double',
     +                 test_nfmpi_put_var1_double)
            call test('nfmpi_put_var1_int8', test_nfmpi_put_var1_int8)
            call test('nfmpi_put_var_text', test_nfmpi_put_var_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_var_int1', test_nfmpi_put_var_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_var_int2', test_nfmpi_put_var_int2)
#endif
            call test('nfmpi_put_var_int', test_nfmpi_put_var_int)
            call test('nfmpi_put_var_real', test_nfmpi_put_var_real)
            call test('nfmpi_put_var_double',
     +                 test_nfmpi_put_var_double)
            call test('nfmpi_put_var_int8', test_nfmpi_put_var_int8)
            call test('nfmpi_put_vara_text', test_nfmpi_put_vara_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_vara_int1', test_nfmpi_put_vara_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_vara_int2', test_nfmpi_put_vara_int2)
#endif
            call test('nfmpi_put_vara_int', test_nfmpi_put_vara_int)
            call test('nfmpi_put_vara_real', test_nfmpi_put_vara_real)
            call test('nfmpi_put_vara_double',
     +                 test_nfmpi_put_vara_double)
            call test('nfmpi_put_vara_int8', test_nfmpi_put_vara_int8)
            call test('nfmpi_put_vars_text', test_nfmpi_put_vars_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_vars_int1', test_nfmpi_put_vars_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_vars_int2', test_nfmpi_put_vars_int2)
#endif
            call test('nfmpi_put_vars_int', test_nfmpi_put_vars_int)
            call test('nfmpi_put_vars_real', test_nfmpi_put_vars_real)
            call test('nfmpi_put_vars_double',
     +                 test_nfmpi_put_vars_double)
            call test('nfmpi_put_vars_int8', test_nfmpi_put_vars_int8)

            call test('nfmpi_put_varm_text', test_nfmpi_put_varm_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_varm_int1', test_nfmpi_put_varm_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_varm_int2', test_nfmpi_put_varm_int2)
#endif
            call test('nfmpi_put_varm_int', test_nfmpi_put_varm_int)
            call test('nfmpi_put_varm_real', test_nfmpi_put_varm_real)
            call test('nfmpi_put_varm_double',
     +                 test_nfmpi_put_varm_double)
            call test('nfmpi_put_varm_int8', test_nfmpi_put_varm_int8)

        if (cdf_format .ne. 4) then
            call test('nfmpi_iput_var1_text', test_nfmpi_iput_var1_text)
#if defined(NF_INT1_T)
            call test('nfmpi_iput_var1_int1', test_nfmpi_iput_var1_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_iput_var1_int2', test_nfmpi_iput_var1_int2)
#endif
            call test('nfmpi_iput_var1_int', test_nfmpi_iput_var1_int)
            call test('nfmpi_iput_var1_real', test_nfmpi_iput_var1_real)
            call test('nfmpi_iput_var1_double',
     +                 test_nfmpi_iput_var1_double)
            call test('nfmpi_iput_var1_int8', test_nfmpi_iput_var1_int8)

            call test('nfmpi_iput_var_text', test_nfmpi_iput_var_text)
#if defined(NF_INT1_T)
            call test('nfmpi_iput_var_int1', test_nfmpi_iput_var_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_iput_var_int2', test_nfmpi_iput_var_int2)
#endif
            call test('nfmpi_iput_var_int', test_nfmpi_iput_var_int)
            call test('nfmpi_iput_var_real', test_nfmpi_iput_var_real)
            call test('nfmpi_iput_var_double',
     +                 test_nfmpi_iput_var_double)
            call test('nfmpi_iput_var_int8', test_nfmpi_iput_var_int8)

            call test('nfmpi_iput_vara_text', test_nfmpi_iput_vara_text)
#if defined(NF_INT1_T)
            call test('nfmpi_iput_vara_int1', test_nfmpi_iput_vara_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_iput_vara_int2', test_nfmpi_iput_vara_int2)
#endif
            call test('nfmpi_iput_vara_int', test_nfmpi_iput_vara_int)
            call test('nfmpi_iput_vara_real', test_nfmpi_iput_vara_real)
            call test('nfmpi_iput_vara_double',
     +                 test_nfmpi_iput_vara_double)
            call test('nfmpi_iput_vara_int8', test_nfmpi_iput_vara_int8)

            call test('nfmpi_iput_vars_text', test_nfmpi_iput_vars_text)
#if defined(NF_INT1_T)
            call test('nfmpi_iput_vars_int1', test_nfmpi_iput_vars_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_iput_vars_int2', test_nfmpi_iput_vars_int2)
#endif
            call test('nfmpi_iput_vars_int', test_nfmpi_iput_vars_int)
            call test('nfmpi_iput_vars_real', test_nfmpi_iput_vars_real)
            call test('nfmpi_iput_vars_double',
     +                 test_nfmpi_iput_vars_double)
            call test('nfmpi_iput_vars_int8', test_nfmpi_iput_vars_int8)

            call test('nfmpi_iput_varm_text', test_nfmpi_iput_varm_text)
#if defined(NF_INT1_T)
            call test('nfmpi_iput_varm_int1', test_nfmpi_iput_varm_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_iput_varm_int2', test_nfmpi_iput_varm_int2)
#endif
            call test('nfmpi_iput_varm_int', test_nfmpi_iput_varm_int)
            call test('nfmpi_iput_varm_real', test_nfmpi_iput_varm_real)
            call test('nfmpi_iput_varm_double',
     +                 test_nfmpi_iput_varm_double)
            call test('nfmpi_iput_varm_int8', test_nfmpi_iput_varm_int8)
        end if

            call test('nfmpi_rename_var', test_nfmpi_rename_var)
            call test('nfmpi_put_att_text', test_nfmpi_put_att_text)
#if defined(NF_INT1_T)
            call test('nfmpi_put_att_int1', test_nfmpi_put_att_int1)
#endif
#if defined(NF_INT2_T)
            call test('nfmpi_put_att_int2', test_nfmpi_put_att_int2)
#endif
            call test('nfmpi_put_att_int', test_nfmpi_put_att_int)
            call test('nfmpi_put_att_real', test_nfmpi_put_att_real)
            call test('nfmpi_put_att_double',
     +                 test_nfmpi_put_att_double)
            call test('nfmpi_put_att_int8', test_nfmpi_put_att_int8)
            call test('nfmpi_copy_att', test_nfmpi_copy_att)
            call test('nfmpi_rename_att', test_nfmpi_rename_att)
            call test('nfmpi_del_att', test_nfmpi_del_att)
            if (cdf_format .ne. 4) then
            call test('nfmpi_set_fill', test_nfmpi_set_fill)
            end if
            call test('nfmpi_set_default_format',
     +                test_nfmpi_set_default_format)
        end if

        call MPI_Info_free(info, err)

        call report_test

        ! if (nfailsTotal .eq. 0) call ud_exit(0)
        call ud_exit(0)
        end
