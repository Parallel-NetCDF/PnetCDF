!
!   Copyright (C) 2015, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

! This program tests if NF_ERANGE is properly returned with a coredump
! when an out-of-range value is used to write to a global attribute.
! When using NAG Fortran compiler, "Arithmetic exception" and coredump
! happens.
!
!    % mpif77 -O2 -o attrf attrf.f -lpnetcdf
!    % mpiexec -n 1 ./attrf /pvfs2/wkliao/testfile.nc
!

       INTEGER FUNCTION XTRIM(STRING)
           CHARACTER(LEN=*) STRING
           INTEGER I, N
           N = LEN(STRING)
           DO I = N, 1, -1
              IF (STRING(I:I) .NE. ' ') GOTO 10
           ENDDO
 10        XTRIM = I
       END ! FUNCTION XTRIM

      subroutine check(err, message, nerrs)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          integer err, nerrs, XTRIM
          character(LEN=*) message
          character(LEN=128) msg

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) message(1:XTRIM(message)), nfmpi_strerror(err)
              msg = '*** TESTING F77 attrf.f for attribute overflow '
              call pass_fail(1, msg)
              nerrs = nerrs + 1
          end if
      end ! subroutine check

      program main
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          real                    buf_flt, one_flt
          double precision        buf_dbl
          integer                 buf_int, XTRIM, cmode
          integer*2 buf_int2
          integer*8 buf_int8, one

          character(LEN=256) filename, cmd, msg
          integer ncid, err, ierr, nerrs, nprocs, rank, get_args
          integer*8 malloc_size, sum_size

          call MPI_Init(ierr)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

          one = 1
          one_flt = 1.0

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              filename = "testfile.nc"
              err = get_args(cmd, filename)
          endif
          call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          if (err .EQ. 0) goto 999

          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0,
     +                   MPI_COMM_WORLD, ierr)

          nerrs = 0

          cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ', nerrs)

          ! use a number > X_INT_MAX (which is 2147483647)
          buf_flt  =  2147483648.0
          buf_dbl  =  2147483648.0
          buf_int8 =  1073741824
          buf_int8 =  buf_int8 * 2

          err = nfmpi_put_att_real(ncid, NF_GLOBAL, "attr1", NF_INT,
     +                             one, buf_flt)
          if (err .NE. NF_ERANGE) then
              print*, "Error: expect NF_ERANGE but got ", err
              if (err .NE. NF_NOERR) print*, nfmpi_strerror(err)
              nerrs = nerrs + 1
              ! Note: even with an error, the attribute is still being created
          endif
          err = nfmpi_put_att_real(ncid, NF_GLOBAL, "attr1", NF_INT,
     +                             one, one_flt)
          call check(err, 'In nfmpi_put_att_real: ', nerrs)

          err = nfmpi_put_att_double(ncid, NF_GLOBAL, "attr2", NF_INT,
     +                               one, buf_dbl)
          if (err .NE. NF_ERANGE) then
              print*, "Error: expect NF_ERANGE but got ", err
              if (err .NE. NF_NOERR) print*, nfmpi_strerror(err)
              nerrs = nerrs + 1
              ! Note: even with an error, the attribute is still being created
              ! in this case, valgrind complains about uninitialised buffer at
              ! nfmpi_enddef when writing header to file.
          endif
          err = nfmpi_put_att_real(ncid, NF_GLOBAL, "attr2", NF_INT,
     +                             one, one_flt)
          call check(err, 'In nfmpi_put_att_real: ', nerrs)

          err = nfmpi_put_att_int8(ncid, NF_GLOBAL, "attr3", NF_INT,
     +                             one, buf_int8)
          if (err .NE. NF_ERANGE) then
              print*, "Error: expect NF_ERANGE but got ", err
              if (err .NE. NF_NOERR) print*, nfmpi_strerror(err)
              nerrs = nerrs + 1
              ! Note: even with an error, the attribute is still being created
              ! in this case, valgrind complains about uninitialised buffer at
              ! nfmpi_enddef when writing header to file.
          endif
          err = nfmpi_put_att_int8(ncid, NF_GLOBAL, "attr3", NF_INT,
     +                             one, one)
          call check(err, 'In nfmpi_put_att_int8: ', nerrs)

          buf_int = 2147483647
          err = nfmpi_put_att_int(ncid, NF_GLOBAL, "attr4", NF_INT,
     +                            one, buf_int)
          call check(err, 'In nfmpi_put_att_int: ', nerrs)

          ! Because of the NF_ERANGE error, the attributes may become
          ! inconsistent among processes, So NC_EMULTIDEFINE_ATTR_VAL
          ! or NF_EMULTIDEFINE may be returned from nfmpi_enddef (1.7.0
          ! and before), or nfmpi_put_att_xxx (after 1.7.0 when in safe
          ! mode).
          ! While in safe mode, the warning message of inconsistent metadata
          ! may appear on the screen. This is expected.
          ! In addition, when running under valgrind, NF_ERANGE can
          ! cause skipping requests and a valgrind warning message of
          ! "uninitialised byte", which is also expected.
          err = nfmpi_enddef(ncid)
          if (err .NE. NF_NOERR .AND. err .NE. NF_EMULTIDEFINE .AND.
     +        err .NE. NF_EMULTIDEFINE_ATTR_VAL)
     +        call check(err, 'In nfmpi_enddef: ', nerrs)

          err = nfmpi_get_att_int2(ncid, NF_GLOBAL, "attr4", buf_int2)
          if (err .NE. NF_ERANGE) then
              print*, "Error: expect NF_ERANGE but got ", err
              if (err .NE. NF_NOERR) print*, nfmpi_strerror(err)
              nerrs = nerrs + 1
          endif

          ! check if can overwrite an attribute in data mode
          err = nfmpi_put_att_int(ncid, NF_GLOBAL, "attr4", NF_INT,
     +                            one, buf_int)
          call check(err, 'In nfmpi_put_att_int: ', nerrs)

          buf_int2 = 32767
          err = nfmpi_put_att_int2(ncid, NF_GLOBAL, "attr4", NF_INT2,
     +                             one, buf_int2)
          call check(err, 'In nfmpi_put_att_int: xxxxxxx', nerrs)

          err = nfmpi_put_att_int8(ncid, NF_GLOBAL, "attr4", NF_INT64,
     +                             one, buf_int8)
          if (err .NE. NF_ENOTINDEFINE)
     +       print*,'Error: expect error code NF_ENOTINDEFINE but got ',
     +       err

          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ', nerrs)

          ! check if there is any PnetCDF internal malloc residue
 998      format(A,I13,A)
          err = nfmpi_inq_malloc_size(malloc_size)
          if (err .EQ. NF_NOERR) then
              call MPI_Reduce(malloc_size, sum_size, 1, MPI_INTEGER8,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)
              if (rank .EQ. 0 .AND. sum_size .GT. 0)
     +            print 998,
     +            'heap memory allocated by PnetCDF internally has ',
     +            sum_size, ' bytes yet to be freed'
          endif

          if (rank .eq. 0) then
              msg = '*** TESTING F77 '//cmd(1:XTRIM(cmd))//
     +              ' for attribute overflow '
              call pass_fail(nerrs, msg)
          endif

 999      call MPI_Finalize(ierr)

          if (nerrs .GT. 0) STOP 2

      end ! program main
