!
!  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

      subroutine check(err, message)
          use mpi
          use pnetcdf
          implicit none
          integer err
          character(len=*) message
          character(len=256) msg

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF90_NOERR) then
              write(6,*) trim(message), trim(nf90mpi_strerror(err))
              msg = '*** TESTING F90 test_attr_int64.f90 '
              call pass_fail(1, msg)
              STOP 2
          end if
      end subroutine check

      program main
          use mpi
          use pnetcdf
          implicit none

          character(LEN=256) filename, cmd, msg
          integer rank, err, ierr, ncid, cmode, get_args, xtype, varid
          integer(kind=MPI_OFFSET_KIND) :: buf
          integer,parameter :: INT2_KIND = selected_int_kind(4)
          integer fillmode

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              filename = 'testfile.nc'
              err = get_args(cmd, filename)
          endif
          call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          if (err .EQ. 0) goto 999

          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

          cmode = IOR(NF90_64BIT_DATA, NF90_CLOBBER)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
                               MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          buf = 5
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'attr_ll', buf)
          call check(err, 'In nf90mpi_put_att: ')

          err = nf90mpi_inquire_attribute(ncid, NF90_GLOBAL, 'attr_ll', xtype)
          call check(err, 'In nf90mpi_inquire_attribute: ')

          if (xtype .NE. NF90_INT64) then
              msg = '*** TESTING F90 test_attr_int64.f90 '
              call pass_fail(1, msg)
              STOP 2
          endif

          err = nf90mpi_def_var(ncid, "var_int64", NF90_INT64, varid)
          call check(err, 'In nf90mpi_def_var var_int64: ')

          buf = -9
          err = nf90mpi_put_att(ncid, varid, '_FillValue', buf)
          call check(err, 'In nf90mpi_put_att var_int64: ')

          err = nf90mpi_def_var(ncid, "var_short", NF90_SHORT, varid)
          call check(err, 'In nf90mpi_def_var var_short: ')

          err = nf90mpi_put_att(ncid, varid, '_FillValue', INT(-999,KIND=INT2_KIND))
          call check(err, 'In nf90mpi_put_att var_short: ')

          err = nf90mpi_def_var(ncid, "var_int", NF90_INT, varid)
          call check(err, 'In nf90mpi_def_var var_int: ')

          err = nf90mpi_put_att(ncid, varid, '_FillValue', -999.9)
          if (err .NE. NF90_EBADTYPE) then
10            FORMAT(A,I3)
              write(msg,10) '*** test_attr_int64.f90 expects NF90_EBADTYPE but got ', err
              call pass_fail(1, msg)
              STOP 2
          endif

          err = nf90mpi_set_fill(ncid, nf90_fill, fillmode)
          call check(err, 'In nf90mpi_set_fill: ')

          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')

          msg = '*** TESTING F90 '//trim(cmd)//' for scalar attr of INT64 '
          if (rank .eq. 0) call pass_fail(0, msg)

 999      call MPI_Finalize(err)

      end program main
