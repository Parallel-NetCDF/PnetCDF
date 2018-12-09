!
!   Copyright (C) 2018, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!

!
! This program tests Fortran 90 for adding attribute _FillValue.
!
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
              msg = '*** TESTING F90 test_fill.f90 '
              call pass_fail(1, msg)
              STOP 2
          end if
      end subroutine check

      integer function tst_fmt(filename, mode)
          use mpi
          use pnetcdf
          implicit none

          character(LEN=256) filename
          integer i, err, ierr, rank
          integer :: ncid, mode, cmode, dimid(1), varid
          integer(kind=MPI_OFFSET_KIND) :: start(1)
          integer(kind=MPI_OFFSET_KIND) :: count(1)
          integer(kind=MPI_OFFSET_KIND), parameter :: len = 3
          integer, parameter :: k = selected_int_kind(18)
          integer(kind=k) :: buf(len)

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

          tst_fmt = 0

          do i=1,len
             buf(i) = rank + i
          end do

          ! create netcdf file
          cmode = IOR(mode, NF90_CLOBBER)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')
          tst_fmt = tst_fmt + err

          err = nf90mpi_def_dim(ncid, "dim", len, dimid(1))
          call check(err, 'In nf90mpi_def_dim: ')
          tst_fmt = tst_fmt + err

          ! Make variable
          err = nf90mpi_def_var(ncid, "var", NF90_INT64, dimid, varid)
          call check(err, 'In nf90mpi_def_var: ')
          tst_fmt = tst_fmt + err

          ! new scalar attribute
          err = nf90mpi_put_att(ncid, varid, 'att', int(2, kind=k))
          call check(err, 'In nf90mpi_put_att: ')

          ! add scalar attribute "_FillValue"
          err = nf90mpi_put_att(ncid, varid, '_FillValue', int(1, kind=k))
          call check(err, 'In nf90mpi_put_att: ')
          tst_fmt = tst_fmt + err

          err = nf90mpi_enddef(ncid)
          call check(err, 'In nf90mpi_enddef: ')
          tst_fmt = tst_fmt + err

          ! Write buf
          start(1) = 1
          count(1) = len
          err = nf90mpi_put_var_all(ncid, varid, buf, start, count)
          call check(err, 'In nf90mpi_put_var_all: ')
          tst_fmt = tst_fmt + err

          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')
          tst_fmt = tst_fmt + err

      end function tst_fmt

      program test
          use mpi
          use pnetcdf
          implicit none
          character(LEN=256) filename, cmd, msg
          integer err, ierr, nerrs, rank, get_args, tst_fmt

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              filename = 'testfile.nc'
              err = get_args(cmd, filename)
          endif
          call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          if (err .EQ. 0) goto 999

          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, err)

          nerrs = 0

          if (PNETCDF_DRIVER_NETCDF4 .EQ. 1) then
              err = tst_fmt(filename, NF90_NETCDF4)
              nerrs = nerrs + err
          endif

          err = tst_fmt(filename, NF90_64BIT_DATA)
          nerrs = nerrs + err

          msg = '*** TESTING F90 '//trim(cmd)//' for _FillValue '
          if (rank .eq. 0) call pass_fail(nerrs, msg)

 999      call MPI_Finalize(err)

      end program test
