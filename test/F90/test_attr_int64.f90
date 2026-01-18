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
              call pass_fail(1, msg, 0)
              STOP 2
          end if
      end subroutine check

      program main
          use mpi
          use pnetcdf
          implicit none

          character(LEN=256) out_path, in_path, cmd, msg
          integer rank, err, ierr, ncid, cmode, get_args, xtype, varid
          integer(kind=MPI_OFFSET_KIND) :: buf
          integer,parameter :: INT2_KIND = selected_int_kind(4)
          integer fillmode
          logical keep_files
          double precision timing

          call MPI_Init(ierr)

          timing = MPI_Wtime()

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          ! take out_path from command-line argument if there is any
          cmd = ' '
          if (rank .EQ. 0) then
              out_path = 'testfile.nc'
              err = get_args(cmd, out_path, in_path, keep_files)
          endif
          call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          if (err .EQ. 0) goto 999

          call MPI_Bcast(out_path, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

          call MPI_Bcast(keep_files, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

          cmode = IOR(NF90_64BIT_DATA, NF90_CLOBBER)
          err = nf90mpi_create(MPI_COMM_WORLD, out_path, cmode, &
                               MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          buf = 5
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'attr_ll', buf)
          call check(err, 'In nf90mpi_put_att: ')

          err = nf90mpi_inquire_attribute(ncid, NF90_GLOBAL, 'attr_ll', xtype)
          call check(err, 'In nf90mpi_inquire_attribute: ')

          if (xtype .NE. NF90_INT64) then
              msg = '*** TESTING F90 test_attr_int64.f90 '
              call pass_fail(1, msg, 0)
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
              call pass_fail(1, msg, 0)
              STOP 2
          endif

          err = nf90mpi_set_fill(ncid, nf90_fill, fillmode)
          call check(err, 'In nf90mpi_set_fill: ')

          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')

 999      timing = MPI_Wtime() - timing
          call MPI_Allreduce(MPI_IN_PLACE, timing, 1, &
                             MPI_DOUBLE_PRECISION, MPI_MAX, &
                             MPI_COMM_WORLD, ierr)

          if (rank .eq. 0) then
              if (.NOT. keep_files) then
                  err = nf90mpi_delete(out_path, MPI_INFO_NULL)
              end if

              msg = '*** TESTING F90 '//trim(cmd)//' - scalar attr of INT64 '
              call pass_fail(0, msg, timing)
          end if

          call MPI_Finalize(err)

      end program main
