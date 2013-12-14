!
!   Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This program tests put_att on Little Endian
!
      subroutine check(err, message)
          use mpi
          use pnetcdf
          implicit none
          integer err
          character(len=*) message

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF90_NOERR) then
              write(6,*) trim(message), trim(nf90mpi_strerror(err))
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end subroutine check

      program main
          use mpi
          use pnetcdf
          implicit none

          character(LEN=128) filename, cmd
          integer argc, IARGC, err, rank
          integer cmode, ncid

          integer buf(3)
          PARAMETER(buf=(/1,2,3/))

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 1) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd),' [filename]'
              goto 999
          endif
          filename = "testfile.nc"
          if (argc .EQ. 1) call getarg(1, filename)

          ! create file, truncate it if exists
          cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
                               MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'attr1', buf)
          call check(err, 'In nf90mpi_put_att: ')

          err = nfmpi_put_att_int(ncid, NF_GLOBAL, 'attr2', NF_INT, 3_8, buf)
          call check(err, 'In nfmpi_put_att_int: ')

          ! close the file
          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')

 999      call MPI_Finalize(err)
      end program main

