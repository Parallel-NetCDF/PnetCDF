!
!  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example shows how to use a single call of nf90mpi_put_varn_all()
! to write a sequence of one-element requests with arbitrary array indices.
!
! The compile and run commands are given below, together with an nfmpidump of
! the output file.
!
!    % mpif90 -g -o put_varn_real put_varn_real.f90 -lpnetcdf
!    % mpiexec -n 4 put_varn_real testfile.nc
!    % ncmpidump testfile.nc
!    netcdf testfile {
!    // file format: CDF-5 (big variables)
!    dimensions:
!             Y = 4 ;
!             X = 10 ;
!    variables:
!             int var(Y, X) ;
!    data:
!
!     var =
!       3, 3, 3, 1, 1, 0, 0, 2, 1, 1,
!       0, 2, 2, 2, 3, 1, 1, 2, 2, 2,
!       1, 1, 2, 3, 3, 3, 0, 0, 1, 1,
!       0, 0, 0, 2, 1, 1, 1, 3, 3, 3 ;
!    }
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

          integer NDIMS
          PARAMETER(NDIMS=2)

          character(LEN=128) filename
          integer i, rank, nprocs, err, num_reqs, argc, iargc
          integer ncid, cmode, varid, dimid(2), y, x
          real buffer(13)
          integer(kind=MPI_OFFSET_KIND) NY, NX
          integer(kind=MPI_OFFSET_KIND) starts(NDIMS, 13)
          integer(kind=MPI_OFFSET_KIND) counts(NDIMS, 13)

          NY = 4
          NX = 10

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          argc = IARGC()
          if (argc .NE. 1) then
              print*,'Usage: block_cyclic filename'
              STOP
          endif
          call getarg(1, filename)

          if (nprocs .NE. 4 .AND. rank .EQ. 0) &
              print*,"Warning: this program is intended to run on ", &
                     "4 processes"

          ! create file, truncate it if exists
          cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
                             MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          ! create a global array of size NY * NX */
          err = nf90mpi_def_dim(ncid, "Y", NY, dimid(2))
          call check(err, 'In nf90mpi_def_dim Y: ')
          err = nf90mpi_def_dim(ncid, "X", NX, dimid(1))
          call check(err, 'In nf90mpi_def_dim X: ')
          err = nf90mpi_def_var(ncid, "var", NF90_FLOAT, dimid, varid)
          call check(err, 'In nf90mpi_def_var var: ')
          err = nf90mpi_enddef(ncid)
          call check(err, 'In nf90mpi_enddef: ')

          ! pick arbitrary numbers of requests for 4 processes
          num_reqs = 1
          if (rank .EQ.  0) then
              num_reqs = 8
          elseif (rank .EQ. 1) then
              num_reqs = 13
          elseif (rank .EQ. 2) then
              num_reqs = 9
          elseif (rank .EQ. 3) then
              num_reqs = 10
          endif

          ! assign arbitrary starts
          y=2
          x=1
          if (rank .EQ. 0) then
              starts(y, 1) = 1
              starts(x, 1) = 6
              starts(y, 2) = 2
              starts(x, 2) = 1
              starts(y, 3) = 3
              starts(x, 3) = 7
              starts(y, 4) = 4
              starts(x, 4) = 1
              starts(y, 5) = 1
              starts(x, 5) = 7
              starts(y, 6) = 3
              starts(x, 6) = 8
              starts(y, 7) = 4
              starts(x, 7) = 2
              starts(y, 8) = 4
              starts(x, 8) = 3
              ! rank 0 is writing the following locations: ("-" means skip)
              !   -  -  -  -  -  0  0  -  -  - 
              !   0  -  -  -  -  -  -  -  -  - 
              !   -  -  -  -  -  -  0  0  -  - 
              !   0  0  0  -  -  -  -  -  -  - 
          elseif (rank .EQ. 1) then
              starts(y,  1) = 1
              starts(x,  1) = 4
              starts(y,  2) = 1
              starts(x,  2) = 9
              starts(y,  3) = 2
              starts(x,  3) = 6
              starts(y,  4) = 3
              starts(x,  4) = 1
              starts(y,  5) = 3
              starts(x,  5) = 9
              starts(y,  6) = 4
              starts(x,  6) = 5
              starts(y,  7) = 1
              starts(x,  7) = 5
              starts(y,  8) = 1
              starts(x,  8) = 10
              starts(y,  9) = 2
              starts(x,  9) = 7
              starts(y, 10) = 3
              starts(x, 10) = 2
              starts(y, 11) = 3
              starts(x, 11) = 10
              starts(y, 12) = 4
              starts(x, 12) = 6
              starts(y, 13) = 4
              starts(x, 13) = 7
              ! rank 1 is writing the following locations: ("-" means skip)
              !   -  -  -  1  1  -  -  -  1  1 
              !   -  -  -  -  -  1  1  -  -  - 
              !   1  1  -  -  -  -  -  -  1  1 
              !   -  -  -  -  1  1  1  -  -  - 
          elseif (rank .EQ. 2) then
              starts(y, 1) = 1
              starts(x, 1) = 8
              starts(y, 2) = 2
              starts(x, 2) = 2
              starts(y, 3) = 2
              starts(x, 3) = 8
              starts(y, 4) = 3
              starts(x, 4) = 3
              starts(y, 5) = 4
              starts(x, 5) = 4
              starts(y, 6) = 2
              starts(x, 6) = 3
              starts(y, 7) = 2
              starts(x, 7) = 9
              starts(y, 8) = 2
              starts(x, 8) = 4
              starts(y, 9) = 2
              starts(x, 9) = 10
              ! rank 2 is writing the following locations: ("-" means skip)
              !   -  -  -  -  -  -  -  2  -  - 
              !   -  2  2  2  -  -  -  2  2  2 
              !   -  -  2  -  -  -  -  -  -  - 
              !   -  -  -  2  -  -  -  -  -  - 
          elseif (rank .EQ. 3) then
              starts(y,  1) = 1
              starts(x,  1) = 1
              starts(y,  2) = 2
              starts(x,  2) = 5
              starts(y,  3) = 3
              starts(x,  3) = 4
              starts(y,  4) = 4
              starts(x,  4) = 8
              starts(y,  5) = 1
              starts(x,  5) = 2
              starts(y,  6) = 3
              starts(x,  6) = 5
              starts(y,  7) = 4
              starts(x,  7) = 9
              starts(y,  8) = 1
              starts(x,  8) = 3
              starts(y,  9) = 3
              starts(x,  9) = 6
              starts(y, 10) = 4
              starts(x, 10) = 10
              ! rank 3 is writing the following locations: ("-" means skip)
              !   3  3  3  -  -  -  -  -  -  - 
              !   -  -  -  -  3  -  -  -  -  - 
              !   -  -  -  3  3  3  -  -  -  - 
              !   -  -  -  -  -  -  -  3  3  3 
          endif
          counts = 1

          ! allocate I/O buffer and initialize its contents
          buffer = rank

          ! set the buffer pointers to different offsets to the I/O buffer
          err = nf90mpi_put_varn_all(ncid, varid, buffer, num_reqs, &
                                     starts, counts)
          call check(err, 'In nf90mpi_put_varn_real_all: ')

          err = nf90mpi_close(ncid);
          call check(err, 'In nf90mpi_close: ')

          call MPI_Finalize(err)

      end program

