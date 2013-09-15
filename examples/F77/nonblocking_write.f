!
!  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!    This example is the Fortran 77 version of nonblocking_write.c
!    It creates a netcdf file in CD-5 format and writes a number of
!    3D integer non-record vaiables.
!    Usage: (for example)
!    To compile:
!        mpif77 -O2 nonblocking_write.f -o nonblocking_write -lpnetcdf
!    To run (for example):
!        mpiexec -n 32 ./nonblocking_write 10 /pvfs2/wkliao/testfile.nc
!    The size of each local array is len x len x len. Each non-record
!    variable is of size len*len*len * nprocs * sizeof(int)
!    All variables are partitioned among all processes in a 3D
!    block-block-block fashion.
!

      subroutine check(err, message)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          integer err
          character(len=*) message

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) trim(message), trim(nfmpi_strerror(err))
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end subroutine check

      program main
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          integer NDIMS, NUM_VARS, BUFSIZE
          PARAMETER(NDIMS=3, NUM_VARS=10, BUFSIZE=512)

          character(LEN=128) filename, cmd, str
          integer i, j, cmode, argc, iargc, err
          integer rank, nprocs, ncid, len
          integer buf(BUFSIZE, NUM_VARS)
          integer psizes(NDIMS), dimids(NDIMS), varids(NUM_VARS)
          integer req(NUM_VARS), st(NUM_VARS)
          integer(kind=MPI_OFFSET_KIND) gsizes(NDIMS)
          integer(kind=MPI_OFFSET_KIND) starts(NDIMS), counts(NDIMS)

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 2) then
              print*,'Usage: ',trim(cmd),' len [filename]'
              goto 999
          endif
          call getarg(1, str)
          read (str,'(I10)') len
          filename = "testfile.nc"
          if (argc .EQ. 1) call getarg(2, filename)

          do i=1,NDIMS
             psizes(i) = 0
          enddo

          ! create a block-block-block data partitioning
          call MPI_Dims_create(nprocs, NDIMS, psizes, err);
          starts(1) = mod(rank, psizes(1))
          starts(2) = mod((rank / psizes(2)), psizes(2))
          starts(3) = mod((rank / (psizes(1) * psizes(2))), psizes(3))

          do i=1,NDIMS
              gsizes(i) = len * psizes(i)
              starts(i) = starts(i) * len + 1
              counts(i) = len
          enddo

          buf = rank

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define dimensions
          do i=1, NDIMS
             write(str,'(I2)') i
             err = nfmpi_def_dim(ncid, "x"//trim(ADJUSTL(str)),
     +                           gsizes(i), dimids(i))
             call check(err, 'In nfmpi_def_dim x'//trim(str))
          enddo

          !define variables
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nfmpi_def_var(ncid, "var"//trim(ADJUSTL(str)),
     +                           NF_INT, NDIMS, dimids, varids(i))
             call check(err, 'In nfmpi_def_var var'//trim(str))
          enddo

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! write one variable at a time
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nfmpi_iput_vara_int(ncid, varids(i), starts,
     +                                 counts, buf(:,i), req(i))
             call check(err, 'In nfmpi_iput_vara_int '//trim(str))
          enddo

          ! wait for the nonblocking I/O to complete
          err = nfmpi_wait_all(ncid, NUM_VARS, req, st)
          call check(err, 'In nfmpi_wait_all')

          ! check the status of each nonblocking request
          do i=1, NUM_VARS
             write(str,'(I2)') i
             call check(st(i), 'In nfmpi_wait_all req '//trim(str))
          enddo

          ! close the file
          err = nfmpi_close(ncid);
          call check(err, 'In nfmpi_close')

 999      call MPI_Finalize(err)

      end program main

