!
!   Copyright (C) 2013, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example shows how to use nfmpi_put_vara_int_all() to write a 2D
! 4-byte integer array in parallel. It first defines a netCDF variable of
! size global_nx * global_ny where
!    global_nx == 5 and
!    global_ny == (4 * number of MPI processes).
! The data partitioning pattern is a column-wise partitioning across all
! proceses. Each process writes a subarray of size nx * ny.
! Note the description above follows the Fortran array index order.
!
! Example commands for MPI run and outputs from running ncmpidump on the
! NC file produced by this example program:
!
!    % mpif77 -O2 -o put_vara put_vara.f -lpnetcdf
!    % mpiexec -n 4 ./put_vara /pvfs2/wkliao/testfile.nc
!
!    % ncmpidump /pvfs2/wkliao/testfile.nc
!    netcdf testfile {
!    // file format: CDF-5 (big variables)
!    dimensions:
!            x = 5 ;
!            y = 16 ;
!    variables:
!            int var(y, x) ;
!    data:
!
!     var =
!      0, 0, 0, 0, 0,
!      0, 0, 0, 0, 0,
!      0, 0, 0, 0, 0,
!      0, 0, 0, 0, 0,
!      1, 1, 1, 1, 1,
!      1, 1, 1, 1, 1,
!      1, 1, 1, 1, 1,
!      1, 1, 1, 1, 1,
!      2, 2, 2, 2, 2,
!      2, 2, 2, 2, 2,
!      2, 2, 2, 2, 2,
!      2, 2, 2, 2, 2,
!      3, 3, 3, 3, 3,
!      3, 3, 3, 3, 3,
!      3, 3, 3, 3, 3,
!      3, 3, 3, 3, 3 ;
!    }
!
      subroutine check(err, message)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'
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
          include 'mpif.h'
          include 'pnetcdf.inc'

          character(LEN=128) filename, cmd
          integer argc, IARGC, err, nprocs, rank
          integer cmode, ncid, varid, dimid(2)
          integer(kind=MPI_OFFSET_KIND) nx, ny, global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) starts(2), counts(2)
          PARAMETER(nx=5, ny=4)
          integer buf(nx,ny)

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 1) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd),' [filename]'
              goto 999
          endif
          filename = "testfile.nc"
          if (argc .EQ. 1) call getarg(1, filename)

          ! set parameters
          global_nx = nx
          global_ny = ny * nprocs

          buf = rank;

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define dimensions x and y
          err = nfmpi_def_dim(ncid, "x", global_nx, dimid(1))
          call check(err, 'In nfmpi_def_dim x: ')
          err = nfmpi_def_dim(ncid, "y", global_ny, dimid(2))
          call check(err, 'In nfmpi_def_dim y: ')

          ! define a 2D variable of integer type
          err = nfmpi_def_var(ncid, "var", NF_INT, 2, dimid, varid)
          call check(err, 'In nfmpi_def_var: ')

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! now we are in data mode

          ! Note that in Fortran, array indices start with 1
          starts(1) = 1
          starts(2) = ny * rank + 1
          counts(1) = nx
          counts(2) = ny

          err = nfmpi_put_vara_int_all(ncid, varid, starts, counts, buf)
          call check(err, 'In nfmpi_put_vara_int_all: ')

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

 999      call MPI_Finalize(err)
      end program main

