!
!   Copyright (C) 2013, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example shows how to use nfmpi_put_vara_int_all().
! It first defines a netCDF variable as a 2D array of size
! global_nx * global_ny of 4-byte integer type where global_nx == 5 and
! global_ny == (4 * number of MPI processes). Each process writes a subarray
! of size nx * ny. Note the description above follows the Fortran order.
!
! Example MPI run and output from ncmpidump:
!
!    % mpiexec -l -n 4 ./put_vara
!
!    % ncmpidump testfile.nc
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
!      1, 2, 3, 4, 5,
!      6, 7, 8, 9, 10,
!      11, 12, 13, 14, 15,
!      16, 17, 18, 19, 20,
!      1001, 1002, 1003, 1004, 1005,
!      1006, 1007, 1008, 1009, 1010,
!      1011, 1012, 1013, 1014, 1015,
!      1016, 1017, 1018, 1019, 1020,
!      2001, 2002, 2003, 2004, 2005,
!      2006, 2007, 2008, 2009, 2010,
!      2011, 2012, 2013, 2014, 2015,
!      2016, 2017, 2018, 2019, 2020,
!      3001, 3002, 3003, 3004, 3005,
!      3006, 3007, 3008, 3009, 3010,
!      3011, 3012, 3013, 3014, 3015,
!      3016, 3017, 3018, 3019, 3020 ;
!    }
!
      subroutine check(err, message)
          use mpi
          use pnetcdf
          implicit none
          integer err
          character(len=*) message

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) trim(message), trim(nfmpi_strerror(err))
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end subroutine check

      program main
          use mpi
          use pnetcdf
          implicit none

          character(LEN=128) filename
          integer i, err, nprocs, rank, cmode, ncid, varid, dimid(2)
          integer buf(20)
          integer(kind=MPI_OFFSET_KIND) nx, ny, global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) starts(2), counts(2)

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! set parameters
          filename  = 'testfile.nc'
          nx        = 5
          ny        = 4
          global_nx = nx
          global_ny = ny * nprocs

          do i = 1, 20
             buf(i) = rank * 1000 + i
          end do

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     &                        MPI_INFO_NULL, ncid)
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

          call MPI_Finalize(err)
      end program main

