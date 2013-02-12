!
!   Copyright (C) 2013, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example shows how to use PnetCDF nonblocking flexible API,
! nfmpi_iput_vara() to write a 2D 4-byte integer array in parallel.
! It first defines a netCDF variable of size global_nx * global_ny where
!    global_nx == 5 and
!    global_ny == (4 * number of MPI processes).
! The data partitioning pattern is a column-wise partitioning across all
! proceses. Each process writes a subarray of size nx * ny.
! The local buffer has a ghost cell of length 3 surrounding the 2D array
!    integer buf(nx+2*ghost_len, ny+2*ghost_len)
! Note the description above follows the Fortran array index order.
!
! Example commands for MPI run and outputs from running ncmpidump on the
! NC file produced by this example program:
!
!    % mpiexec -l -n 4 ./flex_f testfile.nc
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
          integer argc, IARGC, err, nprocs, rank, i, j, ghost_len
          integer cmode, ncid, varid, dimid(2)
          integer(kind=MPI_OFFSET_KIND) nx, ny, global_nx, global_ny
          integer(kind=MPI_OFFSET_KIND) starts(2), counts(2), nReqs
          PARAMETER(nx=5, ny=4, ghost_len=3)
          integer buf(nx+2*ghost_len, ny+2*ghost_len)
          integer subarray, req(1), status(1)
          integer array_of_sizes(2), array_of_subsizes(2)
          integer array_of_starts(2)
          integer nfmpi_put_vara_all
          integer nfmpi_iput_vara

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          argc = IARGC()
          if (argc .GE. 1) then
              call getarg(1, filename)
          else
              filename  = 'testfile.nc'
          endif

          ! set parameters
          global_nx = nx
          global_ny = ny * nprocs

          ! first initialize the entire buffer to -1
          buf = -1;
          ! assign values for non-ghost cells
          do j=ghost_len+1, ny+ghost_len+1
             do i=ghost_len+1, nx+ghost_len+1
                 buf(i, j) = rank
             enddo
          enddo

          ! define an MPI datatype using MPI_Type_create_subarray()
          array_of_sizes(1)    = nx + 2*ghost_len
          array_of_sizes(2)    = ny + 2*ghost_len
          array_of_subsizes(1) = nx
          array_of_subsizes(2) = ny
          array_of_starts(1)   = ghost_len+1
          array_of_starts(2)   = ghost_len+1
          call MPI_Type_create_subarray(2, array_of_sizes,
     &         array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN,
     &         MPI_INT, subarray, err)
          call MPI_Type_commit(subarray, err)

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
          nReqs = 1

          err = nfmpi_iput_vara(ncid, varid, starts, counts, buf,
     &          nReqs, subarray, req)
          call check(err, 'In nfmpi_iput_vara: ')

          err = nfmpi_wait_all(ncid, 1, req, status)
          call check(err, 'In nfmpi_wait_all: ')
          call check(status(1), 'In nfmpi_wait_all status error: ')

          call MPI_Type_free(subarray, err)

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

          call MPI_Finalize(err)
      end program main

