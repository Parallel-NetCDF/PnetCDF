!
!   Copyright (C) 2013, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example shows how to use nfmpi_get_vara_int_all() to read a 2D
! 4-byte integer array in parallel. It is the opposite of program
! 'put_vara.f'. It reads a netCDF variable defined in the input file of
! size global_nx * global_ny where
!    global_nx == 5 and
!    global_ny == (4 * number of MPI processes).
! The data partitioning pattern is a column-wise partitioning across all
! processes. Each process writes a subarray of size nx * ny.
! Note the description above follows the Fortran array index order.
!
! Example commands for MPI run.
!
!    % mpif77 -O2 -o get_vara get_vara.f -lpnetcdf
!    % mpiexec -n 4 ./get_vara /pvfs2/wkliao/testfile.nc
!
      subroutine check(err, message)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'
          integer err
          character message*(*)

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) message//' '//nfmpi_strerror(err)
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end ! subroutine check

      program main
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'

          character(LEN=256) filename, cmd
          integer i, j, err, ierr, nprocs, rank, get_args, dummy
          integer cmode, ncid, varid, dimid(2)
          integer*8 nx, ny, global_nx, global_ny
          integer*8 starts(2), counts(2)
          PARAMETER(nx=5, ny=4)
          integer buf(nx,ny)
          integer*8 malloc_size, sum_size
          logical verbose

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              verbose = .TRUE.
              filename = "testfile.nc"
              ierr = get_args(2, cmd, filename, verbose, dummy)
          endif
          call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
          if (ierr .EQ. 0) goto 999

          call MPI_Bcast(verbose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     +                   err)
          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0,
     +                   MPI_COMM_WORLD, err)

          ! open file
          err = nfmpi_open(MPI_COMM_WORLD, filename, NF_NOWRITE,
     +                     MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_open: ')

          ! inquire dimensions x and y
          err = nfmpi_inq_dimid(ncid, "x", dimid(1))
          call check(err, 'In nfmpi_inq_dim x: ')
          err = nfmpi_inq_dimid(ncid, "y", dimid(2))
          call check(err, 'In nfmpi_inq_dim y: ')

          ! get dimension lengths for dimensions Y and X */
          err = nfmpi_inq_dimlen(ncid, dimid(1), global_nx)
          call check(err, 'In nfmpi_inq_dimlen x: ')
          err = nfmpi_inq_dimlen(ncid, dimid(2), global_ny)
          call check(err, 'In nfmpi_inq_dimlen y: ')

          ! assert local array sizes
          if (global_nx .NE. nx) then
              print *, "Expect x dimension of size ",nx,
     +                  " but got ",global_nx
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          endif
          if (global_ny .NE. ny * nprocs) then
              print *, "Expect y dimension of size ",ny*nprocs,
     +                 " but got ",global_ny
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          endif

          ! get the variable ID
          err = nfmpi_inq_varid(ncid, "var", varid)
          call check(err, 'In nfmpi_inq_varid: ')

          do i=1, ny
          do j=1, nx
             buf(j,i) = -1
          enddo
          enddo

          ! Note that in Fortran, array indices start with 1
          starts(1) = 1
          starts(2) = ny * rank + 1
          counts(1) = nx
          counts(2) = ny

          err = nfmpi_get_vara_int_all(ncid, varid, starts, counts, buf)
          call check(err, 'In nfmpi_get_vara_int_all: ')

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

          ! check contents of read buffer
 997      format(A,I4,A,I4,A,I4,A,I4)
          do i=1, ny
          do j=1, nx
             if (buf(j,i) .NE. rank) then
                 print 997,"expect read buf(",j,",",i,") to be ",rank,
     +                     " but got ",buf(j,i)
                 stop 2
             end if
          enddo
          enddo

          ! check if there is any PnetCDF internal malloc residue
 998      format(A,I13,A)
          err = nfmpi_inq_malloc_size(malloc_size)
          if (err .EQ. NF_NOERR) then
              call MPI_Reduce(malloc_size, sum_size, 1, MPI_INTEGER8,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, err)
              if (rank .EQ. 0 .AND. sum_size .GT. 0)
     +            print 998,
     +            'heap memory allocated by PnetCDF internally has ',
     +            sum_size/1048576, ' MiB yet to be freed'
          endif

 999      call MPI_Finalize(err)
          ! call EXIT(0) ! EXIT() is a GNU extension
      end ! program main

