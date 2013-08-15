!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

!     This program tests PnetCDF parallel I/O functions from fortran 90.
!
!     We are writing 2D data, a 6 x 12 grid, on 4 processors. Each
!     processor will write it's rank to it's quarter of the array. The
!     result will be (in CDL):
!
! netcdf f90tst_parallel {
! dimensions:
! 	x = 16 ;
! 	y = 16 ;
! variables:
! 	int data(x, y) ;
! data:
! 
!  data =
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
!   2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 ;
! }

program f90tst_parallel
  use mpi
  use pnetcdf
  implicit none
  
  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "f90tst_parallel.nc"

  integer, parameter :: MAX_DIMS = 2
  integer, parameter :: NX = 16, NY = 16
  integer, parameter :: NUM_PROC = 4
  integer :: ncid, varid, dimids(MAX_DIMS), chunksizes(MAX_DIMS), chunksizes_in(MAX_DIMS)
  integer :: x_dimid, y_dimid, contig
  integer :: data_out(NY / 2, NX / 2), data_in(NY / 2, NX / 2)
  integer :: mode_flag
  integer :: nvars, ngatts, ndims, unlimdimid, file_format
  integer :: x, y, retval
  integer :: p, my_rank, ierr
  integer(KIND=MPI_OFFSET_KIND) :: start(MAX_DIMS), count(MAX_DIMS)
  integer(KIND=MPI_OFFSET_KIND) :: nx_ll, ny_ll

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

  if (my_rank .eq. 0) then
     write(*,"(A)",advance="no") '*** Testing PnetCDF parallel I/O from Fortran 90.'
  endif

  ! There must be 4 procs for this test.
  if (p .ne. 4) then
     print *, 'Sorry, this test program must be run on four processors.'
     stop 2 
  endif

  ! Create some pretend data.
  do x = 1, NX / 2
     do y = 1, NY / 2
        data_out(y, x) = my_rank
     end do
  end do

  ! Create the netCDF file. 
  mode_flag = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
  call handle_err(nf90mpi_create(MPI_COMM_WORLD, FILE_NAME, mode_flag, MPI_INFO_NULL, ncid))

  ! Define the dimensions.
  nx_ll = NX
  ny_ll = NY
  call handle_err(nf90mpi_def_dim(ncid, "x", nx_ll, x_dimid))
  call handle_err(nf90mpi_def_dim(ncid, "y", ny_ll, y_dimid))
  dimids =  (/ y_dimid, x_dimid /)

  ! Define the variable. 
  call handle_err(nf90mpi_def_var(ncid, "data", NF90_INT, dimids, varid))

  call handle_err(nf90mpi_enddef(ncid))

  ! Determine what part of the variable will be written for this
  ! processor. It's a checkerboard decomposition.
  count = (/ NX / 2, NY / 2 /)
  if (my_rank .eq. 0) then
     start = (/ 1, 1 /)
  else if (my_rank .eq. 1) then
     start = (/ NX / 2 + 1, 1 /)
  else if (my_rank .eq. 2) then
     start = (/ 1, NY / 2 + 1 /)
  else if (my_rank .eq. 3) then
     start = (/ NX / 2 + 1, NY / 2 + 1 /)
  endif

  ! Write this processor's data.
  call handle_err(nf90mpi_put_var_all(ncid, varid, data_out, start = start, count = count))

  ! Close the file. 
  call handle_err(nf90mpi_close(ncid))

  ! Reopen the file.
  call handle_err(nf90mpi_open(MPI_COMM_WORLD, FILE_NAME, nf90_nowrite, MPI_INFO_NULL, ncid))
  
  ! Check some stuff out.
  call handle_err(nf90mpi_inquire(ncid, ndims, nvars, ngatts, unlimdimid, file_format))
  if (ndims /= 2 .or. nvars /= 1 .or. ngatts /= 0 .or. unlimdimid /= -1 .or. &
       file_format /= nf90_format_64bit_data) stop 3

  ! Read this processor's data.
  call handle_err(nf90mpi_get_var_all(ncid, varid, data_in, start = start, count = count))

  ! Check the data.
  do x = 1, NX / 2
     do y = 1, NY / 2
        if (data_in(y, x) .ne. my_rank) stop 4
     end do
  end do

  ! Close the file.
  call handle_err(nf90mpi_close(ncid))

  if (my_rank .eq. 0) write(*,"(A)") '                  ------ pass'

  call MPI_Finalize(ierr)


contains
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
  subroutine handle_err(errcode)
    implicit none
    integer, intent(in) :: errcode
    
    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90mpi_strerror(errcode))
       stop 2
    endif
  end subroutine handle_err
end program f90tst_parallel

