!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$


program tst_f90_nc4
  use mpi
  use pnetcdf
  implicit none
  integer :: fh, ierr, dimid, varid, ndim, nvar
  character (len = *), parameter :: FILE_NAME = "tst_f90_nc4.nc"
  integer :: cmode, err
  integer(KIND=MPI_OFFSET_KIND) :: ten=10

  write(*,"(A)",advance="no") '*** Testing simple netCDF file from Fortran 90.'

  call MPI_Init(err)

  cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
  call check(nf90mpi_create(MPI_COMM_WORLD, FILE_NAME, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_close(fh))
  
  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_WORLD, FILE_NAME, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

  call check(nf90mpi_create(MPI_COMM_WORLD, FILE_NAME, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_close(fh))
  
  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_WORLD, FILE_NAME, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

  write(*,"(A)") '                    ------ pass'
  call MPI_Finalize(err)

contains
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
  subroutine check(errcode)
    implicit none
    integer, intent(in) :: errcode
    
    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90mpi_strerror(errcode))
       stop 2
    endif
  end subroutine check
end program tst_f90_nc4

