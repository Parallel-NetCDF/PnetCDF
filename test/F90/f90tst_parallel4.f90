!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
! $Id$

! This parallel test was contributed by Jim Edwards at UCAR. Thanks Jim!
program f90tst
  use mpi
  use pnetcdf
  implicit none

  character (len = *), parameter :: FILE_NAME = "f90tst_nc4_par.nc"
  integer :: nmode, ierr, fh, my_task, nprocs, i, varid
  integer :: dimid(3)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)
  real :: f(3)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_task, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  if(nprocs/=8)then
     stop 'requires 8 tasks'
  end if

  if (my_task .eq. 0) then
     write(*,"(A)",advance="no") '*** Testing PnetCDF parallel I/O from Fortran 90.'
  endif

  nmode = ior(NF90_CLOBBER,NF90_64BIT_DATA)

  call handle_err(nf90mpi_create(MPI_COMM_WORLD, FILE_NAME, nmode, MPI_INFO_NULL, fh))

  call handle_err(nf90mpi_def_dim(fh, 'dim1', 6_8, dimid(1)))
  call handle_err(nf90mpi_def_dim(fh, 'dim2', 4_8, dimid(2)))
  call handle_err(nf90mpi_def_dim(fh, 'dim3', 1_8, dimid(3)))


  call handle_err(nf90mpi_def_var(fh, 'var1', NF90_DOUBLE, dimid, varid))
  call handle_err(nf90mpi_enddef(fh))


  do i=1,3
     f(i) = my_task*3+i
  end do

  count = (/3,1,1/)
  start(1) = mod(my_task,2)*3+1
  start(2) = my_task/2+1
  start(3) = 1

  call handle_err(nf90mpi_put_var_all(fh, varid, f,start=start,count=count))

  call handle_err(nf90mpi_close(fh))

  ! Reopen the file and check it.
  call handle_err(nf90mpi_open(MPI_COMM_WORLD, FILE_NAME, NF90_NOWRITE, MPI_INFO_NULL, fh))

  call handle_err(nf90mpi_get_var_all(fh, varid, f, start=start, count=count))
  do i=1,3
     if (f(i) .ne. my_task*3+i) stop 3
  end do

  call handle_err(nf90mpi_close(fh))

  if (my_task .eq. 0) write(*,"(A)") '                  ------ pass'

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

end program f90tst

