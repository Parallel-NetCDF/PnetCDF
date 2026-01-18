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
  integer :: nmode, err, ierr, fh, my_rank, nprocs, i, varid, get_args
  integer :: dimid(3)
  integer(KIND=MPI_OFFSET_KIND) :: start(3), count(3)
  real :: f(3)
  character(LEN=256) out_path, in_path, cmd, msg
  logical keep_files
  double precision timing

  call MPI_Init(ierr)

  timing = MPI_Wtime()

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  ! take out_path from command-line argument if there is any
  if (my_rank .EQ. 0) then
      out_path = FILE_NAME
      err = get_args(cmd, out_path, in_path, keep_files)
  endif
  call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (err .EQ. 0) goto 999

  call MPI_Bcast(out_path, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(keep_files, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

  ! if (nprocs .ne. 8 .AND. my_rank .eq. 0) then
  !    print *, 'Warning: ',trim(cmd),' is design to run on 8 processes.'
  ! endif

  nmode = IOR(NF90_CLOBBER,NF90_64BIT_DATA)

  call handle_err(nf90mpi_create(MPI_COMM_WORLD, out_path, nmode, MPI_INFO_NULL, fh))

  call handle_err(nf90mpi_def_dim(fh, 'dim1', 6_MPI_OFFSET_KIND, dimid(1)))
  call handle_err(nf90mpi_def_dim(fh, 'dim2', 4_MPI_OFFSET_KIND, dimid(2)))
  call handle_err(nf90mpi_def_dim(fh, 'dim3', 1_MPI_OFFSET_KIND, dimid(3)))

  call handle_err(nf90mpi_def_var(fh, 'var1', NF90_DOUBLE, dimid, varid))

  ! fill with default fill value
  call handle_err(nf90mpi_def_var_fill(fh, varid, 0, NF90_FILL_DOUBLE))

  call handle_err(nf90mpi_enddef(fh))

  do i=1,3
     f(i) = my_rank*3+i
  end do

  count = (/3,1,1/)
  start(1) = mod(my_rank,2)*3+1
  start(2) = my_rank/2+1
  start(3) = 1
  if (my_rank .GE. 8) then
      start = 1
      count = 0
  endif

  call handle_err(nf90mpi_put_var_all(fh, varid, f,start=start,count=count))

  call handle_err(nf90mpi_close(fh))

  ! Reopen the file and check it.
  call handle_err(nf90mpi_open(MPI_COMM_WORLD, out_path, NF90_NOWRITE, MPI_INFO_NULL, fh))

  call handle_err(nf90mpi_get_var_all(fh, varid, f, start=start, count=count))

  if (my_rank .LE. 8) then
     do i=1,3
        if (f(i) .ne. my_rank*3+i) then
           print *, 'Error: unexpected read value ',f(i),' should be ', my_rank*3+i
           goto 999
        endif
     end do
  endif

  call handle_err(nf90mpi_close(fh))

  999 timing = MPI_Wtime() - timing
  call MPI_Allreduce(MPI_IN_PLACE, timing, 1, &
                     MPI_DOUBLE_PRECISION, MPI_MAX, &
                     MPI_COMM_WORLD, ierr)

  if (my_rank .eq. 0) then
      if (.NOT. keep_files) then
          err = nf90mpi_delete(out_path, MPI_INFO_NULL)
      end if

      msg = '*** TESTING F90 '//trim(cmd)
      call pass_fail(0, msg, timing)
  endif

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

