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
  integer :: fh, cmode, err, ierr, dimid, varid, ndim, nvar, get_args
  character (len = *), parameter :: FILE_NAME = "tst_f90_nc4.nc"
  integer(KIND=MPI_OFFSET_KIND) :: ten=10
  character(LEN=256) out_path, in_path, cmd, msg
  integer my_rank, nprocs, fillmode
  logical keep_files
  double precision timing

  call MPI_Init(ierr)

  timing = MPI_Wtime()

  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! take out_path from command-line argument if there is any
  cmd = ' '
  if (my_rank .EQ. 0) then
      out_path = FILE_NAME
      err = get_args(cmd, out_path, in_path, keep_files)
  endif
  call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (err .EQ. 0) goto 999

  call MPI_Bcast(out_path, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(keep_files, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

  if (my_rank > 0) goto 999;

!  if (nprocs .ne. 1 .AND. my_rank .eq. 0) then
!     print *, 'Warning: ',trim(cmd),' is design to run on 1 process'
!  endif

  cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
  call check(nf90mpi_create(MPI_COMM_SELF, out_path, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_set_fill(fh, NF90_FILL, fillmode))
  call check(nf90mpi_close(fh))

  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_SELF, out_path, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

  call check(nf90mpi_create(MPI_COMM_SELF, out_path, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_set_fill(fh, NF90_FILL, fillmode))
  call check(nf90mpi_close(fh))

  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_SELF, out_path, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

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
  end if

  call MPI_Finalize(ierr)

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

