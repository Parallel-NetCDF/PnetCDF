!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

! This program tests large files (> 4 GB)

program tst_flarge
  use mpi
  use pnetcdf
  implicit none

  integer :: ncFileID, dimID, varID1, varID2 
  integer(KIND=MPI_OFFSET_KIND) :: BIG_DIMENSION = 300000000
  character (len = *), parameter :: FILE_NAME = "tst_flarge.nc"
  character (len = *), parameter :: dimName = "really_big_dimension"
  character (len = *), parameter :: var1Name = "TweedleDum"
  character (len = *), parameter :: var2Name = "TweedleDee"
  double precision, parameter :: VAL1 = 42.5
  double precision, parameter :: VAL2 = -42.5
  double precision :: val1_in
  double precision :: val2_in
  integer :: cmode, err
  double precision dbl_buf(1)
  integer(KIND=MPI_OFFSET_KIND) :: start(1), count(1)
  character(LEN=128) filename, cmd, msg
  integer argc, iargc, my_rank, p

  call MPI_Init(err)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, err)
  call MPI_Comm_size(MPI_COMM_WORLD, p, err)

  ! take filename from command-line argument if there is any
  call getarg(0, cmd)
  argc = IARGC()
  if (argc .GT. 1) then
     if (my_rank .EQ. 0) print*,'Usage: ',trim(cmd),' [filename]'
     goto 999
  endif
  filename = FILE_NAME
  if (argc .EQ. 1) call getarg(1, filename)

  if (p .ne. 1 .AND. my_rank .eq. 0) then
     print *, 'Warning: ',trim(cmd),' is design to run on 1 process'
  endif

  ! Create the file with 2 NF90_DOUBLE vars, each with one really long dimension.
  cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
  call check(nf90mpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, ncFileID))
  call check(nf90mpi_def_dim(ncFileID, dimName, BIG_DIMENSION, dimID))
  call check(nf90mpi_def_var(ncFileID, var1Name, nf90_double, (/ dimID /), varID1) )
  call check(nf90mpi_def_var(ncFileID, var2Name, nf90_double, (/ dimID /), varID2) )

  call check(nf90mpi_enddef(ncFileID))
  call check(nf90mpi_begin_indep_data(ncFileID))

  ! Write a value in each variable.
  dbl_buf(1) = 42.5
  start(1) = 1
  count(1) = 1
  call check(nf90mpi_put_var(ncFileID, VarID1, dbl_buf, start, count))
  dbl_buf(1) = -42.5
  start(1) = BIG_DIMENSION
  count(1) = 1
  call check(nf90mpi_put_var(ncFileID, VarID2, dbl_buf, start, count))

  call check(nf90mpi_close(ncFileID))

  ! Now open the file to read and check a few values
  call check(nf90mpi_open(MPI_COMM_WORLD, filename, NF90_NOWRITE, MPI_INFO_NULL, ncFileID))
  call check(nf90mpi_begin_indep_data(ncFileID))
  start(1) = 1
  call check(nf90mpi_get_var(ncFileID, VarID1, val1_in, start))
  start(1) = BIG_DIMENSION
  call check(nf90mpi_get_var(ncFileID, VarID2, val2_in, start))
  if(val1_in /= VAL1 .or. val2_in /= VAL2) then
     print *, 'Variable value not what was written'
     stop 2
  end if

  call check(nf90mpi_close(ncFileID))

   msg = '*** TESTING F90 '//trim(cmd)//' for large files'
   if (my_rank .eq. 0)   write(*,"(A67,A)") msg,'------ pass'

 999  call MPI_Finalize(err)

contains
  ! Internal subroutine - checks error status after each netcdf, prints out text message each time
  !   an error code is returned. 
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90mpi_strerror(status))
       stop 2
    end if
  end subroutine check
end program tst_flarge
