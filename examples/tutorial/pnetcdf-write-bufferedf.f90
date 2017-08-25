!
!   Copyright (C) 2012, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

      subroutine check(err, message)
          use mpi
          use pnetcdf
          implicit none
          integer err
          character(len=*) message

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF90_NOERR) then
              write(6,*) trim(message), trim(nf90mpi_strerror(err))
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end subroutine check

      program main
      use mpi
      use pnetcdf
      implicit none

      integer i, j, ncid, varid, cmode, err, rank, nprocs
      integer ierr, get_args, dummy
      integer dimid(2), req(2), status(2)
      integer(kind=MPI_OFFSET_KIND) start(2)
      integer(kind=MPI_OFFSET_KIND) count(2)
      integer(kind=MPI_OFFSET_KIND) stride(2)
      integer(kind=MPI_OFFSET_KIND) imap(2)
      integer(kind=MPI_OFFSET_KIND) bufsize
      integer(kind=MPI_OFFSET_KIND) put_size
      real  var(6,4)
      character(len=256) filename, cmd
      logical verbose

      call MPI_INIT(err)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, err)

      ! take filename from command-line argument if there is any
      if (rank .EQ. 0) then
          filename = "testfile.nc"
          ierr = get_args(2, cmd, filename, verbose, dummy)
      endif
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      if (ierr .EQ. 0) goto 999

      call MPI_Bcast(filename, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, err)

      cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
      err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
                         MPI_INFO_NULL, ncid)
      call check(err, 'Error at nf90mpi_create ')

      ! define a variable of a (4*nprocs) x 6 integer array in the nc file
      err = nf90mpi_def_dim(ncid, 'X', 4_MPI_OFFSET_KIND*nprocs, dimid(1))
      call check(err, 'Error at nf90mpi_def_dim ')

      err = nf90mpi_def_dim(ncid, 'Y', 6_MPI_OFFSET_KIND, dimid(2))
      call check(err, 'Error at nf90mpi_def_dim ')

      err = nf90mpi_def_var(ncid, 'var', NF90_INT64, dimid, varid)
      call check(err, 'Error at nf90mpi_def_var ')

      err = nf90mpi_enddef(ncid)
      call check(err, 'Error at nf90mpi_enddef ')

      ! set the contents of the local write buffer var, a 4 x 6 real array
      ! for example, for rank == 2, var(4,6) =
      !     48, 54, 60, 65,
      !     49, 55, 61, 67,
      !     50, 56, 62, 68,
      !     51, 57, 63, 69,
      !     52, 58, 64, 70,
      !     53, 59, 65, 71
      do j = 1, 4
         do i = 1, 6
            var(i,j) = (j-1)*6+(i-1) + rank*24
         enddo
      enddo

      ! bufsize must be max of data type converted before and after
      bufsize = 4*6*8
      err = nf90mpi_buffer_attach(ncid, bufsize)
      call check(err, 'Error at nf90mpi_buffer_attach ')

      ! write var to the NC variable in the matrix transposed way
      count(1)  = 2
      count(2)  = 6
      stride(1) = 1
      stride(2) = 1
      imap(1)   = 6
      imap(2)   = 1

      req(:) = NF90_REQ_NULL  ! actually not necessary, added for testing

      ! write to the 1st two columns of the variable in matrix transposed way
      start(1)  = 1 + rank*4
      start(2)  = 1
      err = nf90mpi_bput_var(ncid, varid, var(1:,1:), req(1), &
                             start, count, stride, imap)
      call check(err, 'Error at nf90mpi_bput_var ')

      ! write to the 2nd two columns of the variable in transposed way
      start(1)  = 3 + rank*4
      start(2)  = 1
      err = nf90mpi_bput_var(ncid, varid, var(1:,3:), req(2), &
                             start, count, stride, imap)
      call check(err, 'Error at nf90mpi_bput_var ')

      err = nf90mpi_wait_all(ncid, 2, req, status)
      call check(err, 'Error at nf90mpi_wait_all ')

      ! check each bput status
      do i = 1, 2
          if (status(i) .ne. NF90_NOERR) then
              print*,'Error at bput status ', &
                     nf90mpi_strerror(status(i))
              stop 2
          endif
      enddo

      err = nf90mpi_buffer_detach(ncid)
      call check(err, 'Error at nf90mpi_buffer_detach ')

      ! The output from command "ncmpidump test.nc" is shown below if run
      ! this example on 4 processes.
      !
      ! netcdf test {
      ! // file format: CDF-5 (big variables)
      ! dimensions:
      !        Y = 6 ;
      !        X = 16 ;
      ! variables:
      !        int64 var(Y, X) ;
      !data:
      !
      ! var =
      !  0,  6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90,
      !  1,  7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85, 91,
      !  2,  8, 14, 20, 26, 32, 38, 44, 50, 56, 62, 68, 74, 80, 86, 92,
      !  3,  9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93,
      !  4, 10, 16, 22, 28, 34, 40, 46, 52, 58, 64, 70, 76, 82, 88, 94,
      !  5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95 ;
      !
      ! note that the display of ncmpidump is in C array dimensional order

      err = nf90mpi_inq_put_size(ncid, put_size)
      call check(err, 'Error at nf90mpi_inq_put_size ')
      ! print*,'pnetcdf reports total put size by this proc =', put_size

      err = nf90mpi_close(ncid)
      call check(err, 'Error at nf90mpi_close ')

 999  CALL MPI_Finalize(err)
      end program

