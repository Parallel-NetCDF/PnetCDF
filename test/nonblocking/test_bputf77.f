!
!   Copyright (C) 2012, Northwestern University and Argonne National Lab
!   See COPYRIGHT notice in top-level directory.
!
!   $Id$
!

      INTEGER FUNCTION XTRIM(STRING)
           CHARACTER*(*) STRING
           INTEGER I, N
           N = LEN(STRING)
           DO I = N, 1, -1
              IF (STRING(I:I) .NE. ' ') GOTO 10
           ENDDO
 10        XTRIM = I
      END ! FUNCTION XTRIM

      subroutine check(err, message)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          integer err
          character message*(*)

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) message//' '//nfmpi_strerror(err)
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          endif
      end ! subroutine check

      program main
      implicit none
      include "mpif.h"
      include "pnetcdf.inc"

      logical verbose
      integer i, j, ncid, varid, err, ierr, rank, nprocs, info
      integer no_err, cmode, get_args, XTRIM
      integer dimid(2), req(2), status(2)
      integer*8 start(2)
      integer*8 count(2)
      integer*8 stride(2)
      integer*8 imap(2)
      integer*8 bufsize, dim_size
      real  var(6,4)
      character*256 filename, cmd, msg

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if (rank .EQ. 0) then
          filename = "testfile.nc"
          err = get_args(cmd, filename)
      endif
      call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (err .EQ. 0) goto 999

      call MPI_Bcast(filename, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD,
     +               ierr)

      verbose = .FALSE.
      if (nprocs .GT. 1 .AND. rank .EQ. 0 .AND. verbose) then
          print*,'Warning: ',cmd(1:XTRIM(cmd)),
     +           ' is designed to run on 1 process'
      endif

      call MPI_Info_create(info, ierr)
      ! call MPI_Info_set(info, "romio_pvfs2_posix_write","enable",ierr)

      cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
      err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                   info, ncid)
      call check(err, 'Error at nfmpi_create ')

      call MPI_Info_free(info, ierr)

      ! define a variable of a 4 x 6 integer array in the nc file
      dim_size = 4
      err = nfmpi_def_dim(ncid, 'X', dim_size, dimid(1))
      call check(err, 'Error at nfmpi_def_dim ')

      dim_size = 6
      err = nfmpi_def_dim(ncid, 'Y', dim_size, dimid(2))
      call check(err, 'Error at nfmpi_def_dim ')

      err = nfmpi_def_var(ncid, 'var', NF_INT64, 2, dimid, varid)
      call check(err, 'Error at nfmpi_def_var ')

      err = nfmpi_enddef(ncid)
      call check(err, 'Error at nfmpi_enddef ')

      ! set the contents of write buffer var, a 6 x 4 real array
      !     50, 56, 62, 68,
      !     51, 57, 63, 69,
      !     52, 58, 64, 70,
      !     53, 59, 65, 71,
      !     54, 60, 66, 72,
      !     55, 61, 67, 73
      do j = 1, 4
         do i = 1, 6
            var(i,j) = (j-1)*6+(i-1) + 50
         enddo
      enddo

      ! bufsize must be max of data type converted before and after
      bufsize = 4*6*8
      err = nfmpi_buffer_attach(ncid, bufsize)
      call check(err, 'Error at nfmpi_buffer_attach ')

      ! write var to the NC variable in the matrix transposed way
      count(1)  = 2
      count(2)  = 6
      stride(1) = 1
      stride(2) = 1
      imap(1)   = 6
      imap(2)   = 1   ! imap would be {1, 4} if not transposing

      if (rank .GT. 0) then
         count(1)  = 0
         count(2)  = 0
      endif

      ! write the first two columns of the NC variable in the matrix transposed way
      start(1)  = 1
      start(2)  = 1
      err = nfmpi_bput_varm_real(ncid, varid, start, count, stride,
     +                           imap, var(1,1), req(1))
      call check(err, 'Error at nfmpi_bput_varm_real ')

      ! write the second two columns of the NC variable in the matrix transposed way
      start(1)  = 3
      start(2)  = 1
      err = nfmpi_bput_varm_real(ncid, varid, start, count, stride,
     +                           imap, var(1,3), req(2))
      call check(err, 'Error at nfmpi_bput_varm_real ')

      err = nfmpi_wait_all(ncid, 2, req, status)
      call check(err, 'Error at nfmpi_wait_all ')

      ! check each bput status
      do i = 1, 2
          if (status(i) .ne. NF_NOERR) then
              print*,'Error at bput status ', nfmpi_strerror(status(i))
          endif
      enddo

      err = nfmpi_buffer_detach(ncid)
      call check(err, 'Error at nfmpi_buffer_detach ')

      ! the output from command "ncmpidump -v var test.nc" should be:
      !      var =
      !       50, 56, 62, 68,
      !       51, 57, 63, 69,
      !       52, 58, 64, 70,
      !       53, 59, 65, 71,
      !       54, 60, 66, 72,
      !       55, 61, 67, 73 ;
      ! note that the display of ncmpidump is in C array dimensional order

      ! check if the contents of write buffer have been altered (should not be)
      no_err = 0
      if (rank .EQ. 0) then
         do j = 1, 4
            do i = 1, 6
               if (var(i,j) .NE. (j-1)*6+(i-1) + 50) then
! #ifdef PRINT_ERR_ON_SCREEN
!                  ! this error is a pnetcdf internal error, if occurs */
!                  print*, &
!                  'Error: bput_varm write buffer has been altered at j=', &
!                  j,' i=',i,' var=',var(i,j)
! #endif
                   no_err = no_err + 1
                endif
            enddo
         enddo
      endif

      err = nfmpi_close(ncid)
      call check(err, 'Error at nfmpi_close ')

      if (rank .EQ. 0) then
         msg = '*** TESTING F77 '//cmd(1:XTRIM(cmd))//
     +         ' for bput_varm_real API'
         call pass_fail(no_err, msg)
      endif

 999  CALL MPI_Finalize(ierr)

      end ! program

