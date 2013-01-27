!
!   Copyright (C) 2012, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!

      program main

      use pnetcdf
      implicit none

#include "mpif.h"
      integer i, j, ncid, varid, retval, err, rank, nprocs
      integer no_err, cmode
      integer dimid(2)
      integer req(2), status(2)
      integer(kind=MPI_OFFSET_KIND) start(2)
      integer(kind=MPI_OFFSET_KIND) count(2)
      integer(kind=MPI_OFFSET_KIND) stride(2)
      integer(kind=MPI_OFFSET_KIND) imap(2)
      integer(kind=MPI_OFFSET_KIND) bufsize
      real  var(6,4)

      call MPI_INIT(err)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, err)

      if (nprocs > 1) then
          print*,"This test is designed to run on one process"
          if (rank > 0) then
              call MPI_Finalize(err)
              Stop
          endif
      endif

      cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
      err = nfmpi_create(MPI_COMM_WORLD, 'test.nc', cmode,
     &                   MPI_INFO_NULL, ncid)
      if (err < NF_NOERR) print*,'Error at nfmpi_create ',
     &                           nfmpi_strerror(err)

      ! define a variable of a 4 x 6 integer array in the nc file
      err = nfmpi_def_dim(ncid, 'X', 4_MPI_OFFSET_KIND, dimid(1))
      if (err < NF_NOERR) print*,'Error at nfmpi_def_dim ',
     &                           nfmpi_strerror(err)

      err = nfmpi_def_dim(ncid, 'Y', 6_MPI_OFFSET_KIND, dimid(2))
      if (err < NF_NOERR) print*,'Error at nfmpi_def_dim ',
     &                           nfmpi_strerror(err)

      err = nfmpi_def_var(ncid, 'var', NF_INT64, 2, dimid, varid)
      if (err < NF_NOERR) print*,'Error at nfmpi_def_var ',
     &                           nfmpi_strerror(err)

      err = nfmpi_enddef(ncid)
      if (err < NF_NOERR) print*,'Error at nfmpi_enddef ',
     &                           nfmpi_strerror(err)

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
      if (err < NF_NOERR) print*,'Error at nfmpi_buffer_attach ',
     &                           nfmpi_strerror(err)

      ! write var to the NC variable in the matrix transposed way
      count(1)  = 2
      count(2)  = 6
      stride(1) = 1
      stride(2) = 1
      imap(1)   = 6
      imap(2)   = 1   ! imap would be {1, 4} if not transposing

      ! write the first two columns of the NC variable in the matrix transposed way
      start(1)  = 1
      start(2)  = 1
      err = nfmpi_bput_varm_real(ncid, varid, start, count, stride,
     &                           imap, var(1,1), req(1))
      if (err < NF_NOERR) print*,'Error at nfmpi_bput_varm_real ',
     &                           nfmpi_strerror(err)

      ! write the second two columns of the NC variable in the matrix transposed way
      start(1)  = 3
      start(2)  = 1
      err = nfmpi_bput_varm_real(ncid, varid, start, count, stride,
     &                           imap, var(1,3), req(2))
      if (err < NF_NOERR) print*,'Error at nfmpi_bput_varm_real ',
     &                           nfmpi_strerror(err)

      err = nfmpi_wait_all(ncid, 2, req, status)
      if (err < NF_NOERR) print*,'Error at nfmpi_wait_all ',
     &                           nfmpi_strerror(err)

      ! check each bput status
      do i = 1, 2
          if (status(i) .ne. NF_NOERR) then
              print*,'Error at bput status ',
     &                           nfmpi_strerror(status(i))
          endif
      enddo

      err = nfmpi_buffer_detach(ncid)
      if (err < NF_NOERR) print*,'Error at nfmpi_buffer_detach ',
     &                           nfmpi_strerror(err)

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
      do j = 1, 4
         do i = 1, 6
            if (var(i,j) .NE. (j-1)*6+(i-1) + 50) then
#ifdef PRINT_ERR_ON_SCREEN
                ! this error is a pntecdf internal error, if occurs */
                print*,
                'Error: bput_varm write buffer has been altered at j=',
     &          j,' i=',i,' var=',var(i,j)
#endif
                no_err = no_err + 1
             endif
         enddo
      enddo

      err = nfmpi_close(ncid)
      if (err < NF_NOERR) print*,'Error at nfmpi_close ',
     &                           nfmpi_strerror(err)

      CALL MPI_Finalize(err)

      if (no_err .GT. 0) then
          print*,'*** Testing Fortran bput_varm failed'
      else
          print*,'*** Testing Fortran bput_varm: Success'
      endif
      Stop
      end program

