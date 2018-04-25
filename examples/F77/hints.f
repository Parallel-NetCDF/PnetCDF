!
!  Copyright (C) 2012, Northwestern University and Argonne National
!  Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This example sets two PnetCDF hints:
!    nc_header_align_size and nc_var_align_size
! and prints the hint values as well as the header size, header extent, and
! two variables' starting file offsets.
!
! The compile and run commands are given below.
!
!    % mpif77 -O2 -o hints hints.f -lpnetcdf
!
!    % mpiexec -n 4 ./hints /pvfs2/wkliao/testfile.nc
!
!    nc_header_align_size      set to = 1024
!    nc_var_align_size         set to = 512
!    nc_header_read_chunk_size set to = 256
!    header size                      = 252
!    header extent                    = 1024
!    var_zy start file offset         = 1024
!    var_yx start file offset         = 3072
!


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
          end if
      end ! subroutine check

      subroutine print_hints(ncid, varid0, varid1)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          integer ncid, varid0, varid1

          ! local variables
          character*(MPI_MAX_INFO_VAL) value
          integer err, len, info_used
          logical flag
          integer*8 header_size, header_extent
          integer*8 var_zy_start, var_yx_start
          integer*8 h_align, v_align, h_chunk

          h_align=-1
          v_align=-1
          h_chunk=-1

          err = nfmpi_inq_header_size  (ncid, header_size)
          call check(err, 'In nfmpi_inq_header_size: ')
          err = nfmpi_inq_header_extent(ncid, header_extent)
          call check(err, 'In nfmpi_inq_header_extent: ')
          err = nfmpi_inq_varoffset(ncid, varid0, var_zy_start)
          call check(err, 'In nfmpi_inq_varoffset varid0: ')
          err = nfmpi_inq_varoffset(ncid, varid1, var_yx_start)
          call check(err, 'In nfmpi_inq_varoffset varid1: ')

          err = nfmpi_inq_file_info(ncid, info_used)
          call check(err, 'In nfmpi_inq_file_info : ')

          call MPI_Info_get_valuelen(info_used, "nc_header_align_size",
     +                               len, flag, err)
          if (flag) then
              call MPI_Info_get(info_used, "nc_header_align_size",
     +                          len+1, value, flag, err)
              read(value, '(i16)') h_align
          endif
          call MPI_Info_get_valuelen(info_used, "nc_var_align_size",
     +                               len, flag, err)
          if (flag) then
              call MPI_Info_get(info_used, "nc_var_align_size", len+1,
     +                          value, flag, err)
              read(value, '(i16)') v_align
          endif
          call MPI_Info_get_valuelen(info_used,
     +                               "nc_header_read_chunk_size",
     +                               len, flag, err)
          if (flag) then
              call MPI_Info_get(info_used, "nc_header_read_chunk_size",
     +                          len+1, value, flag, err)
              read(value, '(i16)') h_chunk
          endif
          call MPI_Info_free(info_used, err)

          if (h_align .EQ. -1) then
              print*,"nc_header_align_size      is NOT set"
          else
              print*,"nc_header_align_size      set to = ", h_align
          endif
          if (v_align .EQ. -1) then
              print*,"nc_var_align_size         is NOT set"
          else
              print*,"nc_var_align_size         set to = ", v_align
          endif
          if (h_chunk .EQ. -1) then
              print*,"nc_header_read_chunk_size is NOT set"
          else
              print*,"nc_header_read_chunk_size set to = ", h_chunk
          endif

          print*,"header size                      = ", header_size
          print*,"header extent                    = ", header_extent
          print*,"var_zy start file offset         = ", var_zy_start
          print*,"var_yx start file offset         = ", var_yx_start
      end ! subroutine print_hints

      subroutine put_vara(ncid, varid, ny, nx, start, count)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          integer ncid, varid, ny, nx
          integer*8 start(2), count(2)

          integer i, j, rank, err
          integer buf(nx, ny)

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          do j=1, ny
             do i=1, nx
                buf(i,j) = rank ! j*nx + i
             enddo
          enddo

          err = nfmpi_put_vara_int_all(ncid, varid, start, count, buf)
          call check(err, 'In nfmpi_put_vara_int_all : ')

      end ! subroutine put_vara

      program main
      implicit none
      include "mpif.h"
      include "pnetcdf.inc"

      character*256 filename, cmd
      integer NZ, NY, NX
      integer ncid, rank, nprocs, info, cmode, err, ierr, get_args
      integer varid0, varid1, dimid(3), dimid2(2)
      integer*8 start(2), count(2), llen
      integer*8 malloc_size, sum_size
      logical verbose
      integer dummy

      call MPI_Init(err)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

      NZ = 5
      NY = 5
      NX = 5

      if (rank .EQ. 0) then
          verbose = .TRUE.
          filename = "testfile.nc"
          ierr = get_args(2, cmd, filename, verbose, dummy)
      endif
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      if (ierr .EQ. 0) goto 999

      call MPI_Bcast(verbose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     +               err)
      call MPI_Bcast(filename, 256, MPI_CHARACTER, 0,
     +               MPI_COMM_WORLD, err)

      call MPI_Info_create(info, err)
      call MPI_Info_set(info, "nc_header_align_size",      "1024", err)
      call MPI_Info_set(info, "nc_var_align_size",         "512",  err)
      call MPI_Info_set(info, "nc_header_read_chunk_size", "256",  err)
      ! note that set the above values to 1 to disable the alignment

      ! create a new file for writing
      cmode = NF_CLOBBER + NF_64BIT_DATA
      err = nfmpi_create(MPI_COMM_WORLD, filename, cmode, info, ncid)
      call check(err, 'In nfmpi_create : ')
      call MPI_Info_free(info, err)

      ! define 3 dimensions
      llen = NZ*nprocs
      err = nfmpi_def_dim(ncid, "Z", llen, dimid(1))
      call check(err, 'In nfmpi_def_dim Z : ')
      llen = NY*nprocs
      err = nfmpi_def_dim(ncid, "Y", llen, dimid(2))
      call check(err, 'In nfmpi_def_dim Y : ')
      llen = NX*nprocs
      err = nfmpi_def_dim(ncid, "X", llen, dimid(3))
      call check(err, 'In nfmpi_def_dim X : ')

      ! define a variable of size (NZ * nprocs) * (NY * nprocs)
      dimid2(1) = dimid(1)
      dimid2(2) = dimid(2)
      err = nfmpi_def_var(ncid, "var_zy", NF_INT,   2, dimid2, varid0)
      call check(err, 'In nfmpi_def_var var_zy : ')
      ! define a variable of size (NY * nprocs) * (NX * nprocs)
      dimid2(1) = dimid(2)
      dimid2(2) = dimid(3)
      err = nfmpi_def_var(ncid, "var_yx", NF_FLOAT, 2, dimid2, varid1)
      call check(err, 'In nfmpi_def_var var_yx : ')

      err = nfmpi_enddef(ncid)
      call check(err, 'In nfmpi_enddef : ')

      ! var_zy is partitioned along Z dimension
      start(1) = NZ * rank + 1
      start(2) = 1
      count(1) = NZ
      count(2) = NY * nprocs
      call put_vara(ncid, varid0, NZ, NY * nprocs, start, count)

      ! var_yx is partitioned along X dimension
      start(1) = 1
      start(2) = NX * rank + 1
      count(1) = NY * nprocs
      count(2) = NX
      call put_vara(ncid, varid1, NY * nprocs, NX, start, count)

      if (rank .EQ. 0 .AND. verbose)
     +    call print_hints(ncid, varid0, varid1)

      err = nfmpi_close(ncid)

      ! check if there is any PnetCDF internal malloc residue
 998  format(A,I13,A)
      err = nfmpi_inq_malloc_size(malloc_size)
      if (err .EQ. NF_NOERR) then
          call MPI_Reduce(malloc_size, sum_size, 1, MPI_INTEGER8,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, err)
      if (rank .EQ. 0 .AND. sum_size .GT. 0)
     +        print 998,
     +        'heap memory allocated by PnetCDF internally has ',
     +        sum_size/1048576, ' MiB yet to be freed'
      endif

 999  call MPI_Finalize(err)
      ! call EXIT(0) ! EXIT() is a GNU extension
      end ! program

