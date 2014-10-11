!
!   Copyright (C) 2014, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This program tests the vard API.
! The write buffer is a 2D array of size NY x NX
! The MPI data type for the buffer is defined by swapping the 1st and 2nd
! rows of the array using a butype constructed by MPI_Type_create_hindex().
! It also writes a fixed-size variable using a buftype constructed by
! MPI_Type_create_subarray(). Both record and foxed-size variables are read
! back using various filetypes and buftypes and check the contents.
!
! The expected reults from the output file contents are:
! (when running on 1 MPI process)
!
!  % ncmpidump testfile.nc
!    netcdf testfile {
!    // file format: CDF-1
!    dimensions:
!           REC_DIM = UNLIMITED ; // (2 currently)
!           X = 5 ;
!           FIX_DIM = 2 ;
!    variables:
!           int rec_var(REC_DIM, X) ;
!           int dummy_rec(REC_DIM, X) ;
!           int fix_var(FIX_DIM, X) ;
!    data:
! 
!    rec_var =
!      10, 11, 12, 13, 14,
!      0, 1, 2, 3, 4 ;
! 
!    dummy_rec =
!      0, 0, 0, 0, 0,
!      0, 0, 0, 0, 0 ;
!
!    fix_var =
!      10, 11, 12, 13, 14,
!      0, 1, 2, 3, 4 ;
! }
!
      subroutine check(err, message)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          integer err
          character(len=*) message

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) trim(message), trim(nfmpi_strerror(err))
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end subroutine check

      subroutine check_value(rank, nx, ny, buf, nerrs)
          implicit none
          include "mpif.h"
          integer(kind=MPI_OFFSET_KIND) NX, NY
          integer rank, buf(nx, ny), nerrs
          integer i, j
 543      format(A,I2,A,I2,A,I3,A,I3)
          do j=1, ny
             do i=1, nx
                if (buf(i,j) .NE. rank*100+j*10+i) then
                    print 543,'expecting buf(',j,',',i,')=',
     +              rank*100+j*10+i,' but got ',buf(j,i)
                    nerrs = nerrs + 1
                endif
             enddo
          enddo
      end subroutine check_value

      subroutine check_value_permuted(rank, nx, ny, buf, nerrs)
          implicit none
          include "mpif.h"
          integer(kind=MPI_OFFSET_KIND) NX, NY
          integer rank, buf(nx, ny), nerrs, val
          integer i, j
 543      format(A,I2,A,I2,A,I3,A,I3)
          do j=1, ny
             val = 1
             if (j .EQ. 1) val = 2
             val = rank * 100 + val * 10
             do i=1, nx
                if (buf(i,j) .NE. val+i) then
                    print 543,'expecting buf(',i,',',j,')=',
     +              val+i,' but got ',buf(j,i)
                    nerrs = nerrs + 1
                endif
             enddo
          enddo
      end subroutine check_value_permuted

      subroutine clear_buf(nx, ny, buf)
          implicit none
          include "mpif.h"
          integer(kind=MPI_OFFSET_KIND) NX, NY
          integer buf(nx, ny)
          integer i, j
          do j=1, NY
             do i=1, NX
                buf(j,i) = -1
             enddo
          enddo
      end subroutine clear_buf

      subroutine get_var_and_verify(ncid,varid,NX,NY,start,count,buf,
     +                       buftype, ghost_buftype, filetype, nerrs)
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          integer(kind=MPI_OFFSET_KIND) NX, NY, start(2), count(2)
          integer i, j, rank, err, ncid, varid, nerrs, buf(NX, NY)
          integer buftype, ghost_buftype, filetype
          integer ncbuf(NX+4, NY+4)

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          ! clear the contents of the read buffer
          call clear_buf(NX, NY, buf)

          ! read back using regular vara API
          err = nfmpi_get_vara_int_all(ncid, varid, start, count, buf)
          call check(err, 'In nfmpi_get_vara_int_all: ')

          ! check if the contents of buf are expected
          call check_value_permuted(rank, NX, NY, buf, nerrs)
          call clear_buf(NX, NY, buf)

          ! read back using flexible vara API
          err = nfmpi_get_vara_all(ncid, varid, start, count, buf(:,2),
     +                             1_8, buftype)
          call check(err, 'In nfmpi_get_vara_all: ')

          ! check if the contents of buf are expected
          call check_value(rank, NX, NY, buf, nerrs)
          call clear_buf(NX, NY, buf)

          ! read back using vard API and permuted buftype
          err = nfmpi_get_vard_all(ncid, varid, filetype, buf(:,2),
     +                             1_8, buftype)
          call check(err, 'In nfmpi_get_vard_all: ')

          ! check if the contents of buf are expected
          call check_value(rank, NX, NY, buf, nerrs)
          call clear_buf(NX, NY, buf)

          ! read back using vard API and no buftype
          err = nfmpi_get_vard_all(ncid, varid, filetype, buf, 0_8,
     +                             MPI_DATATYPE_NULL)
          call check(err, 'In nfmpi_get_vard_all: ')

          ! check if the contents of buf are expected
          call check_value_permuted(rank, NX, NY, buf, nerrs)

          ! read back using ghost buftype
          err = nfmpi_get_vard_all(ncid, varid, filetype, ncbuf, 1_8,
     +                             ghost_buftype)
          call check(err, 'In nfmpi_get_vard_all: ')

 543      format(A,I2,A,I2,A,I3,A,I3)
          do j=1, NY
             do i=1, NX
                if (buf(i,j) .NE. ncbuf(i+2, j+2)) then
                    print 543,'expecting buf(',i,',',j,')=',
     +              buf(i,j),' but got ',ncbuf(i+2,j+2)
                    nerrs = nerrs + 1
                endif
             enddo
          enddo
      end subroutine get_var_and_verify

      program main
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"
          character(LEN=128) filename, cmd, msg
          integer argc, IARGC, err, nprocs, rank, i, j
          integer cmode, ncid, varid0, varid1, varid2, dimid(2), nerrs
          integer(kind=MPI_OFFSET_KIND) NX, NY, len
          integer(kind=MPI_OFFSET_KIND) start(2), count(2)
          PARAMETER(NX=5, NY=2)
          integer buf(NX, NY)
          integer buftype, ghost_buftype, rec_filetype, fix_filetype
          integer array_of_sizes(2), array_of_subsizes(2)
          integer array_of_starts(2), blocklengths(2)
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size, recsize
          integer(kind=MPI_ADDRESS_KIND) a0, a1, disps(2)
          character(len = 4) :: quiet_mode
          logical verbose

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 2) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd),
     +                                ' [-q] [filename]'
              goto 999
          endif
          verbose = .TRUE.
          filename = "testfile.nc"
          call getarg(1, quiet_mode)
          if (quiet_mode(1:2) .EQ. '-q') then
              verbose = .FALSE.
              if (argc .EQ. 2) call getarg(2, filename)
          else
              if (argc .EQ. 1) call getarg(1, filename)
          endif

          nerrs = 0

          start(1) = NX * rank
          start(2) = 0
          count(1) = NX
          count(2) = 2

          ! construct various MPI derived data types
          blocklengths(1) = NX;
          blocklengths(2) = NX;
          call MPI_Get_address(buf(1,2), a0, err)
          call MPI_Get_address(buf(1,1), a1, err)
          disps(1) = 0
          disps(2) = a1 - a0
          call MPI_Type_create_hindexed(2, blocklengths, disps,
     +                  MPI_INTEGER, buftype, err)
          call MPI_Type_commit(buftype, err)

          ! create a file type for the fixed-size variable
          array_of_sizes(1) = NX*nprocs
          array_of_sizes(2) = 2
          array_of_subsizes(1) = count(1)
          array_of_subsizes(2) = count(2)
          array_of_starts(1) = start(1)
          array_of_starts(2) = start(2)
          call MPI_Type_create_subarray(2, array_of_sizes,
     +         array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN,
     +         MPI_INTEGER, fix_filetype, err)
          call MPI_Type_commit(fix_filetype, err)

          ! create a buftype with ghost cells on each side */
          array_of_sizes(1) = count(1)+4
          array_of_sizes(2) = count(2)+4
          array_of_subsizes(1) = count(1)
          array_of_subsizes(2) = count(2)
          array_of_starts(1) = 2
          array_of_starts(2) = 2
          call MPI_Type_create_subarray(2, array_of_sizes,
     +         array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN,
     +         MPI_INTEGER, ghost_buftype, err)
          call MPI_Type_commit(ghost_buftype, err)

          ! initialized buffer contents
          do j=1, count(2)
             do i=1, count(1)
                 buf(i, j) = rank*100 + j*10 + i
             enddo
          enddo

          ! create file, truncate it if exists
          cmode = NF_CLOBBER
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define 2 dimensions
          err = nfmpi_def_dim(ncid, "REC_DIM", NFMPI_UNLIMITED,dimid(2))
          call check(err, 'In nfmpi_def_dim REC_DIM: ')
          len = NX * nprocs
          err = nfmpi_def_dim(ncid, "X", len, dimid(1))
          call check(err, 'In nfmpi_def_dim X: ')

          ! define 2D record variables of integer type
          err = nfmpi_def_var(ncid, "rec_var", NF_INT, 2, dimid, varid0)
          call check(err, 'In nfmpi_def_var: rec_var ')
          err = nfmpi_def_var(ncid, "dummy_rec", NF_INT, 2,dimid,varid2)
          call check(err, 'In nfmpi_def_var: dummy_rec ')

          ! define 2D fixed-size variable of integer type
          err = nfmpi_def_dim(ncid, "FIX_DIM", 2_8, dimid(2))
          call check(err, 'In nfmpi_def_dim RECV_DIM: ')
          err = nfmpi_def_var(ncid, "fix_var", NF_INT, 2, dimid, varid1)
          call check(err, 'In nfmpi_def_var: fix_var ')

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! now we are in data mode

          ! create a file type for the record variable */
          err = nfmpi_inq_recsize(ncid, recsize)
          call check(err, 'In nfmpi_inq_recsize: ')
          blocklengths(1) = count(1)
          blocklengths(2) = count(1)
          disps(1) = start(1)*4
          disps(2) = recsize + start(1)*4
          call MPI_Type_create_hindexed(2, blocklengths, disps,
     +                                  MPI_INTEGER, rec_filetype, err)
          call MPI_Type_commit(rec_filetype, err)

          ! write the record variable */
          err = nfmpi_put_vard_all(ncid, varid0, rec_filetype, buf(:,2),
     +                             1_8, buftype)
          call check(err, 'In nfmpi_put_vard_all: ')
          ! check if the contents of buf are altered
          call check_value(rank, NX, NY, buf, nerrs)

          ! write the fixed-size variable
          err = nfmpi_put_vard_all(ncid, varid1, fix_filetype, buf(:,2),
     +                             1_8, buftype)
          call check(err, 'In nfmpi_put_vard_all: ')
          ! check if the contents of buf are altered
          call check_value(rank, NX, NY, buf, nerrs)

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

          ! open the same file and read back for validate */
          err = nfmpi_open(MPI_COMM_WORLD, filename, NF_NOWRITE,
     +                     MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_open: ')

          err = nfmpi_inq_varid(ncid, "rec_var", varid0)
          call check(err, 'In nfmpi_inq_varid rec_var: ')
          err = nfmpi_inq_varid(ncid, "fix_var", varid1)
          call check(err, 'In nfmpi_inq_varid fix_var: ')

          ! PnetCDF start() argument starts with 1
          start = start + 1
          call get_var_and_verify(ncid, varid0, NX,NY,start, count, buf,
     +                 buftype, ghost_buftype, rec_filetype, nerrs)

          call get_var_and_verify(ncid, varid1, NX,NY,start, count, buf,
     +                 buftype, ghost_buftype, fix_filetype, nerrs)

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

          call MPI_Type_free(rec_filetype, err)
          call MPI_Type_free(fix_filetype, err)
          call MPI_Type_free(buftype, err)
          call MPI_Type_free(ghost_buftype, err)

          ! check if there is any PnetCDF internal malloc residue
 998      format(A,I13,A)
          err = nfmpi_inq_malloc_size(malloc_size)
          if (err == NF_NOERR) then
              call MPI_Reduce(malloc_size, sum_size, 1, MPI_OFFSET,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, err)
              if (rank .EQ. 0 .AND. sum_size .GT. 0_8) print 998,
     +            'heap memory allocated by PnetCDF internally has ',
     +            sum_size/1048576, ' MiB yet to be freed'
          endif

          msg = '*** TESTING F90 '//trim(cmd)//' for vard API '
          if (rank .eq. 0) then
              if (nerrs .EQ. 0) then
                  write(*,"(A67,A)") msg,'------ pass'
              else
                  write(*,"(A67,A)") msg,'------ failed'
              endif
          endif

 999      call MPI_Finalize(err)
      end program main

