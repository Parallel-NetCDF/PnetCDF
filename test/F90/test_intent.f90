!
!   Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This program tests Fortran 90 INTENT modifier used in PnetCDF.
! It also tests put_att on Little Endian when using parameter
! buffer (read-only memory)
!
      subroutine check(err, message)
          use mpi
          use pnetcdf
          implicit none
          integer err
          character(len=*) message
          character(len=256) msg

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF90_NOERR) then
              write(6,*) trim(message), trim(nf90mpi_strerror(err))
              msg = '*** TESTING F90 test_intent.f90 '
              call pass_fail(1, msg, 0)
              STOP 2
          end if
      end subroutine check

      program main
          use mpi
          use pnetcdf
          implicit none

          integer, parameter ::   OneByteInt = selected_int_kind(2), &
                                  TwoByteInt = selected_int_kind(4), &
                                 FourByteInt = selected_int_kind(9), &
                                EightByteInt = selected_int_kind(18)

          character(LEN=256) out_path, in_path, cmd, msg
          integer err, ierr, rank, get_args
          integer cmode, ncid, varid, dimid(1), req(1), status(1)
          integer(kind=MPI_OFFSET_KIND) start(1)
          integer(kind=MPI_OFFSET_KIND) count(1)
          integer(kind=MPI_OFFSET_KIND) bufsize
          logical keep_files

          character(LEN=3)           cbuf
          integer(KIND=OneByteInt)   i1buf(3)
          integer(KIND=TwoByteInt)   sbuf(3)
          integer(KIND=FourByteInt)  ibuf(3)
          real                       rbuf(3)
          double precision           dbuf(3)
          integer(KIND=EightByteInt) i8buf(3)

          PARAMETER( cbuf="123")
          PARAMETER(i1buf=(/1_OneByteInt,2_OneByteInt,3_OneByteInt/))
          PARAMETER( sbuf=(/1_TwoByteInt,2_TwoByteInt,3_TwoByteInt/))
          PARAMETER( ibuf=(/1,2,3/))
          PARAMETER( rbuf=(/1.0,2.0,3.0/))
          PARAMETER( dbuf=(/1.0,2.0,3.0/))
          PARAMETER(i8buf=(/1_EightByteInt,2_EightByteInt,3_EightByteInt/))
          double precision timing

          call MPI_Init(ierr)

          timing = MPI_Wtime()

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

          ! take out_path from command-line argument if there is any
          cmd = ' '
          if (rank .EQ. 0) then
              out_path = 'testfile.nc'
              err = get_args(cmd, out_path, in_path, keep_files)
          endif
          call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          if (err .EQ. 0) goto 999

          call MPI_Bcast(out_path, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

          call MPI_Bcast(keep_files, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

          ! create file, truncate it if exists
          cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
          err = nf90mpi_create(MPI_COMM_WORLD, out_path, cmode, &
                               MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_text', cbuf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_int1', i1buf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_int2', sbuf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_int', ibuf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_real', rbuf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_double', dbuf)
          call check(err, 'In nf90mpi_put_att: ')
          err = nf90mpi_put_att(ncid, NF90_GLOBAL, 'nf90_attr_int8', i8buf)
          call check(err, 'In nf90mpi_put_att: ')

          err = nfmpi_put_att_text(ncid, NF90_GLOBAL, 'nf_attr_text', 3_MPI_OFFSET_KIND, cbuf)
          call check(err, 'In nfmpi_put_att_text: ')
          err = nfmpi_put_att_int1(ncid, NF90_GLOBAL, 'nf_attr_int1', NF90_INT1, 3_MPI_OFFSET_KIND, i1buf)
          call check(err, 'In nfmpi_put_att_int1: ')
          err = nfmpi_put_att_int2(ncid, NF90_GLOBAL, 'nf_attr_int2', NF90_INT2, 3_MPI_OFFSET_KIND, sbuf)
          call check(err, 'In nfmpi_put_att_int2: ')
          err = nfmpi_put_att_int(ncid, NF90_GLOBAL, 'nf_attr_int', NF90_INT, 3_MPI_OFFSET_KIND, ibuf)
          call check(err, 'In nfmpi_put_att_int: ')
          err = nfmpi_put_att_real(ncid, NF90_GLOBAL, 'nf_attr_real', NF90_FLOAT, 3_MPI_OFFSET_KIND, rbuf)
          call check(err, 'In nfmpi_put_att_real: ')
          err = nfmpi_put_att_double(ncid, NF90_GLOBAL, 'nf_attr_double', NF90_DOUBLE, 3_MPI_OFFSET_KIND, dbuf)
          call check(err, 'In nfmpi_put_att_double: ')
          err = nfmpi_put_att_int8(ncid, NF90_GLOBAL, 'nf_attr_int8', NF90_INT64, 3_MPI_OFFSET_KIND, i8buf)
          call check(err, 'In nfmpi_put_att_int8: ')

          ! define a variable of an integer array of size 3 in the nc file
          err = nf90mpi_def_dim(ncid, 'X', 3_MPI_OFFSET_KIND, dimid(1))
          call check(err, 'In nf90mpi_def_dim: ')

          err = nf90mpi_def_var(ncid, 'var', NF90_INT, dimid, varid)
          call check(err, 'In nf90mpi_def_var: ')

          ! fill with default fill value
          err = nf90mpi_def_var_fill(ncid, varid, 0, NF90_FILL_INT)
          call check(err, 'In nf90mpi_def_var_fill: ')

          err = nf90mpi_enddef(ncid)
          call check(err, 'In nf90mpi_enddef: ')

          ! bufsize must be max of data type converted before and after
          bufsize = 3*4
          err = nf90mpi_buffer_attach(ncid, bufsize)
          call check(err, 'In nf90mpi_buffer_attach: ')

          start(1) = 1
          count(1) = 3
          err = nfmpi_bput_vara_int(ncid, varid, start, count, ibuf(1:), req(1))
          call check(err, 'In nfmpi_bput_vara_int: ')

          err = nf90mpi_wait_all(ncid, 1, req, status)
          call check(err, 'In nf90mpi_wait_all: ')

          if (status(1) .ne. NF90_NOERR) then
              print*,'Error at bput status ', nf90mpi_strerror(status(1))
          endif

          err = nf90mpi_buffer_detach(ncid)
          call check(err, 'In nf90mpi_buffer_detach: ')

          ! close the file
          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')

 999      timing = MPI_Wtime() - timing
          call MPI_Allreduce(MPI_IN_PLACE, timing, 1, &
                             MPI_DOUBLE_PRECISION, MPI_MAX, &
                             MPI_COMM_WORLD, ierr)

          if (rank .eq. 0) then
              if (.NOT. keep_files) then
                  err = nf90mpi_delete(out_path, MPI_INFO_NULL)
              end if

              msg = '*** TESTING F90 '//trim(cmd)//' - INTENT modifier '
              call pass_fail(0, msg, timing)
          end if

          call MPI_Finalize(ierr)
      end program main

