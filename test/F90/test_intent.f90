!
!   Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
! This program tests put_att on Little Endian when using parameter
! buffer (read-only memory)
!
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

          character(LEN=128) filename, cmd
          integer argc, IARGC, err, rank
          integer cmode, ncid

          character(LEN=3)  cbuf
          integer*1        i1buf(3)
          integer*2         sbuf(3)
          integer           ibuf(3)
          real              rbuf(3)
          double precision  dbuf(3)
          integer*8        i8buf(3)

          PARAMETER( cbuf="123")
          PARAMETER(i1buf=(/1,2,3/))
          PARAMETER( sbuf=(/1,2,3/))
          PARAMETER( ibuf=(/1,2,3/))
          PARAMETER( rbuf=(/1.0,2.0,3.0/))
          PARAMETER( dbuf=(/1.0,2.0,3.0/))
          PARAMETER(i8buf=(/1_8,2_8,3_8/))

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 1) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd),' [filename]'
              goto 999
          endif
          filename = "testfile.nc"
          if (argc .EQ. 1) call getarg(1, filename)

          ! create file, truncate it if exists
          cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
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

          err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'nf_attr_text', 3_8, cbuf)
          call check(err, 'In nfmpi_put_att_text: ')
          err = nfmpi_put_att_int1(ncid, NF_GLOBAL, 'nf_attr_int1', NF_INT1, 3_8, i1buf)
          call check(err, 'In nfmpi_put_att_int1: ')
          err = nfmpi_put_att_int2(ncid, NF_GLOBAL, 'nf_attr_int2', NF_INT2, 3_8, sbuf)
          call check(err, 'In nfmpi_put_att_int2: ')
          err = nfmpi_put_att_int(ncid, NF_GLOBAL, 'nf_attr_int', NF_INT, 3_8, ibuf)
          call check(err, 'In nfmpi_put_att_int: ')
          err = nfmpi_put_att_real(ncid, NF_GLOBAL, 'nf_attr_real', NF_FLOAT, 3_8, rbuf)
          call check(err, 'In nfmpi_put_att_real: ')
          err = nfmpi_put_att_double(ncid, NF_GLOBAL, 'nf_attr_double', NF_DOUBLE, 3_8, dbuf)
          call check(err, 'In nfmpi_put_att_double: ')
          err = nfmpi_put_att_int8(ncid, NF_GLOBAL, 'nf_attr_int8', NF_INT64, 3_8, i8buf)
          call check(err, 'In nfmpi_put_att_int8: ')

          ! close the file
          err = nf90mpi_close(ncid)
          call check(err, 'In nf90mpi_close: ')

 999      call MPI_Finalize(err)
      end program main

