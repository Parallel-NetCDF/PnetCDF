!
!   Copyright (C) 2022, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! This program tests if PnetCDF handles Fortran predefine datatypes when
! bufcount argument in the flexible APIs is -1 and buftype argument is a
! predefined MPI datatype
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
          integer err, XTRIM
          character*(*) message
          character*128 msg

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) message(1:XTRIM(message)), nfmpi_strerror(err)
              msg = '*** TESTING F77 flexible_api.f for flexible API '
              call pass_fail(1, msg)
              STOP 2
          end if
      end ! subroutine check

      program main
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          character*256 filename, cmd, msg
          integer XTRIM
          integer err, ierr, nerrs, nprocs, rank, i, j
          integer cmode, ncid, varid, dimid(2), ghost_len, get_args
          integer*8 nx, ny, global_nx, global_ny
          integer*8 starts(2), counts(2), bufcount, one
          PARAMETER(nx=5, ny=4, ghost_len=3)
          integer buf(nx+2*ghost_len, ny+2*ghost_len)
          integer subarray, req(1), st(1)
          integer array_of_sizes(2), array_of_subsizes(2)
          integer array_of_starts(2)
          integer*8 malloc_size, sum_size
          logical verbose
          integer dummy

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              verbose = .TRUE.
              filename = "testfile.nc"
              ierr = get_args(cmd, filename)
          endif
          call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
          if (ierr .EQ. 0) goto 999

          call MPI_Bcast(verbose, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     +                   err)
          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0,
     +                   MPI_COMM_WORLD, err)

          nerrs = 0

          ! set parameters
          global_nx = nx
          global_ny = ny * nprocs

          ! first initialize the entire buffer to -1
          do j=1, ny+2*ghost_len
             do i=1, nx+2*ghost_len
                buf(i, j) = -1
             enddo
          enddo
          ! assign values for non-ghost cells
          do j=ghost_len+1, ny+ghost_len
             do i=ghost_len+1, nx+ghost_len
                 buf(i, j) = rank
             enddo
          enddo

          ! define an MPI datatype using MPI_Type_create_subarray()
          array_of_sizes(1)    = nx + 2*ghost_len
          array_of_sizes(2)    = ny + 2*ghost_len
          array_of_subsizes(1) = nx
          array_of_subsizes(2) = ny
          array_of_starts(1)   = ghost_len  ! MPI start index starts with 0
          array_of_starts(2)   = ghost_len
          call MPI_Type_create_subarray(2, array_of_sizes,
     +         array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN,
     +         MPI_INTEGER, subarray, err)
          call MPI_Type_commit(subarray, err)

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define dimensions x and y
          err = nfmpi_def_dim(ncid, "y", global_ny, dimid(2))
          call check(err, 'In nfmpi_def_dim y: ')
          err = nfmpi_def_dim(ncid, "x", global_nx, dimid(1))
          call check(err, 'In nfmpi_def_dim x: ')

          ! define a 2D variable of integer type
          err = nfmpi_def_var(ncid, "var", NF_INT, 2, dimid, varid)
          call check(err, 'In nfmpi_def_var: ')

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! now we are in data mode

          ! Note that in Fortran, array indices start with 1
          starts(1) = 1
          starts(2) = ny * rank + 1
          counts(1) = nx
          counts(2) = ny
          bufcount  = 1

          ! In Fortran, subarray buffer type can also use array index ranges
          ! for example in this case
          !    buf(1+ghost_len:nx+ghost_len, 1+ghost_len:ny+ghost_len)
          ! However, this does not work for nonblocking APIs because
          ! wait/wait_all is called separately from the nonblocking calls
          ! Fortran the subarray indexing is lost in the wait call
          err = nfmpi_put_vara_all(ncid, varid, starts, counts, buf,
     +                             bufcount, subarray)
          call check(err, 'In nfmpi_put_vara_all: ')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          ! test blocking and nonblocking APIs when buftype argument is
          ! a predefined MPI datatype
          bufcount = counts(1) * counts(2)
          one = 1

          err = nfmpi_put_vara_all(ncid, varid, starts, counts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_put_vara_all: ')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_put_var1_all(ncid, varid, starts, buf,
     +                             one, MPI_INTEGER)
          call check(err, 'In nfmpi_put_var1_all: ')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_iput_vara(ncid, varid, starts, counts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iput_vara: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_iput_var1(ncid, varid, starts, buf,
     +                          one, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iput_var1: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_get_vara_all(ncid, varid, starts, counts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_get_vara_all: ')

          err = nfmpi_get_var1_all(ncid, varid, starts, buf,
     +                             one, MPI_INTEGER)
          call check(err, 'In nfmpi_get_var1_all: ')

          err = nfmpi_iget_vara(ncid, varid, starts, counts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iget_vara: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          err = nfmpi_iget_var1(ncid, varid, starts, buf,
     +                          one, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iget_var1: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! test blocking and nonblocking APIs when bufcount argument is
          ! -1 and buftype argument is a predefined MPI datatype
          bufcount = -1

          err = nfmpi_put_vara_all(ncid, varid, starts, counts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_put_vara_all: ')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_put_var1_all(ncid, varid, starts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_put_var1_all: ')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_iput_vara(ncid, varid, starts, counts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iput_vara: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_iput_var1(ncid, varid, starts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iput_var1: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! flush in case burst buffering is enabled
          err = nfmpi_flush(ncid)
          call check(err, 'In nfmpi_flush: ')

          err = nfmpi_get_vara_all(ncid, varid, starts, counts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_get_vara_all: ')

          err = nfmpi_get_var1_all(ncid, varid, starts, buf,
     +                             bufcount, MPI_INTEGER)
          call check(err, 'In nfmpi_get_var1_all: ')

          err = nfmpi_iget_vara(ncid, varid, starts, counts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iget_vara: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          err = nfmpi_iget_var1(ncid, varid, starts, buf,
     +                          bufcount, MPI_INTEGER, req)
          call check(err, 'In nfmpi_iget_var1: ')

          err = nfmpi_wait_all(ncid, 1, req, st)
          call check(err, 'In nfmpi_wait_all:')

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

          call MPI_Type_free(subarray, err)

          ! check if there is any PnetCDF internal malloc residue
 998      format(A,I13,A)
          err = nfmpi_inq_malloc_size(malloc_size)
          if (err .EQ. NF_NOERR) then
              call MPI_Reduce(malloc_size, sum_size, 1, MPI_INTEGER8,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, err)
              if (rank .EQ. 0 .AND. sum_size .GT. 0)
     +            print 998,
     +            'heap memory allocated by PnetCDF internally has ',
     +            sum_size, ' bytes yet to be freed'
          endif

          if (rank .eq. 0) then
              msg = '*** TESTING F77 '//cmd(1:XTRIM(cmd))//
     +              ' for bufcount=-1 & buftype predefined'
              call pass_fail(nerrs, msg)
          endif

 999      call MPI_Finalize(ierr)
          if (nerrs .GT. 0) stop 2

      end ! program main

