!
!  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!    This example is the Fortran 90 version of nonblocking_write.c
!    It creates a netcdf file in CD-5 format and writes a number of
!    3D integer non-record vaiables.
!    Usage: (for example)
!    To compile:
!        mpif90 -O2 nonblocking_write.f90 -o nonblocking_write -lpnetcdf
!    To run:
!        mpiexec -n num_processes ./nonblocking_write [filename] [len]
!    where len decides the size of each local array, which is
!    len x len x len. Each non-record variable is of size
!          len*len*len * nprocs * sizeof(int)
!    All variables are partitioned among all processes in a 3D
!    block-block-block fashion.
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

          integer NDIMS, NUM_VARS
          PARAMETER(NDIMS=3, NUM_VARS=10)

          character(LEN=128) filename, cmd, str
          integer i, cmode, argc, iargc, err
          integer rank, nprocs, len, bufsize, ncid
          integer psizes(NDIMS), dimids(NDIMS), varids(NUM_VARS)
          integer req(NUM_VARS), st(NUM_VARS)
          integer(kind=MPI_OFFSET_KIND) gsizes(NDIMS)
          integer(kind=MPI_OFFSET_KIND) starts(NDIMS), counts(NDIMS)
          integer(kind=MPI_OFFSET_KIND) bbufsize

          ! define an array of pointers
          type intp
              integer, dimension(:), pointer :: p
          end type intp
          type(intp) buf(NUM_VARS)

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 2) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd), &
                                      ' [filename] [len]'
              goto 999
          endif
          filename = "testfile.nc"
          if (argc .GE. 1) call getarg(1, filename)
          len = 10
          if (argc .EQ. 2) then
             call getarg(2, str)
             read (str,'(I10)') len
          endif

          do i=1,NDIMS
             psizes(i) = 0
          enddo

          ! create a block-block-block data partitioning
          call MPI_Dims_create(nprocs, NDIMS, psizes, err);
          starts(1) = mod(rank, psizes(1))
          starts(2) = mod((rank / psizes(2)), psizes(2))
          starts(3) = mod((rank / (psizes(1) * psizes(2))), psizes(3))

          bufsize = 1
          do i=1,NDIMS
              gsizes(i) = len * psizes(i)
              starts(i) = starts(i) * len + 1;
              counts(i) = len;
              bufsize   = bufsize * len;
          enddo

          ! allocate buffer and initialize with rank IDs
          do i=1, NUM_VARS
             allocate(buf(i)%p(bufsize))
             buf(i)%p(:) = rank
          enddo

          ! create file, truncate it if exists
          cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
          err = nf90mpi_create(MPI_COMM_WORLD, filename, cmode, &
                               MPI_INFO_NULL, ncid)
          call check(err, 'In nf90mpi_create: ')

          ! define dimensions
          do i=1, NDIMS
             write(str,'(I2)') i
             err = nf90mpi_def_dim(ncid, "x"//trim(ADJUSTL(str)), &
                                   gsizes(i), dimids(i))
             call check(err, 'In nf90mpi_def_dim x'//trim(str))
          enddo

          !define variables
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nf90mpi_def_var(ncid, "var"//trim(ADJUSTL(str)), &
                                   NF90_INT, dimids, varids(i))
             call check(err, 'In nf90mpi_def_var var'//trim(str))
          enddo

          ! do not forget to exit define mode
          err = nf90mpi_enddef(ncid)
          call check(err, 'In nf90mpi_enddef: ')

          ! write one variable at a time using iput
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nf90mpi_iput_var(ncid, varids(i), buf(i)%p, req(i), &
                                    starts, counts)
             call check(err, 'In nf90mpi_iput_var '//trim(str))
          enddo

          ! wait for the nonblocking I/O to complete
          err = nf90mpi_wait_all(ncid, NUM_VARS, req, st)
          call check(err, 'In nf90mpi_wait_all')

          ! check the status of each nonblocking request
          do i=1, NUM_VARS
             write(str,'(I2)') i
             call check(st(i), 'In nf90mpi_wait_all req '//trim(str))
          enddo

          ! write one variable at a time using bput

          ! bbufsize must be max of data type converted before and after
          bbufsize = bufsize * NUM_VARS * 4  ! 4 is size of integer
          err = nf90mpi_buffer_attach(ncid, bbufsize)
          call check(err, 'In nf90mpi_buffer_attach')

          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nf90mpi_bput_var(ncid, varids(i), buf(i)%p, req(i), &
                                    starts, counts)
             call check(err, 'In nf90mpi_bput_var '//trim(str))

             ! can safely free up the buffer here
             deallocate(buf(i)%p)
          enddo

          ! wait for the nonblocking I/O to complete
          err = nf90mpi_wait_all(ncid, NUM_VARS, req, st)
          call check(err, 'In nf90mpi_wait_all')

          ! check the status of each nonblocking request
          do i=1, NUM_VARS
             write(str,'(I2)') i
             call check(st(i), 'In nf90mpi_wait_all req '//trim(str))
          enddo

          ! detach the temporary buffer
          err = nf90mpi_buffer_detach(ncid)
          call check(err, 'In nf90mpi_buffer_detach')

          ! close the file
          err = nf90mpi_close(ncid);
          call check(err, 'In nf90mpi_close')

 999      call MPI_Finalize(err)

      end program main

