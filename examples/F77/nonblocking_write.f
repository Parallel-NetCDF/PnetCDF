!
!  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!    This example is the Fortran 77 version of nonblocking_write.c
!    It creates a netcdf file in CD-5 format and writes a number of
!    3D integer non-record vaiables.
!    Usage: (for example)
!    To compile:
!        mpif77 -O2 nonblocking_write.f -o nonblocking_write -lpnetcdf
!    To run (for example):
!        mpiexec -n 32 ./nonblocking_write /pvfs2/wkliao/testfile.nc 10
!    The size of each local array is len x len x len. Each non-record
!    variable is of size len*len*len * nprocs * sizeof(int)
!    All variables are partitioned among all processes in a 3D
!    block-block-block fashion.
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

      program main
          implicit none
          include "mpif.h"
          include "pnetcdf.inc"

          integer NDIMS, NUM_VARS, BUFSIZE
          PARAMETER(NDIMS=3, NUM_VARS=10, BUFSIZE=1000)

          character(LEN=128) filename, cmd, str
          integer i, cmode, argc, iargc, err
          integer rank, nprocs, ncid, len
          integer buf(BUFSIZE, NUM_VARS)
          integer psizes(NDIMS), dimids(NDIMS), varids(NUM_VARS)
          integer req(NUM_VARS), st(NUM_VARS)
          integer(kind=MPI_OFFSET_KIND) gsizes(NDIMS)
          integer(kind=MPI_OFFSET_KIND) starts(NDIMS), counts(NDIMS)
          integer(kind=MPI_OFFSET_KIND) bbufsize
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size
          character(len = 4) :: quiet_mode
          logical verbose

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! take filename from command-line argument if there is any
          call getarg(0, cmd)
          argc = IARGC()
          if (argc .GT. 3) then
              if (rank .EQ. 0) print*,'Usage: ',trim(cmd),
     +                         ' [-q] [filename] [len]'
              goto 999
          endif
          verbose = .TRUE.
          filename = "testfile.nc"
          len = 10
          call getarg(1, quiet_mode)
          if (quiet_mode(1:2) .EQ. '-q') then
              verbose = .FALSE.
              if (argc .GE. 2) call getarg(2, filename)
              if (argc .EQ. 3) then
                  call getarg(3, str)
                  read (str,'(I10)') len
              endif
          else
              if (argc .GE. 1) call getarg(1, filename)
              if (argc .EQ. 2) then
                  call getarg(2, str)
                  read (str,'(I10)') len
              endif
          endif

          if (len .GT. BUFSIZE) then
             print*,'Maximum len is 10, change it to 10'
             len = BUFSIZE
          endif

          do i=1,NDIMS
             psizes(i) = 0
          enddo

          ! create a block-block-block data partitioning
          call MPI_Dims_create(nprocs, NDIMS, psizes, err);
          starts(1) = mod(rank, psizes(1))
          starts(2) = mod((rank / psizes(2)), psizes(2))
          starts(3) = mod((rank / (psizes(1) * psizes(2))), psizes(3))

          bbufsize = 1
          do i=1,NDIMS
              gsizes(i) = len * psizes(i)
              starts(i) = starts(i) * len + 1
              counts(i) = len
              bbufsize  = bbufsize * len
          enddo

          buf = rank

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define dimensions
          do i=1, NDIMS
             write(str,'(I2)') i
             err = nfmpi_def_dim(ncid, "x"//trim(ADJUSTL(str)),
     +                           gsizes(i), dimids(i))
             call check(err, 'In nfmpi_def_dim x'//trim(str))
          enddo

          !define variables
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nfmpi_def_var(ncid, "var"//trim(ADJUSTL(str)),
     +                           NF_INT, NDIMS, dimids, varids(i))
             call check(err, 'In nfmpi_def_var var'//trim(str))
          enddo

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! write one variable at a time using iput
          do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nfmpi_iput_vara_int(ncid, varids(i), starts,
     +                                 counts, buf(:,i), req(i))
             call check(err, 'In nfmpi_iput_vara_int '//trim(str))
          enddo

          ! wait for the nonblocking I/O to complete
          err = nfmpi_wait_all(ncid, NUM_VARS, req, st)
          call check(err, 'In nfmpi_wait_all')

          ! check the status of each nonblocking request
          do i=1, NUM_VARS
             write(str,'(I2)') i
             call check(st(i), 'In nfmpi_wait_all req '//trim(str))
          enddo

          ! write one variable at a time using bput

          ! bbufsize must be max of data type converted before and after
          bbufsize = bbufsize * NUM_VARS * 4  ! 4 is size of integer
          err = nfmpi_buffer_attach(ncid, bbufsize)
          call check(err, 'In nfmpi_buffer_attach')

           do i=1, NUM_VARS
             write(str,'(I2)') i
             err = nfmpi_bput_vara_int(ncid, varids(i), starts,
     +                                 counts, buf(:,i), req(i))
             call check(err, 'In nfmpi_bput_vara_int '//trim(str))

             ! can safely change the contents of buf(:,i) now
           enddo

          ! wait for the nonblocking I/O to complete
          err = nfmpi_wait_all(ncid, NUM_VARS, req, st)
          call check(err, 'In nfmpi_wait_all')

          ! check the status of each nonblocking request
          do i=1, NUM_VARS
             write(str,'(I2)') i
             call check(st(i), 'In nfmpi_wait_all req '//trim(str))
          enddo

          ! detach the temporary buffer
          err = nfmpi_buffer_detach(ncid)
          call check(err, 'In nfmpi_buffer_detach')

          ! close the file
          err = nfmpi_close(ncid);
          call check(err, 'In nfmpi_close')

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

 999      call MPI_Finalize(err)
          call EXIT(0)

      end program main

