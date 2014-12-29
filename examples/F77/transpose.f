!
!   Copyright (C) 2014, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! $Id$

!
!    This example shows how to use varm API to write six 3D integer
!    array variables into a file. Each variable in the file is a
!    dimensional transposed array from the one stored in memory. In
!    memory, a 3D array is partitioned among all processes in a
!    block-block-block fashion and in XYZ (i.e. Fortran) order. Note
!    the variable and dimension naming below is in Fortran order. The
!    dimension structures of the transposed six arrays are
!       integer XYZ_var(X, Y, Z)      XYZ -> XYZ (no transpose)
!       integer XZY_var(X, Z, Y)      XYZ -> XZY
!       integer YXZ_var(Y, X, Z)      XYZ -> YXZ
!       integer YZX_var(Y, Z, X)      XYZ -> YZX
!       integer ZXY_var(Z, X, Y)      XYZ -> ZXY
!       integer ZYX_var(Z, Y, X)      XYZ -> ZYX
!
!    To compile:
!        mpif77 -O2 transpose.f -o transpose -lpnetcdf
!    To run:
!        mpiexec -n num_processes ./transpose [filename] [len]
!    where len decides the size of local array, which is
!    (len+2) x (len+1) x len. So, each variable is of size
!    (len+2)*(len+1)*len * nprocs * sizeof(int)
!
      subroutine check(err, message)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'
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
          include 'mpif.h'
          include 'pnetcdf.inc'

          character(LEN=128) filename, cmd
          integer argc, IARGC, err, nprocs, rank, lower_dims
          integer cmode, ncid, psizes(3), dimids(3), dimidsT(3)
          integer XYZ_id, XZY_id, YZX_id, YXZ_id, ZYX_id, ZXY_id
          integer(kind=MPI_OFFSET_KIND) gsizes(3), starts(3)
          integer(kind=MPI_OFFSET_KIND) counts(3), strides(3), imap(3)
          integer(kind=MPI_OFFSET_KIND) startsT(3), countsT(3)
          integer(kind=MPI_OFFSET_KIND) malloc_size, sum_size
          integer i, j, k, nx, ny, nz
          PARAMETER(nx=4, ny=3, nz=2)
          integer buf(nx,ny,nz)
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
     +                         ' [-q] [filename]'
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

          ! calculate number of processes along each dimension
          psizes = 0
          call MPI_Dims_create(nprocs, 3, psizes, err)
          if (verbose .AND. rank .EQ. 0) print*, "psizes= ",psizes(1:3)

          ! for each MPI rank, find its local rank IDs along each dimension in
          ! starts()
          lower_dims = 1
          do i=1, 3
              starts(i) = MOD(rank / lower_dims, psizes(i))
              lower_dims = lower_dims * psizes(i)
          enddo
          if (verbose)
     +        print*, "proc ",rank,": dim rank= ",starts(1:3)

          strides = 1
          gsizes  = (/ nx, ny, nz /)
          do i=1, 3
             starts(i) = starts(i) * gsizes(i) + 1 ! start indices
             counts(i) = gsizes(i)                 ! array elements
             gsizes(i) = gsizes(i) * psizes(i)     ! global array size
          enddo

          if (verbose) then
              print*, "proc ",rank,": starts= ",starts(1:3)
              print*, "proc ",rank,": counts= ",counts(1:3)
          endif

          ! initialize buffer with contiguous numbers
          do k=1, nz
          do j=1, ny
          do i=1, nx
              buf(i, j, k) = INT((starts(3)+k-2)*gsizes(2)*gsizes(1) +
     +                           (starts(2)+j-2)*gsizes(1) +
     +                           (starts(1)+i-2))
          enddo
          enddo
          enddo
          if (verbose .AND. rank .EQ. 0) print*, "buf= ",buf

          ! create file, truncate it if exists
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
          cmode = NF_CLOBBER
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! define dimensions X, Y, Z
          err = nfmpi_def_dim(ncid, "X", gsizes(1), dimids(1))
          call check(err, 'In nfmpi_def_dim X: ')
          err = nfmpi_def_dim(ncid, "Y", gsizes(2), dimids(2))
          call check(err, 'In nfmpi_def_dim Y: ')
          err = nfmpi_def_dim(ncid, "Z", gsizes(3), dimids(3))
          call check(err, 'In nfmpi_def_dim Z: ')

          ! define variable with no transposed file layout: XYZ
          err = nfmpi_def_var(ncid, "XYZ_var", NF_INT, 3, dimids,
     +                        XYZ_id)
          call check(err, 'In nfmpi_def_var XYZ_var: ')

          ! define variable with transposed file layout: XYZ -> XZY
          dimidsT = (/ dimids(1), dimids(3), dimids(2) /)
          err = nfmpi_def_var(ncid, "XZY_var", NF_INT, 3, dimidsT,
     +                        XZY_id)
          call check(err, 'In nfmpi_def_var XZY_var: ')

          ! define variable with transposed file layout: XYZ -> YXZ
          dimidsT = (/ dimids(2), dimids(1), dimids(3) /)
          err = nfmpi_def_var(ncid, "YXZ_var", NF_INT, 3, dimidsT,
     +                        YXZ_id)
          call check(err, 'In nfmpi_def_var YXZ_var: ')

          ! define variable with transposed file layout: XYZ -> YZX
          dimidsT = (/ dimids(2), dimids(3), dimids(1) /)
          err = nfmpi_def_var(ncid, "YZX_var", NF_INT, 3, dimidsT,
     +                        YZX_id)
          call check(err, 'In nfmpi_def_var YZX_var: ')

          ! define variable with transposed file layout: XYZ -> ZXY
          dimidsT = (/ dimids(3), dimids(1), dimids(2) /)
          err = nfmpi_def_var(ncid, "ZXY_var", NF_INT, 3, dimidsT,
     +                        ZXY_id)
          call check(err, 'In nfmpi_def_var ZXY_var: ')

          ! define variable with transposed file layout: XYZ -> ZYX
          dimidsT = (/ dimids(3), dimids(2), dimids(1) /)
          err = nfmpi_def_var(ncid, "ZYX_var", NF_INT, 3, dimidsT,
     +                        ZYX_id)
          call check(err, 'In nfmpi_def_var ZYX_var: ')

          ! do not forget to exit define mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! now we are in data mode

          ! write the whole variable in file: XYZ
          err = nfmpi_put_vara_int_all(ncid, XYZ_id, starts, counts,
     +                                 buf(:,:,:))
          call check(err, 'In nfmpi_put_vara_int_all XYZ: ')

          ! write the transposed variable:  XYZ -> XZY
          imap    = (/ 1_8, counts(1)*counts(2), counts(1) /)
          startsT = (/ starts(1), starts(3), starts(2) /)
          countsT = (/ counts(1), counts(3), counts(2) /)
          err = nfmpi_put_varm_int_all(ncid, XZY_id, startsT, countsT,
     +                                 strides, imap, buf)
          call check(err, 'In nfmpi_put_varm_int_all XZY: ')

          ! write the transposed variable:  XYZ -> YXZ
          imap    = (/ counts(1), 1_8, counts(1)*counts(2) /)
          startsT = (/ starts(2), starts(1), starts(3) /)
          countsT = (/ counts(2), counts(1), counts(3) /)
          err = nfmpi_put_varm_int_all(ncid, YXZ_id, startsT, countsT,
     +                                 strides, imap, buf)
          call check(err, 'In nfmpi_put_varm_int_all YZX: ')

          ! write the transposed variable:  XYZ -> YZX
          imap    = (/ counts(1), counts(1)*counts(2), 1_8 /)
          startsT = (/ starts(2), starts(3), starts(1) /)
          countsT = (/ counts(2), counts(3), counts(1) /)
          err = nfmpi_put_varm_int_all(ncid, YZX_id, startsT, countsT,
     +                                 strides, imap, buf)
          call check(err, 'In nfmpi_put_varm_int_all YZX: ')

          ! write the transposed variable:  XYZ -> ZXY
          imap    = (/ counts(1)*counts(2), 1_8, counts(1) /)
          startsT = (/ starts(3), starts(1), starts(2) /)
          countsT = (/ counts(3), counts(1), counts(2) /)
          err = nfmpi_put_varm_int_all(ncid, ZXY_id, startsT, countsT,
     +                                 strides, imap, buf)
          call check(err, 'In nfmpi_put_varm_int_all ZXY: ')

          ! write the transposed variable:  XYZ -> ZYX
          imap    = (/ counts(1)*counts(2), counts(1), 1_8 /)
          startsT = (/ starts(3), starts(2), starts(1) /)
          countsT = (/ counts(3), counts(2), counts(1) /)
          err = nfmpi_put_varm_int_all(ncid, ZYX_id, startsT, countsT,
     +                                 strides, imap, buf)
          call check(err, 'In nfmpi_put_varm_int_all ZYX: ')

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

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

