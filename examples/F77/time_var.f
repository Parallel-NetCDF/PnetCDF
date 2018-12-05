!
!   Copyright (C) 2018, Northwestern University
!   See COPYRIGHT notice in top-level directory.
!
! This example shows how to create and write record and fixed-size variables.
! It first defines a record 2D variable of size time * global_nx where
!    time is a expandable dimension and
!    global_nx == (NX * number of MPI processes).
! The data partitioning pattern is a column-wise partitioning across all
! processes. Each process writes a subarray of size 1 * nx.
! It then defines a fixed-size 2D variable of size global_ny * global_nx where
!    global_ny == NY and
!    global_nx == (NX * number of MPI processes).
! The data partitioning pattern is a column-wise partitioning across all
! processes. Each process writes a subarray of size ny * nx.
!
!    To compile:
!        mpif77 -O2 time_var.f -o time_var -lpnetcdf
!
! Example commands for MPI run and outputs from running ncmpidump on the
! output NetCDF file produced by this example program:
!
!    % mpiexec -n 4 ./time_var /pvfs2/wkliao/testfile.nc
!
!    % ncmpidump /pvfs2/wkliao/testfile.nc
!    netcdf testfile {
!    // file format: CDF-5
!    dimensions:
!            time = UNLIMITED ; // (2 currently)
!            Y = 4 ;
!            X = 12 ;
!    variables:
!            float rec_var(time, X) ;
!            float fix_var(Y, X) ;
!
!    // global attributes:
!                    :history = "Mon Aug 13 21:27:48 2018" ;
!    data:
!
!     rec_var =
!      100, 100, 100, 101, 101, 101, 102, 102, 102, 103, 103, 103,
!      200, 200, 200, 201, 201, 201, 202, 202, 202, 203, 203, 203 ;
!
!     fix_var =
!      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
!      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
!      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
!      0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3 ;
!    }
!

      subroutine check(err, message)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'
          integer err
          character message*(*)

          ! It is a good idea to check returned value for possible error
          if (err .NE. NF_NOERR) then
              write(6,*) message//' '//nfmpi_strerror(err)
              call MPI_Abort(MPI_COMM_WORLD, -1, err)
          end if
      end ! subroutine check

      subroutine pnetcdf_write(filename, mode)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'

          character*256 filename
          integer mode

          ! local variables
          integer i, j, err, nprocs, rank, cmode
          integer ncid, rec_var, fix_var, dimid(2), dim_t, dim_y, dim_x
          integer*8 nx, ny, global_nx, global_ny, attr_len
          integer*8 start(2), count(2)
          PARAMETER(nx=3, ny=4)
          double precision buf(nx,ny)

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! create file, clobber it if exists
          cmode = IOR(NF_CLOBBER, mode)
          err = nfmpi_create(MPI_COMM_WORLD, filename, cmode,
     +                       MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_create: ')

          ! set set the global dimensions ny and (nx * nprocs)
          global_nx = nx * nprocs
          global_ny = ny

          ! add a global attribute
          attr_len = 24
          err = nfmpi_put_att_text(ncid, NF_GLOBAL, "history", attr_len,
     +                             "Mon Aug 13 21:27:48 2018")
          call check(err, 'In nfmpi_put_att_text: ')

          ! define dimensions time, Y, and X
          err = nfmpi_def_dim(ncid, "time", NFMPI_UNLIMITED, dim_t)
          call check(err, 'In nfmpi_def_dim time: ')
          err = nfmpi_def_dim(ncid, "Y",    global_ny, dim_y)
          call check(err, 'In nfmpi_def_dim Y: ')
          err = nfmpi_def_dim(ncid, "X",    global_nx, dim_x)
          call check(err, 'In nfmpi_def_dim X: ')

          ! define a 2D record variable of float type
          dimid(1) = dim_x
          dimid(2) = dim_t
          err = nfmpi_def_var(ncid,"rec_var",NF_FLOAT,2,dimid,rec_var)
          call check(err, 'In nfmpi_def_var rec_var: ')

          ! define a 2D fixed-size variable of float type
          dimid(1) = dim_x
          dimid(2) = dim_y
          err = nfmpi_def_var(ncid,"fix_var",NF_FLOAT,2,dimid,fix_var)
          call check(err, 'In nfmpi_def_var fix_var: ')

          ! exit define mode and enter data mode
          err = nfmpi_enddef(ncid)
          call check(err, 'In nfmpi_enddef: ')

          ! write to the 1st record of record variable (the variable
          ! with expandable dimension)
          do j=1, nx
             buf(j,1) = 1.0 * rank + 100.0
          enddo

          ! Note that in Fortran, array indices start with 1
          start(1) = nx * rank + 1
          start(2) = 1
          count(1) = nx
          count(2) = 1

          err = nfmpi_put_vara_double_all(ncid,rec_var,start,count, buf)
          call check(err, 'In nfmpi_put_vara_double_all: ')

          ! write to fixed-size variable (the variable with no expandable
          ! dimension)
          do i=1, ny
          do j=1, nx
             buf(j,i) = 1.0 * rank
          enddo
          enddo

          start(1) = nx * rank + 1
          start(2) = 1
          count(1) = nx
          count(2) = ny
          err = nfmpi_put_vara_double_all(ncid,fix_var,start,count, buf)
          call check(err, 'In nfmpi_put_vara_double_all: ')

          ! write a new record: 2nd record
          do j=1, nx
             buf(j,1) = 1.0 * rank + 200.0
          enddo

          start(1) = nx * rank + 1
          start(2) = 2  ! 2nd record
          count(1) = nx
          count(2) = 1

          err = nfmpi_put_vara_double_all(ncid,rec_var,start,count, buf)
          call check(err, 'In nfmpi_put_vara_double_all: ')

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

      end ! subroutine pnetcdf_write

      subroutine pnetcdf_read(filename)
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'

          character*256 filename

          ! local variables
          integer i, j, err, nprocs, rank, str_len
          integer ncid, rec_var, fix_var, dim_t, dim_y, dim_x
          integer*8 nx, ny, global_nx, global_ny, time_len, local_nx
          integer*8 start(2), count(2)
          PARAMETER(nx=3, ny=4)
          character*256 str_att
          double precision buf(nx,ny), expect

          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, err)

          ! open file for reading
          err = nfmpi_open(MPI_COMM_WORLD, filename, NF_NOWRITE,
     +                     MPI_INFO_NULL, ncid)
          call check(err, 'In nfmpi_open: ')

          ! nfmpi_open automatically enters data mode

          ! inquire dimensions time, Y, and X
          err = nfmpi_inq_dimid(ncid, "time", dim_t)
          call check(err, 'In nfmpi_inq_dimid time: ')
          err = nfmpi_inq_dimid(ncid, "Y",    dim_y)
          call check(err, 'In nfmpi_inq_dimid Y: ')
          err = nfmpi_inq_dimid(ncid, "X",    dim_x)
          call check(err, 'In nfmpi_inq_dimid X: ')
          err = nfmpi_inq_dimlen(ncid, dim_t, time_len)
          call check(err, 'In nfmpi_inq_dimlen time: ')
          err = nfmpi_inq_dimlen(ncid, dim_y, global_ny)
          call check(err, 'In nfmpi_inq_dimlen Y: ')
          err = nfmpi_inq_dimlen(ncid, dim_x, global_nx)
          call check(err, 'In nfmpi_inq_dimlen X: ')

          local_nx = global_nx / nprocs

          ! inquire global attribute history
          err = nfmpi_inq_attlen(ncid, NF_GLOBAL, "history", str_len)
          call check(err, 'In nfmpi_inq_attlen: ')
          err = nfmpi_get_att_text(ncid, NF_GLOBAL, "history", str_att)
          call check(err, 'In nfmpi_get_att_text: ')

          ! inquire variable IDs
          err = nfmpi_inq_varid(ncid, "rec_var", rec_var)
          call check(err, 'In nfmpi_inq_varid rec_var: ')
          err = nfmpi_inq_varid(ncid, "fix_var", fix_var)
          call check(err, 'In nfmpi_inq_varid fix_var: ')

          ! read the 1st record of the record variable
          do j=1, nx
             buf(j,1) = -1.0
          enddo

          ! Note that in Fortran, array indices start with 1
          start(1) = nx * rank + 1
          start(2) = 1
          count(1) = nx
          count(2) = 1

          err = nfmpi_get_vara_double_all(ncid,rec_var,start,count, buf)
          call check(err, 'In nfmpi_get_vara_double_all: ')

          ! check read contents */
 100      format(A,I3,A,F4.1,A,F4.1)
          do i=1, local_nx
              expect=1.0 * rank + 100.0
              if (buf(i,1) .NE. expect) then
                  print 100,"Read error buf(",i,",1) expect ", expect,
     +            " but got ", buf(i,1)
                  STOP 10
              endif
          enddo

          ! read fixed-size variable (the variable with no expandable dimension)
          do i=1, ny
          do j=1, nx
             buf(j,i) = -1.0
          enddo
          enddo

          start(1) = nx * rank + 1
          start(2) = 1
          count(1) = nx
          count(2) = ny
          err = nfmpi_get_vara_double_all(ncid,fix_var,start,count, buf)
          call check(err, 'In nfmpi_get_vara_double_all: ')

          ! check read contents */
          do j=1, ny
          do i=1, local_nx
              expect=1.0 * rank
              if (buf(j,i) .NE. expect) then
                  print 100,"Read error buf(",j,",",i,") expect ",
     +            expect, " but got ", buf(j,i)
                  STOP 20
              endif
          enddo
          enddo

          ! read the 2nd record of the record variable
          do j=1, nx
             buf(j,1) = -1.0
          enddo

          start(1) = nx * rank + 1
          start(2) = 2  ! 2nd record
          count(1) = nx
          count(2) = 1

          err = nfmpi_get_vara_double_all(ncid,rec_var,start,count, buf)
          call check(err, 'In nfmpi_get_vara_double_all: ')

          ! check read contents */
 200      format(A,I3,A,I3,A,F4.1,A,F4.1)
          do i=1, local_nx
              expect=1.0 * rank + 200.0
              if (buf(i,1) .NE. expect) then
                  print 200,"Read error buf(",i,",2) expect ", expect,
     +            " but got ", buf(i,1)
                  STOP 30
              endif
          enddo

          ! close the file
          err = nfmpi_close(ncid)
          call check(err, 'In nfmpi_close: ')

      end ! subroutine pnetcdf_read

      program main
          implicit none
          include 'mpif.h'
          include 'pnetcdf.inc'

          character*256 filename, cmd
          integer err, ierr, rank, get_args, dummy
          integer*8 malloc_size, sum_size
          logical verbose

          call MPI_Init(err)
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, err)

          ! take filename from command-line argument if there is any
          if (rank .EQ. 0) then
              filename = "testfile.nc"
              ierr = get_args(2, cmd, filename, verbose, dummy)
          endif
          call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
          if (ierr .EQ. 0) goto 999

          call MPI_Bcast(filename, 256, MPI_CHARACTER, 0,
     +                   MPI_COMM_WORLD, err)

          ! test classic CDF-1 file format
          call pnetcdf_write(filename, 0)
          call pnetcdf_read(filename)

          ! test classic CDF-2 file format
          call pnetcdf_write(filename, NF_64BIT_OFFSET)
          call pnetcdf_read(filename)

          ! test classic CDF-5 file format
          call pnetcdf_write(filename, NF_64BIT_DATA)
          call pnetcdf_read(filename)

          ! check if there is any PnetCDF internal malloc residue
 998      format(A,I13,A)
          err = nfmpi_inq_malloc_size(malloc_size)
          if (err .EQ. NF_NOERR) then
              call MPI_Reduce(malloc_size, sum_size, 1, MPI_INTEGER8,
     +                        MPI_SUM, 0, MPI_COMM_WORLD, err)
              if (rank .EQ. 0 .AND. sum_size .GT. 0)
     +            print 998,
     +            'heap memory allocated by PnetCDF internally has ',
     +            sum_size/1048576, ' MiB yet to be freed'
          endif

 999      call MPI_Finalize(err)
          ! call EXIT(0) ! EXIT() is a GNU extension
      end ! program main

