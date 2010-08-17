! this test case came from Annette Koonts at PNNL
!
! test creates a dataset with several (large) record variables, but also a
! smaller 'time' record variable.
!
! on 64 bit systems or on BlueGene (where MPI_AINT is 32 bits), this test will
! read the smaller time variable OK, but on systems with a 32 bit MPI_AINT, the
! time variable will be garbled.

      program main
      use pnetcdf
      implicit none
      include "mpif.h"


! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids
      integer  time_dim
      integer  cells_dim
      integer  interfaces_dim

! dimension lengths
      integer  time_len
      integer*8  cells_len
      integer  interfaces_len

      parameter (time_len = NF_UNLIMITED)
      parameter (cells_len = 41943042)
      parameter (interfaces_len = 26)

! variable ids
      integer  time_id
      integer  interfaces_id
      integer  pressure_id

      integer (kind=MPI_OFFSET_KIND) :: longlen

      integer (kind=MPI_OFFSET_KIND) :: start1d(1)
      integer (kind=MPI_OFFSET_KIND) :: count1d(1)

! rank (number of dimensions) for each variable
      integer  time_rank
      integer  interfaces_rank
      integer  pressure_rank

      parameter (time_rank = 1)
      parameter (interfaces_rank = 1)
      parameter (pressure_rank = 3)

! variable shapes
      integer  time_dims(time_rank)
      integer  interfaces_dims(interfaces_rank)
      integer  pressure_dims(pressure_rank)

! data variables
      real  interfaces(interfaces_len)

      integer  myid, ierr, n
      integer  numprocs

      integer*8 i8_size

      integer(kind=MPI_OFFSET_KIND) :: time_start(1), time_count(1)

      double precision  time(4)
!      data time /0., 20., 40., 60./

      data interfaces /2685.8359, 671.81, 495.91, 425.10001, 393.42999, 
     1377.5, 367.59, 360.06, 353.85999, 348.66, 342.5, 336, 328.5, 320, 
     2310, 300, 290, 280, 270, 260, 250, 240, 230, 220, 210, 199.10001/


! attribute vectors
! enter define mode
!      iret = nf_create('pressure_19010101_000000.nc', OR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

      iret = nfmpi_create( MPI_COMM_WORLD,
     2                       'pressure_test.nc', 
     3                       IOR(NF_CLOBBER,NF_64BIT_DATA), 
     4                       MPI_INFO_NULL,
     4                       ncid)

        call check_err(iret)
  
! define dimensions

        iret = nfmpi_def_dim(ncid, 'time', NFMPI_UNLIMITED, time_dim)
        call check_err(iret)

        i8_size = 41943042

        iret = nfmpi_def_dim(ncid, 'cells', i8_size, cells_dim)
        call check_err(iret)

        i8_size = 26
        iret = nfmpi_def_dim(ncid, 'interfaces', 
     1                       i8_size, interfaces_dim)
        call check_err(iret)
! define variables
        time_dims(1) = time_dim

        iret = nfmpi_def_var(ncid, 'time', NF_DOUBLE, 
     1                       time_rank, time_dims, 
     1                       time_id)
        call check_err(iret)
        interfaces_dims(1) = interfaces_dim

        iret = nfmpi_def_var(ncid, 'interfaces', NF_REAL, 
     1                       interfaces_rank, 
     1                       interfaces_dims, interfaces_id)
        call check_err(iret)

        pressure_dims(3) = time_dim
        pressure_dims(2) = cells_dim
        pressure_dims(1) = interfaces_dim
        iret = nfmpi_def_var(ncid, 
     1                     'pressure', 
     2                     NF_REAL, 
     3                     pressure_rank, 
     4                     pressure_dims, 
     5                     pressure_id)

        call check_err(iret)
! assign attributes

        longlen = 4
        iret = nfmpi_put_att_text(ncid, time_id, 'long_name', 
     1                            longlen, 'Time')
        call check_err(iret)
        longlen = 21
        iret = nfmpi_put_att_text(ncid, time_id, 'units', 
     1                            longlen, 
     2                           'days since 1901-01-01')
        call check_err(iret)

        longlen = 41
        iret = nfmpi_put_att_text(ncid, interfaces_id, 'long_name', 
     1                            longlen,
     1                     'Vertical interfaces, in terms of pressure')
        call check_err(iret)

        longlen = 2
        iret = nfmpi_put_att_text(ncid, interfaces_id, 'units', 
     1                            longlen, 'Pa')
        call check_err(iret)

        longlen = 8
        iret = nfmpi_put_att_text(ncid, pressure_id, 'long_name', 
     1                            longlen, 
     1                       'Pressure')
        call check_err(iret)

        longlen = 2
        iret = nfmpi_put_att_text(ncid, pressure_id, 'units', 
     1                            longlen, 'Pa')
        call check_err(iret)

! leave define mode
        iret = nfmpi_enddef(ncid)
        call check_err(iret)

        start1d(1) = 1
        count1d(1) = 26

! store interfaces
        iret = nfmpi_put_vara_real_all(ncid, interfaces_id,
     1                                 start1d, count1d, interfaces)
        call check_err(iret)

        time(1) = 0.0
        time(2) = 20.0
        time(3) = 40.0
        time(4) = 60.0
      
! this test is tricky because it writes out the time variable one at a time.
! This element-at-a-time workload does not actually exercise the tricky 32 bit
! MPI_AINT problem, so the issue only shows up at read time.

        do n = 1, 4
          time_start(1) = n
          if(myid .eq. 0) then
            time_count(1) = 1
          else
            time_count(1) = 0
          endif

          iret = nfmpi_put_vara_double_all(ncid, time_id, 
     1                                     time_start, time_count,
     1                                     time(n))
          call check_err(iret)
        enddo

        iret = nfmpi_close(ncid)

        call MPI_Barrier (MPI_COMM_WORLD, iret)

! todo: insert code to re-open dataset, read time variable all at onece
! 
      iret = nfmpi_open ( MPI_COMM_SELF,
     1                   'pressure_test.nc',
     1                   IOR(NF_CLOBBER,NF_64BIT_DATA),
     1                   MPI_INFO_NULL,
     1                   ncid)
      call check_err(iret)

      iret = nfmpi_inq_varid(ncid, 'time', time_id);
      call check_err(iret)

      ! deliberately want all processes to end up with the full time array
      time_start(1) = 1
      time_count(1) = 4
      iret = nfmpi_get_vara_double_all(ncid, time_id,
     1                               time_start, time_count, time);
      call check_err(iret)

      iret = nfmpi_close(ncid)
      call check_err(iret)

      write(6,*) "Time array: ", time
!      if ( (time(1) .eq. 0) .and. (time(2) .eq. 20.0)
!     &             .and. (time(3) .eq. 40.0) .and. (time(4) .eq. 60))
!           write(6,*) " No Errors"
!      else
!           write(6,*) "Error: time array was ", time
!      endif

      call MPI_FINALIZE(ierr)
      end program main
       
      subroutine writerecs(ncid,time_id,pressure_id)

      include "mpif.h"
#     include "pnetcdf.inc"
       
! netCDF id
      integer  ncid
! variable ids
      integer  time_id
      integer  pressure_id

! error status return
      integer  iret
       
! netCDF dimension sizes for dimensions used with record variables
      integer  cells_len
      parameter (cells_len = 41943042)
      integer  interfaces_len
      parameter (interfaces_len = 26)
       
! rank (number of dimensions) for each variable
      integer  time_rank
      integer  pressure_rank
      parameter (time_rank = 1)
      parameter (pressure_rank = 3)
! starts and counts for array sections of record variables
      integer(kind=MPI_OFFSET_KIND) :: time_start(1), time_count(1)
      integer  pressure_start(pressure_rank), pressure_count(pressure_ra
     1nk)
       
! data variables
       
      integer  time_nr
      parameter (time_nr = 4)
       
      integer  pressure_nr
      parameter (pressure_nr = 1)
!      real  pressure(interfaces_len, cells_len, pressure_nr)
       
      double precision  time(time_nr)
      data time /0., 20., 40., 60./
      

!      pressure = NF_FILL_FLOAT

! store time

      do n = 1, 4
        time_start(1) = n
        time_count(1) = 1
        iret = nfmpi_put_vara_double_all(ncid, time_id, 
     1                                   time_start, time_count,
     1                                   time)

        call check_err(iret)
      enddo
       
      end subroutine writerecs
       
      subroutine check_err(iret)

#     include "pnetcdf.inc"

      integer iret

      if (iret .ne. NF_NOERR) then
          print *, nfmpi_strerror(iret)
          stop
      endif
      end subroutine check_err
