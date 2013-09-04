      program flash_benchmark_io
!
! This is a sample program that setups the FLASH data structures and 
! drives the I/O routines.  It is intended for benchmarking the I/O
! performance.
! 

! the main data structures are contained in common blocks, defined in the
! include files

#include "common.fh"

      integer ierr
      integer i

      double precision time_io(3), time_begin
      double precision chk_io, corner_io, nocorner_io
      double precision checkpoint_wr_ncmpi_par
      double precision plotfile_ncmpi_par

      integer, parameter :: local_blocks = INT(0.8*maxblocks)

! initialize MPI and get the rank and size
      call MPI_INIT(ierr)
      
      call MPI_Comm_Rank (MPI_Comm_World, MyPE, ierr)
      call MPI_Comm_Size (MPI_Comm_World, NumPEs, ierr)

      MasterPE = 0

! put ~100 blocks on each processor -- make it vary a little, since it does
! in the real application.  This is the maximum that we can fit on Blue 
! Pacific comfortably.
      lnblocks = local_blocks + mod(MyPE,3)

! just fill the tree stucture with dummy information -- we are just going to
! dump it out
      size(:,:) = 0.5e0
      lrefine(:) = 1
      nodetype(:) = 1
      refine(:) = .FALSE.
      derefine(:) = .FALSE.
      parent(:,:) = -1
      child(:,:,:) = -1
      coord(:,:) = 0.25e0
      bnd_box(:,:,:) = 0.e0
      neigh(:,:,:) = -1
      empty(:) = 0

! initialize the unknowns with the index of the variable
      do i = 1, nvar
        unk(i,:,:,:,:) = float(i)
      enddo

! setup the file properties
      basenm = "flash_io_test_"

!---------------------------------------------------------------------------
! netCDF checkpoint file
!---------------------------------------------------------------------------
      time_begin = MPI_Wtime()
      chk_io = checkpoint_wr_ncmpi_par(0,0.e0)
      time_io(1) = MPI_Wtime() - time_begin

!---------------------------------------------------------------------------
! netCDF plotfile -- no corners
!---------------------------------------------------------------------------
      time_begin = MPI_Wtime()
      nocorner_io = plotfile_ncmpi_par(0,0.e0,.false.)
      time_io(2) = MPI_Wtime() - time_begin
    
!---------------------------------------------------------------------------
! netCDF plotfile -- corners
!---------------------------------------------------------------------------
      time_begin = MPI_Wtime()
      corner_io = plotfile_ncmpi_par(0,0.e0,.true.)
      time_io(3) = MPI_Wtime() - time_begin
    
      call report_io_performance(local_blocks, time_io, chk_io, &
                                 corner_io, nocorner_io)

      call MPI_Finalize(ierr)
      end


!---------------------------------------------------------------------------
! print I/O performance numbers
!---------------------------------------------------------------------------
      subroutine report_io_performance(local_blocks, time_io, chk_io, &
                                       corner_io, nocorner_io)

#include "common.fh"

       integer local_blocks
       double precision time_io(3)
       double precision chk_io, corner_io, nocorner_io

       ! local variables
       integer ierr
       double precision tmax(3), time_total, io_amount, bw

       call MPI_Reduce(chk_t, tmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       chk_t(:) = tmax(:)

       call MPI_Reduce(corner_t, tmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       corner_t(:) = tmax(:)

       call MPI_Reduce(nocorner_t, tmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       nocorner_t(:) = tmax(:)

       call MPI_Reduce(time_io, tmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MasterPE, MPI_COMM_WORLD, ierr)

       call MPI_Reduce(chk_io, bw, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       chk_io = bw
       call MPI_Reduce(nocorner_io, bw, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       nocorner_io = bw
       call MPI_Reduce(corner_io, bw, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MasterPE, MPI_COMM_WORLD, ierr)
       corner_io = bw

      if (MyPE .EQ. MasterPE) then
          time_total = tmax(1) + tmax(2) + tmax(3)

          io_amount = chk_io + nocorner_io + corner_io
          bw = io_amount / 1048576.0
          io_amount = bw
          bw = bw / time_total

          print *, 'File base name = ', trim(basenm)

          print 2008, nguard
          print 2009, local_blocks
          print 2010, nvar
          print 2011, tmax(1)
          print 2018, chk_t(1)
          print 2019, chk_t(2)
          print 2020, chk_t(3)
          print 2021, chk_io/1048576
          print 2012, tmax(2)
          print 2018, nocorner_t(1)
          print 2019, nocorner_t(2)
          print 2020, nocorner_t(3)
          print 2021, nocorner_io/1048576
          print 2013, tmax(3)
          print 2018, corner_t(1)
          print 2019, corner_t(2)
          print 2020, corner_t(3)
          print 2021, corner_io/1048576
          print 2014, io_amount
          print 2015
          print 2016
          print 2017, NumPEs, nxb, nyb, nzb, time_total, bw
      endif

 2008 format(' number of guards      : ',I6)
 2009 format(' number of blocks      : ',I6)
 2010 format(' number of variables   : ',I6)
 2011 format(' checkpoint time       : ',F16.2, '  sec')
 2012 format(' plot no corner        : ',F16.2, '  sec')
 2013 format(' plot    corner        : ',F16.2, '  sec')
 2014 format(' Total I/O amount      : ',F16.2, '  MiB')
 2015 format(' -------------------------------------------------------')
 2016 format(' nproc    array size      exec (sec)   bandwidth (MiB/s)')
 2017 format(I5, 3x, i3,' x ',i3,' x ',i3, 3x, F7.2 , 2x,F10.2 /)
 2018 format('        max header     : ',F16.2, '  sec')
 2019 format('        max unknown    : ',F16.2, '  sec')
 2020 format('        max close      : ',F16.2, '  sec')
 2021 format('            I/O amount : ',F16.2, '  MiB')
 2022 format(' procs    Local  array size  exec(sec)  write(MiB/s)   IO    req    comm   misc  clean  #sends  #recvs')
 2023 format('-------  ------------------  ---------  -----------   --    ---    ----   ----  -----  ------  ------')
 2024 format(I5,3x,i5,' x',i5,' x',i5,F9.2,F12.2,F9.2,F7.2,F7.2,F7.2,F7.2,I7,I8)
 2025 format(I5,3x,i5,' x',i5,' x',i5,F9.2,F12.2,F9.2,F7.2,F7.2,F7.2,F7.2 /)

      end subroutine report_io_performance




