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

      double precision time_io, time_begin

      integer, parameter :: local_blocks = INT(0.8*maxblocks)

! initialize MPI and get the rank and size
      call MPI_INIT(ierr)
      
      call MPI_Comm_Rank (MPI_Comm_World, MyPE, ierr)
      call MPI_Comm_Size (MPI_Comm_World, NumPEs, ierr)

      MasterPE = 0

      if (MyPE .EQ. MasterPE) then

! print some initialization information
#ifdef BGL
         print *, 'Parallel netCDF test on IBM BG/L'
#endif
         print *, NumPEs, ' processors'
         print *, 'nxb, nyb, nzb = ', nxb, nyb, nzb
         print *, 'nguard = ', nguard
         print *, 'number of blocks = ', local_blocks
         print *, 'nvar = ', nvar
      endif


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

#ifdef CHIBA
      basenm = "pvfs:/stopvfsvol/flash_io_test_"
#endif      

#ifdef VIRGO
      basenm = "/mnt/pvfs/flash_io_bench/ncmpi/flash_io_test_"
#endif

#ifdef BGL
	basenm = "/pvfs-surveyor/robl/flash_io_test_ncmpi_"
#endif


!---------------------------------------------------------------------------
! netCDF checkpoint file
!---------------------------------------------------------------------------
      if (MyPE .EQ. MasterPE) then
         print *, ' '
         print *, 'doing  parallel I/O to a single file'
      endif

      time_begin = MPI_Wtime()
      call checkpoint_wr_ncmpi_par(0,0.e0)
      time_io = MPI_Wtime() - time_begin

      if (MyPE .EQ. MasterPE) then
         print *, 'time to output = ', time_io
      endif

!---------------------------------------------------------------------------
! netCDF plotfile -- no corners
!---------------------------------------------------------------------------
      if (MyPE .EQ. MasterPE) then
         print *, ' '
         print *, 'writing an netCDF plotfile --  no corners'
      endif

      time_begin = MPI_Wtime()
      call plotfile_ncmpi_par(0,0.e0,.false.)
      time_io = MPI_Wtime() - time_begin
    
      if (MyPE .EQ. MasterPE) then
         print *, 'time to output = ', time_io
      endif

!---------------------------------------------------------------------------
! netCDF plotfile -- corners
!---------------------------------------------------------------------------
      if (MyPE .EQ. MasterPE) then
         print *, ' '
         print *, 'writing an netCDF plotfile --  corners'
      endif

      time_begin = MPI_Wtime()
      call plotfile_ncmpi_par(0,0.e0,.true.)
      time_io = MPI_Wtime() - time_begin
    
      if (MyPE .EQ. MasterPE) then
         print *, 'time to output = ', time_io
      endif


      call MPI_Finalize(ierr)
      end





