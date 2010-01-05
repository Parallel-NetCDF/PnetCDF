      program Pnf_Test
C Test program thanks to From: John Tannahill <tannahill1@llnl.gov> 

      implicit none

      include "mpif.h"
#      include "pnetcdf.inc"


      INTEGER(KIND=MPI_OFFSET_KIND) TOTSIZ_3D(3)


      logical reorder

      logical isperiodic(3)

      integer comm_cart
      integer ierr
      integer lat_id, lev_id, lon_id
      integer ncid
      integer totpes
      integer tt_id

      integer dim_id(3)

      integer numpes(3)                   
C    number of PEs along axes;
C    determined by MPI where a zero is specified

      data totsiz_3d/10,20,30/
      data reorder /.false./
      data isperiodic /.false., .false., .false./
      data numpes/1,1,0/

      call MPI_Init (ierr)

      call MPI_Comm_Size (MPI_COMM_WORLD, totpes, ierr)

      call MPI_Dims_Create (totpes, 3, numpes, ierr)

      call MPI_Cart_Create (MPI_COMM_WORLD, 3, numpes, isperiodic,
     &                      reorder, comm_cart, ierr)


      ierr = Nfmpi_Create (comm_cart, "pnf_test.nc", 
     &                   NF_CLOBBER,
     &                   MPI_INFO_NULL, ncid)

      ierr = Nfmpi_Def_Dim
     &  (ncid, "level",     totsiz_3d(1), lev_id)
      ierr = Nfmpi_Def_Dim
     &  (ncid, "latitude",  totsiz_3d(2), lat_id)
      ierr = Nfmpi_Def_Dim
     &  (ncid, "longitude", totsiz_3d(3), lon_id)

      dim_id(1) = lev_id
      dim_id(2) = lat_id
      dim_id(3) = lon_id

      ierr = Nfmpi_Def_Var
     &  (ncid, "tt", NF_FLOAT, 3, dim_id, tt_id)

      ierr = Nfmpi_Enddef (ncid)

      ierr = Nfmpi_Close (ncid)


      call MPI_Comm_Free (comm_cart, ierr)

      call MPI_Finalize  (ierr)

      Write (6,10)

 10   format(" No Errors")

      Stop

      end program Pnf_Test
