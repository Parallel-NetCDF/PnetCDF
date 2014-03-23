!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$
!
! This file is taken from netcdf_file.f90 with changes for PnetCDF use
!  
!

! This is part of the netCDF F90 API, or. Copyright 2006 UCAR. See COPYRIGHT file
! for details.

! This file contains the netcdf file functions that are shared by
! netcdf-3 and netcdf-4.

! $Id$
! -------
function nf90mpi_inq_libvers()
  character(len = 80) :: nf90mpi_inq_libvers

  nf90mpi_inq_libvers = nfmpi_inq_libvers()
end function nf90mpi_inq_libvers
! -------
function nf90mpi_strerror(ncerr)
  integer, intent( in) :: ncerr
  character(len = 80)  :: nf90mpi_strerror

  nf90mpi_strerror = nfmpi_strerror(ncerr)
end function nf90mpi_strerror
! -------
!
! File level control routines:
!
! function nf90mpi_inq_base_pe(ncid, pe)
!   integer, intent( in) :: ncid
!   integer, intent(out) :: pe
!   integer              :: nf90mpi_inq_base_pe
! 
!   nf90mpi_inq_base_pe = nfmpi_inq_base_pe(ncid, pe)
! end function nf90mpi_inq_base_pe
! -------
! function nf90mpi_set_base_pe(ncid, pe)
!   integer, intent( in) :: ncid, pe
!   integer              :: nf90mpi_set_base_pe
! 
!   nf90mpi_set_base_pe = nfmpi_set_base_pe(ncid, pe)
! end function nf90mpi_set_base_pe
! -------
function nf90mpi_create(mpi_comm, path, cmode, mpi_info, ncid)
  character (len = *), intent( in) :: path
  integer,             intent( in) :: mpi_comm, cmode, mpi_info
  integer,             intent(out) :: ncid
  integer                          :: nf90mpi_create

  nf90mpi_create = nfmpi_create(mpi_comm, path, cmode, mpi_info, ncid)
end function nf90mpi_create
! -------
function nf90mpi_open(mpi_comm, path, omode, mpi_info, ncid)
  character (len = *), intent( in) :: path
  integer,             intent( in) :: mpi_comm, omode, mpi_info
  integer,             intent(out) :: ncid
  integer                          :: nf90mpi_open

  nf90mpi_open = nfmpi_open(mpi_comm, path, omode, mpi_info, ncid)
end function nf90mpi_open
! -------
function nf90mpi_set_fill(ncid, fillmode, old_mode)
  integer, intent( in) :: ncid, fillmode 
  integer, intent(out) :: old_mode
  integer              :: nf90mpi_set_fill

  nf90mpi_set_fill = nfmpi_set_fill(ncid, fillmode, old_mode)
end function nf90mpi_set_fill
! -------
function nf90mpi_redef(ncid)
  integer, intent( in) :: ncid
  integer              :: nf90mpi_redef

  nf90mpi_redef = nfmpi_redef(ncid)
end function nf90mpi_redef
! -------
function nf90mpi_enddef(ncid)
  integer, intent( in) :: ncid
  integer              :: nf90mpi_enddef

  nf90mpi_enddef = nfmpi_enddef(ncid)
end function nf90mpi_enddef
! -------
function nf90mpi_sync(ncid)
  integer, intent( in) :: ncid
  integer              :: nf90mpi_sync

  nf90mpi_sync = nfmpi_sync(ncid)
end function nf90mpi_sync
! -------
function nf90mpi_abort(ncid)
  integer, intent( in) :: ncid
  integer              :: nf90mpi_abort

  nf90mpi_abort = nfmpi_abort(ncid)
end function nf90mpi_abort
! -------
function nf90mpi_close(ncid)
  integer, intent( in) :: ncid
  integer              :: nf90mpi_close

  nf90mpi_close = nfmpi_close(ncid)
end function nf90mpi_close
! -------
function nf90mpi_delete(name, mpi_info)
  character(len = *), intent( in) :: name
  integer                         :: nf90mpi_delete
  integer,            intent( in) :: mpi_info

  nf90mpi_delete = nfmpi_delete(name, mpi_info)
end function nf90mpi_delete

!
! A single file level inquiry routine 
! 
function nf90mpi_inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum)
  integer,           intent( in) :: ncid
  integer, optional, intent(out) :: nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum
  integer                        :: nf90mpi_inquire

  integer :: nDims, nVars, nGAtts, unlimDimId, frmt

  nf90mpi_inquire = nfmpi_inq(ncid, nDims, nVars, nGAtts, unlimDimId)
  if(present(nDimensions))    nDimensions    = nDims 
  if(present(nVariables))     nVariables     = nVars
  if(present(nAttributes))    nAttributes    = nGAtts
  if(present(unlimitedDimId)) unlimitedDimId = unlimDimId
  if(present(formatNum)) then
     nf90mpi_inquire = nfmpi_inq_format(ncid, frmt)
     formatNum = frmt
  endif
end function nf90mpi_inquire

! -------
function nf90mpi_inq_striping(ncid, striping_size, striping_count)
  integer,            intent( in) :: ncid
  integer,            intent(out) :: striping_size
  integer,            intent(out) :: striping_count
  integer                         :: nf90mpi_inq_striping

  nf90mpi_inq_striping = nfmpi_inq_striping(ncid, striping_size, striping_count)
end function nf90mpi_inq_striping

! -------
function nf90mpi_begin_indep_data(ncid)
  integer,            intent( in) :: ncid
  integer                         :: nf90mpi_begin_indep_data

  nf90mpi_begin_indep_data = nfmpi_begin_indep_data(ncid)
end function nf90mpi_begin_indep_data

! -------
function nf90mpi_end_indep_data(ncid)
  integer,            intent( in) :: ncid
  integer                         :: nf90mpi_end_indep_data

  nf90mpi_end_indep_data = nfmpi_end_indep_data(ncid)
end function nf90mpi_end_indep_data

! -------
function nf90mpi_inq_put_size(ncid, put_size)
  integer,                          intent( in) :: ncid
  integer (kind = MPI_OFFSET_KIND), intent(out) :: put_size
  integer                                       :: nf90mpi_inq_put_size

  nf90mpi_inq_put_size = nfmpi_inq_put_size(ncid, put_size)
end function nf90mpi_inq_put_size

! -------
function nf90mpi_inq_get_size(ncid, get_size)
  integer,                          intent( in) :: ncid
  integer (kind = MPI_OFFSET_KIND), intent(out) :: get_size
  integer                                       :: nf90mpi_inq_get_size

  nf90mpi_inq_get_size = nfmpi_inq_get_size(ncid, get_size)
end function nf90mpi_inq_get_size

! -------
function nf90mpi_inq_header_size(ncid, h_size)
  integer,                          intent( in) :: ncid
  integer (kind = MPI_OFFSET_KIND), intent(out) :: h_size
  integer                                       :: nf90mpi_inq_header_size

  nf90mpi_inq_header_size = nfmpi_inq_header_size(ncid, h_size)
end function nf90mpi_inq_header_size

! -------
function nf90mpi_inq_header_extent(ncid, h_extent)
  integer,                          intent( in) :: ncid
  integer (kind = MPI_OFFSET_KIND), intent(out) :: h_extent
  integer                                       :: nf90mpi_inq_header_extent

  nf90mpi_inq_header_extent = nfmpi_inq_header_extent(ncid, h_extent)
end function nf90mpi_inq_header_extent

! -------
function nf90mpi_inq_varoffset(ncid, varid, offset)
  integer,                          intent( in) :: ncid, varid
  integer (kind = MPI_OFFSET_KIND), intent(out) :: offset
  integer                                       :: nf90mpi_inq_varoffset

  nf90mpi_inq_varoffset = nfmpi_inq_varoffset(ncid, varid, offset)
end function nf90mpi_inq_varoffset

! -------
function nf90mpi_inq_nreqs(ncid, nreqs)
  integer, intent( in) :: ncid
  integer, intent(out) :: nreqs
  integer              :: nf90mpi_inq_nreqs

  nf90mpi_inq_nreqs = nfmpi_inq_nreqs(ncid, nreqs)
end function nf90mpi_inq_nreqs

! -------
function nf90mpi_get_file_info(ncid, mpi_info)
  integer, intent( in) :: ncid
  integer, intent(out) :: mpi_info
  integer              :: nf90mpi_get_file_info

  nf90mpi_get_file_info = nfmpi_get_file_info(ncid, mpi_info)
end function nf90mpi_get_file_info

