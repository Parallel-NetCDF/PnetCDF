!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$
!
! This file is taken from netcdf_variables.f90 with changes for PnetCDF use
!  
!

  ! ----- 
  ! Variable definitions and inquiry
  ! ----- 
  function nf90mpi_def_var_Scalar(ncid, name, xtype, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer,               intent(out) :: varid
    integer                            :: nf90mpi_def_var_Scalar
    
    ! Dummy - shouldn't get used
    integer, dimension(1) :: dimids
    
    nf90mpi_def_var_Scalar = nfmpi_def_var(ncid, name, xtype, 0, dimids, varid)
  end function nf90mpi_def_var_Scalar
  ! ----- 
  function nf90mpi_def_var_oneDim(ncid, name, xtype, dimids, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer,               intent( in) :: dimids
    integer,               intent(out) :: varid
    integer                            :: nf90mpi_def_var_oneDim
    
    integer, dimension(1) :: dimidsA
    dimidsA(1) = dimids
    nf90mpi_def_var_oneDim = nfmpi_def_var(ncid, name, xtype, 1, dimidsA, varid)
  end function nf90mpi_def_var_oneDim
  ! ----- 
  function nf90mpi_def_var_ManyDims(ncid, name, xtype, dimids, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer, dimension(:), intent( in) :: dimids
    integer,               intent(out) :: varid
    integer                            :: nf90mpi_def_var_ManyDims
    
    nf90mpi_def_var_ManyDims = nfmpi_def_var(ncid, name, xtype, size(dimids), dimids, varid)
  end function nf90mpi_def_var_ManyDims
  ! ----- 
  function nf90mpi_inq_varid(ncid, name, varid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(out) :: varid
    integer                          :: nf90mpi_inq_varid
    
    nf90mpi_inq_varid = nfmpi_inq_varid(ncid, name, varid)
  end function nf90mpi_inq_varid
  ! ----- 
  function nf90mpi_inquire_variable(ncid, varid, name, xtype, ndims, dimids, nAtts)
    integer,                         intent( in) :: ncid, varid
    character (len = *),   optional, intent(out) :: name
    integer,               optional, intent(out) :: xtype, ndims 
    integer, dimension(:), optional, intent(out) :: dimids
    integer,               optional, intent(out) :: nAtts
    integer                                      :: nf90mpi_inquire_variable
    
    ! Local variables
    character (len = nf90_max_name)       :: varName
    integer                               :: externalType, numDimensions
    integer, dimension(nf90_max_var_dims) :: dimensionIDs
    integer                               :: numAttributes
    
    nf90mpi_inquire_variable = nfmpi_inq_var(ncid, varid, varName, externalType, &
                                       numDimensions, dimensionIDs, numAttributes)
    if (nf90mpi_inquire_variable == nf90_noerr) then
        if(present(name))   name                   = trim(varName)
        if(present(xtype))  xtype                  = externalType
        if(present(ndims))  ndims                  = numDimensions
        if(present(dimids)) then
            if (size(dimids) .ge. numDimensions) then
                dimids(:numDimensions) = dimensionIDs(:numDimensions)
            else
                nf90mpi_inquire_variable = nf90_einval
            endif
        endif
        if(present(nAtts))  nAtts                  = numAttributes
    endif
  end function nf90mpi_inquire_variable
  ! ----- 
  function nf90mpi_rename_var(ncid, varid, newname)
    integer,             intent( in) :: ncid, varid
    character (len = *), intent( in) :: newname
    integer                          :: nf90mpi_rename_var
    
    nf90mpi_rename_var = nfmpi_rename_var(ncid, varid, newname)
  end function nf90mpi_rename_var
  ! ----- 
