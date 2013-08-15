!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$
!
! This file is taken from netcdf_constants.f90 with changes for PnetCDF use
!  
!
  !
  ! external netcdf data types:
  !
  integer, parameter, public :: &
    nf90_byte   = 1,            &
    nf90_int1   = nf90_byte,    &
    nf90_char   = 2,            &
    nf90_short  = 3,            &
    nf90_int2   = nf90_short,   &
    nf90_int    = 4,            &
    nf90_int4   = nf90_int,     &
    nf90_float  = 5,            &
    nf90_real   = nf90_float,   &
    nf90_real4  = nf90_float,   &
    nf90_double = 6,            &
    nf90_real8  = nf90_double,  &
    nf90_int64  = 10
                        
  !
  ! default fill values:
  !
  character (len = 1),           parameter, public :: &
    nf90_fill_char  = achar(0)
  integer (kind =  OneByteInt),  parameter, public :: &
    nf90_fill_byte  = -127,                           &
    nf90_fill_int1  = nf90_fill_byte
  integer (kind =  TwoByteInt),  parameter, public :: &
    nf90_fill_short = -32767,                         &
    nf90_fill_int2  = nf90_fill_short
  integer (kind = FourByteInt),  parameter, public :: &
    nf90_fill_int   = -2147483647
  real   (kind =  FourByteReal), parameter, public :: &
    nf90_fill_float = 9.9692099683868690e+36,         &
    nf90_fill_real  = nf90_fill_float,                &
    nf90_fill_real4 = nf90_fill_float
  real   (kind = EightByteReal), parameter, public :: &
    nf90_fill_double = 9.9692099683868690e+36,        &
    nf90_fill_real8  = nf90_fill_double
  integer (kind = EightByteInt), parameter, public :: &
    nf90_fill_int64  = -9223372036854775806_8

  !
  ! mode flags for opening and creating a netcdf dataset:
  !
  integer, parameter, public :: &
    nf90_nowrite   = 0,         &
    nf90_write     = 1,         &
    nf90_clobber   = 0,         &
    nf90_noclobber = 4,         &
    nf90_fill      = 0,         &
    nf90_nofill    = 256,       &
    nf90_64bit_offset    = 512,       &
    nf90_64bit_data    = 16,       &
    nf90_lock      = 1024,      &
    nf90_share     = 2048 
  
  integer, parameter, public ::  &
    nf90_sizehint_default = 0,   & 
    nf90_align_chunk      = -1 

  !
  ! size argument for defining an unlimited dimension:
  !
  integer (kind = EightByteInt), parameter, public :: nf90_unlimited = 0

  ! NULL request for non-blocking I/O APIs
  integer, parameter, public :: nf90_req_null = -1

  !
  ! global attribute id:
  !
  integer, parameter, public :: nf90_global = 0

  !
  ! implementation limits:
  !
  integer, parameter, public :: &
    nf90_max_dims     = 1024,    &
    nf90_max_attrs    = 8192,   &
    nf90_max_vars     = 8192,   &
    nf90_max_name     = 256,    &
    nf90_max_var_dims = 1024
  
  !
  ! error codes:
  !
  integer, parameter, public :: &
    nf90_noerr        = 0,      &
    nf90_ebadid       = -33,    &
    nf90_eexist       = -35,    &
    nf90_einval       = -36,    &
    nf90_eperm        = -37,    &
    nf90_enotindefine = -38,    &
    nf90_eindefine    = -39,    &
    nf90_einvalcoords = -40,    &
    nf90_emaxdims     = -41,    &
    nf90_enameinuse   = -42,    &
    nf90_enotatt      = -43,    &
    nf90_emaxatts     = -44,    &
    nf90_ebadtype     = -45,    &
    nf90_ebaddim      = -46,    &
    nf90_eunlimpos    = -47,    &
    nf90_emaxvars     = -48,    &
    nf90_enotvar      = -49,    &
    nf90_eglobal      = -50,    &
    nf90_enotnc       = -51,    &
    nf90_ests         = -52,    &
    nf90_emaxname     = -53,    &
    nf90_eunlimit     = -54,    &
    nf90_enorecvars   = -55,    &
    nf90_echar        = -56,    &
    nf90_eedge        = -57,    &
    nf90_estride      = -58,    &
    nf90_ebadname     = -59,    &
    nf90_erange       = -60,    &
    nf90_enomem       = -61,    &
    nf90_evarsize     = -62,    &
    nf90_edimsize     = -63,    &
    nf90_etrunc       = -64

  !
  ! error handling modes:
  !
  integer, parameter, public :: &
    nf90_fatal   = 1,           &
    nf90_verbose = 2

  !
  ! format version numbers:
  !
  integer, parameter, public :: &
    nf90_format_classic = 1,    &
    nf90_format_64bit = 2,      &
    nf90_format_netcdf4 = 3,    &
    nf90_format_netcdf4_classic = 4, &
    nf90_format_64bit_data = 5

  ! PnetCDF error codes start here
  integer, parameter, public :: &
      NF90_ESMALL                   = NF_ESMALL                   , & ! size of off_t too small for format
      NF90_ENOTINDEP                = NF_ENOTINDEP                , & ! Operation not allowed in collective data mode
      NF90_EINDEP                   = NF_EINDEP                   , & ! Operation not allowed in independent data mode
      NF90_EFILE                    = NF_EFILE                    , & ! Unknown error in file operation
      NF90_EREAD                    = NF_EREAD                    , & ! Unknown error in reading file
      NF90_EWRITE                   = NF_EWRITE                   , & ! Unknown error in writting to file
      NF90_EMULTIDEFINE             = NF_EMULTIDEFINE             , & ! NC definitions on multiprocesses conflict
      NF90_EOFILE                   = NF_EOFILE                   , & ! file open/creation failed
      NF90_EMULTITYPES              = NF_EMULTITYPES              , & ! Multiple types used in memory data
      NF90_EIOMISMATCH              = NF_EIOMISMATCH              , & ! Input/Output data amount mismatch
      NF90_ENEGATIVECNT             = NF_ENEGATIVECNT             , & ! Negative count is specified
      NF90_EUNSPTETYPE              = NF_EUNSPTETYPE              , & ! Unsupported etype in memory MPI datatype
      NF90_EDIMS_NELEMS_MULTIDEFINE = NF_EDIMS_NELEMS_MULTIDEFINE , & ! Different number of dim defines on multiprocesses conflict
      NF90_EDIMS_SIZE_MULTIDEFINE   = NF_EDIMS_SIZE_MULTIDEFINE   , & ! Different size of dim defines on multiprocesses conflict
      NF90_EVARS_NELEMS_MULTIDEFINE = NF_EVARS_NELEMS_MULTIDEFINE , & ! Different number of var defines on multiprocesses conflict
      NF90_EVARS_NDIMS_MULTIDEFINE  = NF_EVARS_NDIMS_MULTIDEFINE  , & ! Different dim number of var defines on multiprocesses conflict
      NF90_EVARS_DIMIDS_MULTIDEFINE = NF_EVARS_DIMIDS_MULTIDEFINE , & ! Different dimid defines on multiprocesses conflict
      NF90_EVARS_TYPE_MULTIDEFINE   = NF_EVARS_TYPE_MULTIDEFINE   , & ! Different type of var defines on multiprocesses conflict
      NF90_EVARS_LEN_MULTIDEFINE    = NF_EVARS_LEN_MULTIDEFINE    , & ! Different var lenght defines size on multiprocesses conflict
      NF90_EVARS_BEGIN_MULTIDEFINE  = NF_EVARS_BEGIN_MULTIDEFINE  , & ! Different var begin defines size on multiprocesses conflict
      NF90_ENUMRECS_MULTIDEFINE     = NF_ENUMRECS_MULTIDEFINE     , & ! Different number records on multiprocesses conflict
      NF90_EINVAL_REQUEST           = NF_EINVAL_REQUEST           , & ! invalid nonblocking request ID
      NF90_EAINT_TOO_SMALL          = NF_EAINT_TOO_SMALL          , & ! MPI_Aint not large enough to hold requested value
      NF90_ECMODE                   = NF_ECMODE                   , & ! file create modes are inconsistent among processes
      NF90_ENOTSUPPORT              = NF_ENOTSUPPORT              , & ! feature is not yet supported
      NF90_ENULLBUF                 = NF_ENULLBUF                 , & ! trying to attach a NULL buffer
      NF90_EPREVATTACHBUF           = NF_EPREVATTACHBUF           , & ! previous attached buffer is found
      NF90_ENULLABUF                = NF_ENULLABUF                , & ! no attached buffer is found
      NF90_EPENDINGBPUT             = NF_EPENDINGBPUT             , & ! pending bput is found, cannot detach buffer
      NF90_EINSUFFBUF               = NF_EINSUFFBUF               , & ! attached buffer is too small
      NF90_ENOENT                   = NF_ENOENT                   , & ! File does not exist when calling nfmpi_open()
      NF90_EINTOVERFLOW             = NF_EINTOVERFLOW                 ! Overflow when type cast to 4-byte integer



