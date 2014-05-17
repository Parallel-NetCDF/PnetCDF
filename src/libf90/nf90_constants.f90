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
    nf90_ubyte  = 7,            &
    nf90_ushort = 8,            &
    nf90_uint   = 9,            &
    nf90_int64  = 10,           &
    nf90_uint64 = 11
                        
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

  ! netCDF-4 error codes (copied from netCDF release)
  integer, parameter, public :: &
    NF90_EHDFERR      = -101, & ! Error at HDF5 layer. 
    NF90_ECANTREAD    = -102, & ! Can't read. 
    NF90_ECANTWRITE   = -103, & ! Can't write. 
    NF90_ECANTCREATE  = -104, & ! Can't create. 
    NF90_EFILEMETA    = -105, & ! Problem with file metadata. 
    NF90_EDIMMETA     = -106, & ! Problem with dimension metadata. 
    NF90_EATTMETA     = -107, & ! Problem with attribute metadata. 
    NF90_EVARMETA     = -108, & ! Problem with variable metadata. 
    NF90_ENOCOMPOUND  = -109, & ! Not a compound type. 
    NF90_EATTEXISTS   = -110, & ! Attribute already exists. 
    NF90_ENOTNC4      = -111, & ! Attempting netcdf-4 operation on netcdf-3 file.   
    NF90_ESTRICTNC3   = -112, & ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
    NF90_ENOTNC3      = -113, & ! Attempting netcdf-3 operation on netcdf-4 file.   
    NF90_ENOPAR       = -114, & ! Parallel operation on file opened for non-parallel access.   
    NF90_EPARINIT     = -115, & ! Error initializing for parallel access.   
    NF90_EBADGRPID    = -116, & ! Bad group ID.   
    NF90_EBADTYPID    = -117, & ! Bad type ID.   
    NF90_ETYPDEFINED  = -118, & ! Type has already been defined and may not be edited. 
    NF90_EBADFIELD    = -119, & ! Bad field ID.   
    NF90_EBADCLASS    = -120, & ! Bad class.   
    NF90_EMAPTYPE     = -121, & ! Mapped access for atomic types only.   
    NF90_ELATEFILL    = -122, & ! Attempt to define fill value when data already exists. 
    NF90_ELATEDEF     = -123, & ! Attempt to define var properties, like deflate, after enddef.
    NF90_EDIMSCALE    = -124, & ! Probem with HDF5 dimscales.
    NF90_ENOGRP       = -125, & ! No group found.
    NF90_ESTORAGE     = -126, & ! Can't specify both contiguous and chunking.
    NF90_EBADCHUNK    = -127, & ! Bad chunksize.
    NF90_ENOTBUILT    = -128, & ! Attempt to use feature that was not turned on when netCDF was built.
    NF90_EDISKLESS    = -129    ! Error in using diskless  access.

  ! This is the position of NC_NETCDF4 in cmode, counting from the
  ! right, starting (uncharacteristically for fortran) at 0. It's needed
  ! for the BTEST function calls.
  integer, parameter, private :: NETCDF4_BIT = 12


  ! PnetCDF error codes start here
  integer, parameter, public :: &
      NF90_ESMALL                   = NF_ESMALL                   , & ! size of off_t too small for format
      NF90_ENOTINDEP                = NF_ENOTINDEP                , & ! Operation not allowed in collective data mode
      NF90_EINDEP                   = NF_EINDEP                   , & ! Operation not allowed in independent data mode
      NF90_EFILE                    = NF_EFILE                    , & ! Unknown error in file operation
      NF90_EREAD                    = NF_EREAD                    , & ! Unknown error in reading file
      NF90_EWRITE                   = NF_EWRITE                   , & ! Unknown error in writting to file
      NF90_EOFILE                   = NF_EOFILE                   , & ! file open/creation failed
      NF90_EMULTITYPES              = NF_EMULTITYPES              , & ! Multiple types used in memory data
      NF90_EIOMISMATCH              = NF_EIOMISMATCH              , & ! Input/Output data amount mismatch
      NF90_ENEGATIVECNT             = NF_ENEGATIVECNT             , & ! Negative count is specified
      NF90_EUNSPTETYPE              = NF_EUNSPTETYPE              , & ! Unsupported etype in memory MPI datatype
      NF90_EINVAL_REQUEST           = NF_EINVAL_REQUEST           , & ! invalid nonblocking request ID
      NF90_EAINT_TOO_SMALL          = NF_EAINT_TOO_SMALL          , & ! MPI_Aint not large enough to hold requested value
      NF90_ENOTSUPPORT              = NF_ENOTSUPPORT              , & ! feature is not yet supported
      NF90_ENULLBUF                 = NF_ENULLBUF                 , & ! trying to attach a NULL buffer
      NF90_EPREVATTACHBUF           = NF_EPREVATTACHBUF           , & ! previous attached buffer is found
      NF90_ENULLABUF                = NF_ENULLABUF                , & ! no attached buffer is found
      NF90_EPENDINGBPUT             = NF_EPENDINGBPUT             , & ! pending bput is found, cannot detach buffer
      NF90_EINSUFFBUF               = NF_EINSUFFBUF               , & ! attached buffer is too small
      NF90_ENOENT                   = NF_ENOENT                   , & ! File does not exist when calling nfmpi_open()
      NF90_EINTOVERFLOW             = NF_EINTOVERFLOW             , & ! Overflow when type cast to 4-byte integer
      NF90_ENOTENABLED              = NF_ENOTENABLED              , & ! Overflow when type cast to 4-byte integer
      NF90_EBAD_FILE                = NF_EBAD_FILE                , & ! Invalid file name (e.g., path name too long)
      NF90_ENO_SPACE                = NF_ENO_SPACE                , & ! Not enough space
      NF90_EQUOTA                   = NF_EQUOTA                   , & ! Quota exceeded
      NF90_EMULTIDEFINE             = NF_EMULTIDEFINE             , & ! NC definitions on multiprocesses conflict
      NF90_EMULTIDEFINE_OMODE       = NF_EMULTIDEFINE_OMODE       , & ! file create/open modes are inconsistent among processes
      NF90_EMULTIDEFINE_DIM_NUM     = NF_EMULTIDEFINE_DIM_NUM     , & ! inconsistent number of dimensions
      NF90_EMULTIDEFINE_DIM_SIZE    = NF_EMULTIDEFINE_DIM_SIZE    , & ! inconsistent size of dimension
      NF90_EMULTIDEFINE_DIM_NAME    = NF_EMULTIDEFINE_DIM_NAME    , & ! inconsistent dimension names
      NF90_EMULTIDEFINE_VAR_NUM     = NF_EMULTIDEFINE_VAR_NUM     , & ! inconsistent number of variables
      NF90_EMULTIDEFINE_VAR_NAME    = NF_EMULTIDEFINE_VAR_NAME    , & ! inconsistent variable name
      NF90_EMULTIDEFINE_VAR_NDIMS   = NF_EMULTIDEFINE_VAR_NDIMS   , & ! inconsistent variable's number of dimensions
      NF90_EMULTIDEFINE_VAR_DIMIDS  = NF_EMULTIDEFINE_VAR_DIMIDS  , & ! inconsistent variable's dimid
      NF90_EMULTIDEFINE_VAR_TYPE    = NF_EMULTIDEFINE_VAR_TYPE    , & ! inconsistent variable's data type
      NF90_EMULTIDEFINE_VAR_LEN     = NF_EMULTIDEFINE_VAR_LEN     , & ! inconsistent variable's size
      NF90_EMULTIDEFINE_NUMRECS     = NF_EMULTIDEFINE_NUMRECS     , & ! inconsistent number of records
      NF90_EMULTIDEFINE_VAR_BEGIN   = NF_EMULTIDEFINE_VAR_BEGIN   , & ! inconsistent variable file begin offset (internal use)
      NF90_EMULTIDEFINE_ATTR_NUM    = NF_EMULTIDEFINE_ATTR_NUM    , & ! inconsistent number of attributes
      NF90_EMULTIDEFINE_ATTR_SIZE   = NF_EMULTIDEFINE_ATTR_SIZE   , & ! inconsistent memory space used by attribute (internal use)
      NF90_EMULTIDEFINE_ATTR_NAME   = NF_EMULTIDEFINE_ATTR_NAME   , & ! inconsistent attribute name
      NF90_EMULTIDEFINE_ATTR_TYPE   = NF_EMULTIDEFINE_ATTR_TYPE   , & ! inconsistent attribute type
      NF90_EMULTIDEFINE_ATTR_LEN    = NF_EMULTIDEFINE_ATTR_LEN    , & ! inconsistent attribute length
      NF90_EMULTIDEFINE_ATTR_VAL    = NF_EMULTIDEFINE_ATTR_VAL    , & ! inconsistent attribute value
      NF90_EMULTIDEFINE_FNC_ARGS    = NF_EMULTIDEFINE_FNC_ARGS    , & ! inconsistent function arguments used in collective API
      NF90_ECMODE                   = NF_EMULTIDEFINE_OMODE       , &
      NF90_EDIMS_NELEMS_MULTIDEFINE = NF_EMULTIDEFINE_DIM_NUM     , &
      NF90_EDIMS_SIZE_MULTIDEFINE   = NF_EMULTIDEFINE_DIM_SIZE    , &
      NF90_EDIMS_NAME_MULTIDEFINE   = NF_EMULTIDEFINE_DIM_NAME    , &
      NF90_EVARS_NELEMS_MULTIDEFINE = NF_EMULTIDEFINE_VAR_NUM     , &
      NF90_EVARS_NAME_MULTIDEFINE   = NF_EMULTIDEFINE_VAR_NAME    , &
      NF90_EVARS_NDIMS_MULTIDEFINE  = NF_EMULTIDEFINE_VAR_NDIMS   , &
      NF90_EVARS_DIMIDS_MULTIDEFINE = NF_EMULTIDEFINE_VAR_DIMIDS  , &
      NF90_EVARS_TYPE_MULTIDEFINE   = NF_EMULTIDEFINE_VAR_TYPE    , &
      NF90_EVARS_LEN_MULTIDEFINE    = NF_EMULTIDEFINE_VAR_LEN     , &
      NF90_ENUMRECS_MULTIDEFINE     = NF_EMULTIDEFINE_NUMRECS     , &
      NF90_EVARS_BEGIN_MULTIDEFINE  = NF_EMULTIDEFINE_VAR_BEGIN
