!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$
!
! This file lists the Fortran constants used by PnetCDF
!
!

!
! external netcdf data types:
!
      integer, parameter, public :: &
      nf_byte = 1, &
      nf_int1 = nf_byte, &
      nf_char = 2, &
      nf_short = 3, &
      nf_int2 = nf_short, &
      nf_int = 4, &
      nf_float = 5, &
      nf_real = nf_float, &
      nf_double = 6, &
      nf_int64 = 10

!
! default fill values:
!
      integer, parameter, public :: &
      nf_fill_byte = -127, &
      nf_fill_int1 = nf_fill_byte, &
      nf_fill_char = 0, &
      nf_fill_short = -32767, &
      nf_fill_int2 = nf_fill_short, &
      nf_fill_int = -2147483647

      real, parameter, public :: &
      nf_fill_float = 9.9692099683868690e+36, &
      nf_fill_real = nf_fill_float

      double precision, parameter, public :: &
      nf_fill_double = 9.9692099683868690e+36

      integer*8, parameter, public :: &
      nf_fill_int64 = -9223372036854775806_8

!
! mode flags for opening and creating a netcdf dataset:
!
      integer, parameter, public :: &
      nf_nowrite = 0, &
      nf_write = 1, &
      nf_clobber = 0, &
      nf_noclobber = 4, &
      nf_fill = 0, &
      nf_nofill = 256, &
      nf_lock = 1024, &
      nf_share = 2048, &
      nf_64bit_offset = 512, &
      nf_64bit_data = 16, &
      nf_32bit = 16777216, &
      nf_sizehint_default = 0, &
      nf_align_chunk = -1, &
      nf_format_classic = 1, &
      nf_format_64bit = 2, &
      nf_format_64bit_data = 5

!
! size argument for defining an unlimited dimension:
!
      integer, parameter, public :: &
      nf_unlimited = 0

!
! global attribute id:
!
      integer, parameter, public :: &
      nf_global = 0

!
! implementation limits:
!
      integer, parameter, public :: &
      nf_max_dims = 1024, &
      nf_max_attrs = 8192, &
      nf_max_vars = 8192, &
      nf_max_name = 256, &
      nf_max_var_dims = nf_max_dims

!
! error codes:
!
      integer, parameter, public :: &
      NF_NOERR        =  0, &  ! No Error
      NF2_ERR         = -1, &  ! Returned for all errors in the v2 API
      NF_EBADID       = -33, & ! Not a netcdf id
      NF_ENFILE       = -34, & ! Too many netcdfs open
      NF_EEXIST       = -35, & ! netcdf file exists and NF_NOCLOBBER
      NF_EINVAL       = -36, & ! Invalid Argument
      NF_EPERM        = -37, & ! Write to read only
      NF_ENOTINDEFINE = -38, & ! Operation not allowed in data mode
      NF_EINDEFINE    = -39, & ! Operation not allowed in define mode
      NF_EINVALCOORDS = -40, & ! Index exceeds dimension bound
      NF_EMAXDIMS     = -41, & ! NF_MAX_DIMS exceeded
      NF_ENAMEINUSE   = -42, & ! String match to name in use
      NF_ENOTATT      = -43, & ! Attribute not found
      NF_EMAXATTS     = -44, & ! NF_MAX_ATTRS exceeded
      NF_EBADTYPE     = -45, & ! Not a netcdf data type
      NF_EBADDIM      = -46, & ! Invalid dimension id or name
      NF_EUNLIMPOS    = -47, & ! NF_UNLIMITED in the wrong index
      NF_EMAXVARS     = -48, & ! NF_MAX_VARS exceeded
      NF_ENOTVAR      = -49, & ! Variable not found
      NF_EGLOBAL      = -50, & ! Action prohibited on NF_GLOBAL varid
      NF_ENOTNC       = -51, & ! Not a netcdf file
      NF_ESTS         = -52, & ! In Fortran, string too short
      NF_EMAXNAME     = -53, & ! NF_MAX_NAME exceeded
      NF_EUNLIMIT     = -54, & ! NF_UNLIMITED size already in use
      NF_ENORECVARS   = -55, & ! nc_rec op when there are no record vars
      NF_ECHAR        = -56, & ! Attempt to convert between text & numbers
      NF_EEDGE        = -57, & ! Edge+start exceeds dimension bound
      NF_ESTRIDE      = -58, & ! Illegal stride
      NF_EBADNAME     = -59, & ! Attribute or variable name contains illegal characters
      NF_ERANGE       = -60, & ! Math result not representable
      NF_ENOMEM       = -61, & ! Memory allocation (malloc) failure
      NF_EVARSIZE     = -62, & ! One or more variable sizes violate format constraints
      NF_EDIMSIZE     = -63    ! Invalid dimension size

! PnetCDF error codes start here
      integer, parameter, public :: &
      NF_ESMALL                   = -201, & ! size of off_t too small for format
      NF_ENOTINDEP                = -202, & ! Operation not allowed in collective data mode
      NF_EINDEP                   = -203, & ! Operation not allowed in independent data mode
      NF_EFILE                    = -204, & ! Unknown error in file operation
      NF_EREAD                    = -205, & ! Unknown error in reading file
      NF_EWRITE                   = -206, & ! Unknown error in writting to file
      NF_EMULTIDEFINE             = -207, & ! NC definitions on multiprocesses conflict
      NF_EOFILE                   = -208, & ! file open/creation failed
      NF_EMULTITYPES              = -209, & ! Multiple types used in memory data
      NF_EIOMISMATCH              = -210, & ! Input/Output data amount mismatch
      NF_ENEGATIVECNT             = -211, & ! Negative count is specified
      NF_EUNSPTETYPE              = -212, & ! Unsupported etype in memory MPI datatype
      NF_EDIMS_NELEMS_MULTIDEFINE = -213, & ! Different number of dim defines on multiprocesses conflict
      NF_EDIMS_SIZE_MULTIDEFINE   = -214, & ! Different size of dim defines on multiprocesses conflict
      NF_EVARS_NELEMS_MULTIDEFINE = -215, & ! Different number of var defines on multiprocesses conflict
      NF_EVARS_NDIMS_MULTIDEFINE  = -216, & ! Different dim number of var defines on multiprocesses conflict
      NF_EVARS_DIMIDS_MULTIDEFINE = -217, & ! Different dimid defines on multiprocesses conflict
      NF_EVARS_TYPE_MULTIDEFINE   = -218, & ! Different type of var defines on multiprocesses conflict
      NF_EVARS_LEN_MULTIDEFINE    = -219, & ! Different var lenght defines size on multiprocesses conflict
      NF_EVARS_BEGIN_MULTIDEFINE  = -220, & ! Different var begin defines size on multiprocesses conflict
      NF_ENUMRECS_MULTIDEFINE     = -221, & ! Different number records on multiprocesses conflict
      NF_EINVAL_REQUEST           = -222, & ! invalid nonblocking request ID
      NF_EAINT_TOO_SMALL          = -223, & ! MPI_Aint not large enough to hold requested value
      NF_ECMODE                   = -224, & ! file create modes are inconsistent among processes
      NF_ENOTSUPPORT              = -225, & ! feature is not yet supported
      NF_ENULLBUF                 = -226, & ! trying to attach a NULL buffer
      NF_EPREVATTACHBUF           = -227, & ! previous attached buffer is found
      NF_ENULLABUF                = -228, & ! no attached buffer is found
      NF_EPENDINGBPUT             = -229, & ! pending bput is found, cannot detach buffer
      NF_EINSUFFBUF               = -230, & ! attached buffer is too small
      NF_ENOENT                   = -231, & ! File does not exist when calling nfmpi_open()
      NF_EINTOVERFLOW             = -232    ! Overflow when type cast to 4-byte integer


! error handling modes:
!
      integer, parameter, public :: &
      nf_fatal = 1, &
      nf_verbose = 2

! now we can define nfmpi_unlimited properly
      integer(KIND=MPI_OFFSET_KIND), parameter, public :: &
      nfmpi_unlimited = 0

! NULL request for non-blocking I/O APIs
      integer, parameter, public :: &
      NF_REQ_NULL = -1


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!
! functions in the fortran interface
!

!
! netcdf data types:
!
      integer, parameter, public :: &
      ncbyte = 1, &
      ncchar = 2, &
      ncshort = 3, &
      nclong = 4, &
      ncfloat = 5, &
      ncdouble = 6

!
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!
      integer, parameter, public :: &
      ncrdwr = 1, &     ! read/write, 0 => readonly
      nccreat = 2, &    ! in create phase, cleared by ncendef
      ncexcl = 4, &     ! on create destroy existing file
      ncindef = 8, &    ! in define mode, cleared by ncendef
      ncnsync = 16, &   ! synchronise numrecs on change (x'10')
      nchsync = 32, &   ! synchronise whole header on change (x'20')
      ncndirty = 64, &  ! numrecs has changed (x'40')
      nchdirty = 128, & ! header info has changed (x'80')
      ncfill = 0, &     ! prefill vars on endef and increase of record, the default behavior
      ncnofill = 256, & ! do not fill vars on endef and increase of record (x'100')
      nclink = 32768    ! isa link (x'8000')

!
!     'mode' arguments for nccreate and ncopen
!
      integer, parameter, public :: &
      ncnowrit = 0, &
      ncwrite = ncrdwr, &
      ncclob = nf_clobber, &
      ncnoclob = nf_noclobber

!
!     'size' argument to ncdimdef for an unlimited dimension
!
      integer, parameter, public :: &
      ncunlim = 0

!
!     attribute id to put/get a global attribute
!
      integer, parameter, public :: &
      ncglobal  = 0

!
!     advisory maximums:
!
      integer, parameter, public :: &
      maxncop = 32, &
      maxncdim = 100, &
      maxncatt = 2000, &
      maxncvar = 2000

!     not enforced
      integer, parameter, public :: &
      maxncnam = 128, &
      maxvdims = maxncdim

!
!     global netcdf error status variable
!     initialized in error.c
!
      integer, parameter, public :: &
      ncnoerr = nf_noerr, &         ! no error
      ncebadid = nf_ebadid, &       ! not a netcdf id
      ncenfile = -31, &             ! nc_syserr too many netcdfs open
      nceexist = nf_eexist, &       ! netcdf file exists && ncnoclob
      nceinval = nf_einval, &       ! invalid argument
      nceperm = nf_eperm, &         ! write to read only
      ncenotin = nf_enotindefine, & ! operation not allowed in data mode
      nceindef = nf_eindefine, &    ! operation not allowed in define mode
      ncecoord = nf_einvalcoords, & ! coordinates out of domain
      ncemaxds = nf_emaxdims, &     ! maxncdims exceeded
      ncename = nf_enameinuse, &    ! string match to name in use
      ncenoatt = nf_enotatt, &      ! attribute not found
      ncemaxat = nf_emaxatts, &     ! maxncattrs exceeded
      ncebadty = nf_ebadtype, &     ! not a netcdf data type
      ncebadd = nf_ebaddim, &       ! invalid dimension id
      nceunlim = nf_eunlimpos, &    ! ncunlimited in the wrong index
      ncemaxvs = nf_emaxvars, &     ! maxncvars exceeded
      ncenotvr = nf_enotvar, &      ! variable not found
      nceglob = nf_eglobal, &       ! action prohibited on ncglobal varid
      ncenotnc = nf_enotnc, &       ! not a netcdf file
      ncests = nf_ests, &
      ncentool = nf_emaxname, &
      ncfoobar = 32, &
      ncsyserr = -31

!
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!
      integer, parameter, public :: &
      ncfatal = 1, &
      ncverbos = 2

!
!     default fill values.  these must be the same as in the c interface.
!
      integer, parameter, public :: &
      filbyte = -127, &
      filchar = 0, &
      filshort = -32767, &
      fillong = -2147483647

      real, parameter, public :: &
      filfloat = 9.9692099683868690e+36

      double precision, parameter, public :: &
      fildoub = 9.9692099683868690e+36

