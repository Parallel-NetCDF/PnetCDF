/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STRERROR
#include <string.h> /* contains prototype for ansi libc function strerror() */
#else
/* provide a strerror function for older unix systems */
static char *
strerror(int errnum)
{
    extern int sys_nerr;
    extern char *sys_errlist[];

    if(errnum < 0 || errnum >= sys_nerr) return NULL;
    /* else */
    return sys_errlist[errnum];
}
#endif /* HAVE_STRERROR */

#include <pnetcdf.h>

#ifdef vms
/* UNTESTED */
/*
 * On the vms system, when a system error occurs which is not
 * mapped into the unix styled errno values, errno is set EVMSERR
 * and a VMS error code is set in vaxc$errno.
 * This routine prints the systems message associated with status return
 * from a system services call.
 */

#include <errno.h>
#include <descrip.h>
#include <ssdef.h>

static const char *
vms_strerror( int status )
{
    short msglen;
    static char msgbuf[256];
    $DESCRIPTOR(message, msgbuf);
    register ret;

    msgbuf[0] = 0;
    ret = SYS$GETMSG(status, &msglen, &message, 15, 0);

    if(ret != SS$_BUFFEROVF && ret != SS$_NORMAL) {
        (void) strcpy(msgbuf, "EVMSERR");
    }
    return(msgbuf);
}
#endif /* vms */


/* must copy nc_strerror() from netCDF release */
static const char* nc_strerror(int ncerr1);

static char nc_unknown_err_msg[128];

const char *
ncmpi_strerror(int err)
{
    sprintf(nc_unknown_err_msg,"Unknown Error: Unrecognized error code %5d\n",err);

#ifdef vms
    if (err == EVMSERR) return vms_strerror(err);
#endif /* vms */

    if (NC_ISSYSERR(err)) {
        const char *cp = (const char *) strerror(err);
        if (cp == NULL) return nc_unknown_err_msg;
        return cp;
    }

    switch (err) {
        /* PnetCDF errors */
        case NC_ESMALL:
            return "Size of MPI_Offset or MPI_Aint too small for requested format";
            /* this usually happens on 32-bit machines where MPI_Offset and
             * MPI_Aint may be 4-byte integers and when the I/O request amount
             * or accessing file offset is > 2GB.
             */
        case NC_ENOTINDEP:
            return "Operation not allowed in collective data mode";
            /* this means your program is now in independent data mode, but
             * making a call to a collective API
             */
        case NC_EINDEP:
            return "Operation not allowed in independent data mode";
            /* this means your program is now in collective data mode, but
             * making a call to an independent API
             */
        case NC_EFILE:
            return "Unknown error in file operation";
            /* this error is caused by an unsuccessful MPI-IO call and usually
             * accompany with additional MPI error messages
             */
        case NC_EREAD:
            return "Unknown error occurs in reading file";
            /* this error is caused by an unsuccessful call to MPI_File_read or
             * MPI_File_read_all and usually accompany with additional MPI
             * error messages
             */
        case NC_EWRITE:
            return "Unknown error occurs in writing file";
            /* this error is caused by an unsuccessful call to MPI_File_write or
             * MPI_File_write_all and usually accompany with additional MPI
             * error messages
             */
        case NC_EOFILE:
            return "Fail to open/create file";
        case NC_EMULTITYPES:
            return "Multiple etypes used in MPI datatype";
            /* when using flexible APIs, the argument MPI derived datatype is
             * not allowed to contain more than one basic data type
             */
        case NC_EIOMISMATCH:
            return "Input/Output data amount mismatch";
            /* this error indicates the request amount is mismatched between
             * bufcount and the value calculated from argument count[]
             */
        case NC_ENEGATIVECNT:
            return "Negative count is prohibited";
            /* In netCDF, count (or edge) argument is of type size_t, which is
             * an unsigned integer. A negative count will make the value very
             * large, causing NC_EEDGE error and hence netCDF need no such
             * error code.
             */
        case NC_EUNSPTETYPE:
            return "Unsupported etype in the MPI datatype describing the I/O buffer";
            /* when using flexible APIs, the argument MPI derived datatype is
             * only allowed to be constructed from the MPI basic data types
             * known to PnetCDF
             */
        case NC_EINVAL_REQUEST:
            return "Invalid nonblocking request ID.";
        case NC_EAINT_TOO_SMALL:
            return "MPI_Aint not large enough to hold requested value.";
            /* this usually happens on 32-bit machines where MPI_Aint is a
             * 4-byte integer and when the I/O request amount or accessing
             * file offset is > 2GB.
             */
        case NC_ENOTSUPPORT:
            return "Feature is not yet supported.";
        case NC_ENULLBUF:
            return "Trying to attach a NULL buffer or the buffer size is negative.";
            /* an error returned from ncmpi_buffer_attach()
             */
        case NC_EPREVATTACHBUF:
            return "Previous attached buffer is found.";
            /* an error returned from ncmpi_buffer_attach() indicating a
             * buffer has been attached previously
             */
        case NC_ENULLABUF:
            return "No attached buffer is found.";
            /* an error when calling bput APIs and no buffer has been attached
             */
        case NC_EPENDINGBPUT:
            return "Cannot detach buffer due to pending bput request is found.";
            /* an error returned from ncmpi_buffer_detach()
             */
        case NC_EINSUFFBUF:
            return "Attached buffer is too small.";
            /* an error when calling bput APIs
             */
        case NC_ENOENT:
            return "Specified netCDF file does not exist.";
            /* this error code corresponds to MPI error class
             * MPI_ERR_NO_SUCH_FILE, an error generated from MPI_File_open(),
             * MPI_File_delete() or others
             */
        case NC_EINTOVERFLOW:
            return "Integer type casting overflow.";
            /* this usually happens on 32-bit machines where MPI_Offset is a
             * 4-byte integer and when the I/O request amount or accessing
             * file offset is > 2GB.
             */
        case NC_ENOTENABLED:
            return "Feature is not enabled at configure time.";
            /* Some APIs require a specific feature enabled at the configure
             * time, for example, ncmpi_inq_malloc_size() works only
             * --enable-debug is used when configuring PnetCDF
             */
        case NC_EBAD_FILE:
            return "Invalid file name (e.g., path name too long).";
            /* this error code corresponds to MPI error class
             * MPI_ERR_BAD_FILE, an error generated from MPI_File_open()
             */
        case NC_ENO_SPACE:
            return "Not enough space.";
            /* this error code corresponds to MPI error class
             * MPI_ERR_NO_SPACE, an error generated from MPI_File_open()
             */
        case NC_EQUOTA:
            return "Quota exceeded.";
            /* this error code corresponds to MPI error class
             * MPI_ERR_QUOTA, an error generated from MPI_File_open()
             */
        case NC_ENULLSTART:
            return "Argument start is a NULL pointer";
            /* Some APIs require argument start not be a NULL pointer
             */
        case NC_ENULLCOUNT:
            return "Argument count is a NULL pointer";
            /* Some APIs require argument count cannot be a NULL pointer
             */
        case NC_EINVAL_CMODE:
            return "Invalid file create mode";
        case NC_ETYPESIZE:
            return "MPI datatype size error (bigger than the variable size)";
        case NC_ETYPE_MISMATCH:
            return "etype of the MPI datatype mismatches the variable data type";
        case NC_ETYPESIZE_MISMATCH:
            return "MPI filetype size mismatches buftype size * bufcount";
        case NC_ESTRICTCDF2:
            return "Attempting CDF-5 operation on strict CDF or CDF-2 file";
        case NC_ENOTRECVAR:
            return "Attempting operation only for record variables";
        case NC_ENOTFILL:
            return "Attempting to fill a record when its variable fill mode is off";
        case NC_EINVAL_OMODE:
            return "Invalid or unsupported file open mode";
        case NC_EPENDING:
            return "Pending nonblocking request is found at file close";
        case NC_EMAX_REQ:
            return "Size of I/O request exceeds INT_MAX";
        case NC_EMULTIDEFINE:
            return "File header is inconsistent among processes";
            /* this error means the metadata (dimension names, variable names,
             * variable's dimensions, attributes, and whatever will be stored
             * in the file header) is inconsistent among all MPI processes.
             */
        case NC_EMULTIDEFINE_OMODE:
            return "File open mode is inconsistent among processes.";
        case NC_EMULTIDEFINE_DIM_NUM:
            return "Number of dimensions is inconsistent among processes.";
        case NC_EMULTIDEFINE_DIM_SIZE:
            return "Dimension size is inconsistent among processes.";
        case NC_EMULTIDEFINE_DIM_NAME:
            return "Dimension name is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_NUM:
            return "Number of variables is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_NAME:
            return "Variable name is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_NDIMS:
            return "Dimensionality of this variable is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_DIMIDS:
            return "Dimension IDs used to define this variable are inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_TYPE:
            return "Data type of this variable is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_LEN:
            return "Number of elements of this variable is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_BEGIN:
            return "Starting file offset of this variable is inconsistent among processes.";
        case NC_EMULTIDEFINE_NUMRECS:
            return "Number of records is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_NUM:
            return "Number of attributes is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_SIZE:
            return "Memory space used by attribute (internal use) is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_NAME:
            return "Attribute name is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_TYPE:
            return "Attribute type is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_LEN:
            return "Attribute length is inconsistent among processes.";
        case NC_EMULTIDEFINE_ATTR_VAL:
            return "Attribute value is inconsistent among processes.";
        case NC_EMULTIDEFINE_FNC_ARGS:
            return "Arguments in collective API are inconsistent among processes.";
        case NC_EMULTIDEFINE_FILL_MODE:
            return "File fill mode is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_FILL_MODE:
            return "Variable fill mode is inconsistent among processes.";
        case NC_EMULTIDEFINE_VAR_FILL_VALUE:
            return "Variable fill value is inconsistent among processes.";
        case NC_EMULTIDEFINE_CMODE:
            return "File create mode is inconsistent among processes.";
        case NC_EBADLOG:
            return "Unrecognized burst buffering log file format.";
        case NC_EFLUSHED:
            return "Nonblocking requests already flushed.";
        case NC_EADIOS:
            return "unknown ADIOS error.";

        default:
            /* check netCDF-3 and netCDF-4 errors */
            return nc_strerror(err);
    }
}

/*-----------------------------------------------------------------------------
 *     Copy nc_strerror() from netCDF release to ensure the same error
 *     message is returned.
 *----------------------------------------------------------------------------*/

/*! NetCDF Error Handling

\addtogroup error NetCDF Error Handling

NetCDF functions return a non-zero status codes on error.

Each netCDF function returns an integer status value. If the returned
status value indicates an error, you may handle it in any way desired,
from printing an associated error message and exiting to ignoring the
error indication and proceeding (not recommended!). For simplicity,
the examples in this guide check the error status and call a separate
function, handle_err(), to handle any errors. One possible definition
of handle_err() can be found within the documentation of
nc_strerror().

The nc_strerror() function is available to convert a returned integer
error status into an error message string.

Occasionally, low-level I/O errors may occur in a layer below the
netCDF library. For example, if a write operation causes you to exceed
disk quotas or to attempt to write to a device that is no longer
available, you may get an error from a layer below the netCDF library,
but the resulting write error will still be reflected in the returned
status value.

*/

/** \{ */

/*! Given an error number, return an error message.

This function returns a static reference to an error message string
corresponding to an integer netCDF error status or to a system error
number, presumably returned by a previous call to some other netCDF
function. The error codes are defined in netcdf.h.

\param ncerr1 error number

\returns short string containing error message.

Here is an example of a simple error handling function that uses
nc_strerror() to print the error message corresponding to the netCDF
error status returned from any netCDF function call and then exit:

\code
     #include <netcdf.h>
        ...
     void handle_error(int status) {
     if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(-1);
        }
     }
\endcode
*/

static const char *
nc_strerror(int ncerr1)
{
   /* System error? */
   if(NC_ISSYSERR(ncerr1))
   {
      const char *cp = (const char *) strerror(ncerr1);
      if(cp == NULL)
	 return "Unknown Error";
      return cp;
   }

   /* If we're here, this is a netcdf error code. */
   switch(ncerr1)
   {
      case NC_NOERR:
	 return "No error";
      case NC_EBADID:
	 return "NetCDF: Not a valid ID";
      case NC_ENFILE:
	 return "NetCDF: Too many files open";
      case NC_EEXIST:
	 return "NetCDF: File exists && NC_NOCLOBBER";
      case NC_EINVAL:
	 return "NetCDF: Invalid argument";
      case NC_EPERM:
	 return "NetCDF: Write to read only";
      case NC_ENOTINDEFINE:
	 return "NetCDF: Operation not allowed in data mode";
      case NC_EINDEFINE:
	 return "NetCDF: Operation not allowed in define mode";
      case NC_EINVALCOORDS:
	 return "NetCDF: Index exceeds dimension bound";
      case NC_EMAXDIMS:
	 return "NetCDF: NC_MAX_DIMS or NC_MAX_VAR_DIMS exceeded"; /* not enforced after 4.5.0 */
      case NC_ENAMEINUSE:
	 return "NetCDF: String match to name in use";
      case NC_ENOTATT:
	 return "NetCDF: Attribute not found";
      case NC_EMAXATTS:
	 return "NetCDF: NC_MAX_ATTRS exceeded"; /* not enforced after 4.5.0 */
      case NC_EBADTYPE:
	 return "NetCDF: Not a valid data type or _FillValue type mismatch";
      case NC_EBADDIM:
	 return "NetCDF: Invalid dimension ID or name";
      case NC_EUNLIMPOS:
	 return "NetCDF: NC_UNLIMITED in the wrong index";
      case NC_EMAXVARS:
	 return "NetCDF: NC_MAX_VARS exceeded"; /* not enforced after 4.5.0 */
      case NC_ENOTVAR:
	 return "NetCDF: Variable not found";
      case NC_EGLOBAL:
	 return "NetCDF: Action prohibited on NC_GLOBAL varid";
      case NC_ENOTNC:
	 return "NetCDF: Unknown file format (file format violates CDF specification)";
      case NC_ESTS:
	 return "NetCDF: In Fortran, string too short";
      case NC_EMAXNAME:
	 return "NetCDF: NC_MAX_NAME exceeded";
      case NC_EUNLIMIT:
	 return "NetCDF: NC_UNLIMITED size already in use";
      case NC_ENORECVARS:
	 return "NetCDF: nc_rec op when there are no record vars";
      case NC_ECHAR:
	 return "NetCDF: Attempt to convert between text & numbers";
      case NC_EEDGE:
	 return "NetCDF: Start+count exceeds dimension bound";
      case NC_ESTRIDE:
	 return "NetCDF: Illegal stride";
      case NC_EBADNAME:
	 return "NetCDF: Name contains illegal characters";
      case NC_ERANGE:
	 return "NetCDF: Numeric conversion not representable";
      case NC_ENOMEM:
	 return "NetCDF: Memory allocation (malloc) failure";
      case NC_EVARSIZE:
	 return "NetCDF: One or more variable sizes violate format constraints";
      case NC_EDIMSIZE:
	 return "NetCDF: Invalid dimension size";
      case NC_ETRUNC:
	 return "NetCDF: File likely truncated or possibly corrupted";
      case NC_EAXISTYPE:
	 return "NetCDF: Illegal axis type";
      case NC_EDAP:
	 return "NetCDF: DAP failure";
      case NC_ECURL:
	 return "NetCDF: libcurl failure";
      case NC_EIO:
	 return "NetCDF: I/O failure";
      case NC_ENODATA:
	 return "NetCDF: Variable has no data";
      case NC_EDAPSVC:
	 return "NetCDF: DAP server error";
      case NC_EDAS:
	 return "NetCDF: Malformed or inaccessible DAP DAS";
      case NC_EDDS:
	 return "NetCDF: Malformed or inaccessible DAP2 DDS or DAP4 DMR response";
      case NC_EDATADDS:
	 return "NetCDF: Malformed or inaccessible DAP2 DATADDS or DAP4 DAP response";
      case NC_EDAPURL:
	 return "NetCDF: Malformed URL";
      case NC_EDAPCONSTRAINT:
	 return "NetCDF: Malformed or unexpected Constraint";
      case NC_ETRANSLATION:
	 return "NetCDF: Untranslatable construct";
      case NC_EACCESS:
	 return "NetCDF: Access failure";
      case NC_EAUTH:
	 return "NetCDF: Authorization failure";
      case NC_ENOTFOUND:
	 return "NetCDF: file not found";
      case NC_ECANTREMOVE:
	 return "NetCDF: cannot delete file";
      case NC_EINTERNAL:
	 return "NetCDF: internal library error; Please contact Unidata support";
      case NC_EPNETCDF:
	 return "NetCDF: PnetCDF error";
      case NC_EHDFERR:
	 return "NetCDF: HDF error";
      case NC_ECANTREAD:
	 return "NetCDF: Can't read file";
      case NC_ECANTWRITE:
	 return "NetCDF: Can't write file";
      case NC_ECANTCREATE:
	 return "NetCDF: Can't create file";
      case NC_EFILEMETA:
	 return "NetCDF: Can't add HDF5 file metadata";
      case NC_EDIMMETA:
	 return "NetCDF: Can't define dimensional metadata";
      case NC_EATTMETA:
	 return "NetCDF: Can't open HDF5 attribute";
      case NC_EVARMETA:
	 return "NetCDF: Problem with variable metadata.";
      case NC_ENOCOMPOUND:
	 return "NetCDF: Can't create HDF5 compound type";
      case NC_EATTEXISTS:
	 return "NetCDF: Attempt to create attribute that already exists";
      case NC_ENOTNC4:
	 return "NetCDF: Attempting netcdf-4 operation on netcdf-3 file";
      case NC_ESTRICTNC3:
	 return "NetCDF: Attempting netcdf-4 operation on strict nc3 netcdf-4 file";
      case NC_ENOTNC3:
	 return "NetCDF: Attempting netcdf-3 operation on netcdf-4 file";
      case NC_ENOPAR:
	 return "NetCDF: Parallel operation on file opened for non-parallel access";
      case NC_EPARINIT:
	 return "NetCDF: Error initializing for parallel access";
      case NC_EBADGRPID:
	 return "NetCDF: Bad group ID";
      case NC_EBADTYPID:
	 return "NetCDF: Bad type ID";
      case NC_ETYPDEFINED:
	 return "NetCDF: Type has already been defined and may not be edited";
      case NC_EBADFIELD:
	 return "NetCDF: Bad field ID";
      case NC_EBADCLASS:
	 return "NetCDF: Bad class";
      case NC_EMAPTYPE:
	 return "NetCDF: Mapped access for atomic types only";
      case NC_ELATEFILL:
	 return "NetCDF: Attempt to define fill value when data already exists.";
      case NC_ELATEDEF:
	 return "NetCDF: Attempt to define var properties, like deflate, after enddef.";
      case NC_EDIMSCALE:
	 return "NetCDF: Problem with HDF5 dimscales.";
      case NC_ENOGRP:
	 return "NetCDF: No group found.";
      case NC_ESTORAGE:
	 return "NetCDF: Cannot specify both contiguous and chunking.";
      case NC_EBADCHUNK:
	 return "NetCDF: Bad chunk sizes.";
      case NC_ENOTBUILT:
	 return "NetCDF: Attempt to use feature that was not turned on "
	    "when netCDF was built.";
      case NC_EDISKLESS:
	 return "NetCDF: Error in using diskless access";
      case NC_EFILTER:
	 return "NetCDF: Filter error: bad id or parameters or duplicate filter";
      case NC_ENOFILTER:
	 return "NetCDF: Filter error: unimplemented filter encountered";
      case NC_ECANTEXTEND:
	return "NetCDF: Attempt to extend dataset during NC_INDEPENDENT I/O operation. Use nc_var_par_access to set mode NC_COLLECTIVE before extending variable.";
      case NC_EMPI:
	 return "NetCDF: MPI operation failed.";
      case NC_ERCFILE:
	 return "NetCDF: RC File Failure.";
     case NC_ENULLPAD:
       return "NetCDF: File fails strict Null-Byte Header check.";
     case NC_EINMEMORY:
       return "NetCDF: In-memory File operation failed.";
      case NC_ENCZARR:
	 return "NetCDF: NCZarr error";
      case NC_ES3:
	 return "NetCDF: AWS S3 error";
      case NC_EEMPTY:
	 return "NetCDF: Attempt to read empty NCZarr map key";
      case NC_EOBJECT:
	 return "NetCDF: Some object exists when it should not";
      case NC_ENOOBJECT:
	 return "NetCDF: Some object not found";
      case NC_EPLUGIN:
	 return "NetCDF: Unclassified failure in accessing a dynamically loaded plugin";
      default:
	 return nc_unknown_err_msg;
   }
}
/** \} */

/*----< ncmpi_strerrno() >---------------------------------------------------*/
/* print the NC error code name */
const char *
ncmpi_strerrno(int err)
{
    static char unknown_str[64];

    if (err > 0) { /* system error */
        const char *cp = (const char *) strerror(err);
        if (cp == NULL)
            sprintf(unknown_str,"Unknown error code %d",err);
        else
            sprintf(unknown_str,"System error code %d (%s)",err,cp);
        return unknown_str;
    }

    switch (err) {
        case (NC_NOERR):			return "NC_NOERR";
        case (NC_EBADID):			return "NC_EBADID";
        case (NC_ENFILE):			return "NC_ENFILE";
        case (NC_EEXIST):			return "NC_EEXIST";
        case (NC_EINVAL):			return "NC_EINVAL";
        case (NC_EPERM):			return "NC_EPERM";
        case (NC_ENOTINDEFINE):			return "NC_ENOTINDEFINE";
        case (NC_EINDEFINE):			return "NC_EINDEFINE";
        case (NC_EINVALCOORDS):			return "NC_EINVALCOORDS";
        case (NC_EMAXDIMS):			return "NC_EMAXDIMS";
        case (NC_ENAMEINUSE):			return "NC_ENAMEINUSE";
        case (NC_ENOTATT):			return "NC_ENOTATT";
        case (NC_EMAXATTS):			return "NC_EMAXATTS";
        case (NC_EBADTYPE):			return "NC_EBADTYPE";
        case (NC_EBADDIM):			return "NC_EBADDIM";
        case (NC_EUNLIMPOS):			return "NC_EUNLIMPOS";
        case (NC_EMAXVARS):			return "NC_EMAXVARS";
        case (NC_ENOTVAR):			return "NC_ENOTVAR";
        case (NC_EGLOBAL):			return "NC_EGLOBAL";
        case (NC_ENOTNC):			return "NC_ENOTNC";
        case (NC_ESTS):				return "NC_ESTS";
        case (NC_EMAXNAME):			return "NC_EMAXNAME";
        case (NC_EUNLIMIT):			return "NC_EUNLIMIT";
        case (NC_ENORECVARS):			return "NC_ENORECVARS";
        case (NC_ECHAR):			return "NC_ECHAR";
        case (NC_EEDGE):			return "NC_EEDGE";
        case (NC_ESTRIDE):			return "NC_ESTRIDE";
        case (NC_EBADNAME):			return "NC_EBADNAME";
        case (NC_ERANGE):			return "NC_ERANGE";
        case (NC_ENOMEM):			return "NC_ENOMEM";
        case (NC_EVARSIZE):			return "NC_EVARSIZE";
        case (NC_EDIMSIZE):			return "NC_EDIMSIZE";
        case (NC_ETRUNC):			return "NC_ETRUNC";
        case (NC_EAXISTYPE):			return "NC_EAXISTYPE";
        case (NC_EDAP):				return "NC_EDAP";
        case (NC_ECURL):			return "NC_ECURL";
        case (NC_EIO):				return "NC_EIO";
        case (NC_ENODATA):			return "NC_ENODATA";
        case (NC_EDAPSVC):			return "NC_EDAPSVC";
        case (NC_EDAS):				return "NC_EDAS";
        case (NC_EDDS):				return "NC_EDDS";
        case (NC_EDATADDS):			return "NC_EDATADDS";
        case (NC_EDAPURL):			return "NC_EDAPURL";
        case (NC_EDAPCONSTRAINT):		return "NC_EDAPCONSTRAINT";
        case (NC_ETRANSLATION):			return "NC_ETRANSLATION";
        case (NC_EACCESS):			return "NC_EACCESS";
        case (NC_EAUTH):			return "NC_EAUTH";
        case (NC_ENOTFOUND):			return "NC_ENOTFOUND";
        case (NC_ECANTREMOVE):			return "NC_ECANTREMOVE";
        case (NC_EINTERNAL):			return "NC_EINTERNAL";
        case (NC_EPNETCDF):			return "NC_EPNETCDF";
        case (NC_EHDFERR):			return "NC_EHDFERR";
        case (NC_ECANTREAD):			return "NC_ECANTREAD";
        case (NC_ECANTWRITE):			return "NC_ECANTWRITE";
        case (NC_ECANTCREATE):			return "NC_ECANTCREATE";
        case (NC_EFILEMETA):			return "NC_EFILEMETA";
        case (NC_EDIMMETA):			return "NC_EDIMMETA";
        case (NC_EATTMETA):			return "NC_EATTMETA";
        case (NC_EVARMETA):			return "NC_EVARMETA";
        case (NC_ENOCOMPOUND):			return "NC_ENOCOMPOUND";
        case (NC_EATTEXISTS):			return "NC_EATTEXISTS";
        case (NC_ENOTNC4):			return "NC_ENOTNC4";
        case (NC_ESTRICTNC3):			return "NC_ESTRICTNC3";
        case (NC_ENOTNC3):			return "NC_ENOTNC3";
        case (NC_ENOPAR):			return "NC_ENOPAR";
        case (NC_EPARINIT):			return "NC_EPARINIT";
        case (NC_EBADGRPID):			return "NC_EBADGRPID";
        case (NC_EBADTYPID):			return "NC_EBADTYPID";
        case (NC_ETYPDEFINED):			return "NC_ETYPDEFINED";
        case (NC_EBADFIELD):			return "NC_EBADFIELD";
        case (NC_EBADCLASS):			return "NC_EBADCLASS";
        case (NC_EMAPTYPE):			return "NC_EMAPTYPE";
        case (NC_ELATEFILL):			return "NC_ELATEFILL";
        case (NC_ELATEDEF):			return "NC_ELATEDEF";
        case (NC_EDIMSCALE):			return "NC_EDIMSCALE";
        case (NC_ENOGRP):			return "NC_ENOGRP";
        case (NC_ESTORAGE):			return "NC_ESTORAGE";
        case (NC_EBADCHUNK):			return "NC_EBADCHUNK";
        case (NC_ENOTBUILT):			return "NC_ENOTBUILT";
        case (NC_EDISKLESS):			return "NC_EDISKLESS";
        case (NC_ECANTEXTEND):			return "NC_ECANTEXTEND";
        case (NC_EMPI):				return "NC_EMPI";
        case (NC_EFILTER):			return "NC_EFILTER";
        case (NC_ERCFILE):			return "NC_ERCFILE";
        case (NC_ENULLPAD):			return "NC_ENULLPAD";
        case (NC_EINMEMORY):			return "NC_EINMEMORY";
        case (NC_ENOFILTER):			return "NC_ENOFILTER";
/*
        case (NC_EURL):				return "NC_EURL";
        case (NC_ECONSTRAINT):			return "NC_ECONSTRAINT";
*/
        case (NC_ESMALL):			return "NC_ESMALL";
        case (NC_ENOTINDEP):			return "NC_ENOTINDEP";
        case (NC_EINDEP):			return "NC_EINDEP";
        case (NC_EFILE):			return "NC_EFILE";
        case (NC_EREAD):			return "NC_EREAD";
        case (NC_EWRITE):			return "NC_EWRITE";
        case (NC_EOFILE):			return "NC_EOFILE";
        case (NC_EMULTITYPES):			return "NC_EMULTITYPES";
        case (NC_EIOMISMATCH):			return "NC_EIOMISMATCH";
        case (NC_ENEGATIVECNT):			return "NC_ENEGATIVECNT";
        case (NC_EUNSPTETYPE):			return "NC_EUNSPTETYPE";
        case (NC_EINVAL_REQUEST):		return "NC_EINVAL_REQUEST";
        case (NC_EAINT_TOO_SMALL):		return "NC_EAINT_TOO_SMALL";
        case (NC_ENOTSUPPORT):			return "NC_ENOTSUPPORT";
        case (NC_ENULLBUF):			return "NC_ENULLBUF";
        case (NC_EPREVATTACHBUF):		return "NC_EPREVATTACHBUF";
        case (NC_ENULLABUF):			return "NC_ENULLABUF";
        case (NC_EPENDINGBPUT):			return "NC_EPENDINGBPUT";
        case (NC_EINSUFFBUF):			return "NC_EINSUFFBUF";
        case (NC_ENOENT):			return "NC_ENOENT";
        case (NC_EINTOVERFLOW):			return "NC_EINTOVERFLOW";
        case (NC_ENOTENABLED):			return "NC_ENOTENABLED";
        case (NC_EBAD_FILE):			return "NC_EBAD_FILE";
        case (NC_ENO_SPACE):			return "NC_ENO_SPACE";
        case (NC_EQUOTA):			return "NC_EQUOTA";
        case (NC_ENULLSTART):			return "NC_ENULLSTART";
        case (NC_ENULLCOUNT):			return "NC_ENULLCOUNT";
        case (NC_EINVAL_CMODE):			return "NC_EINVAL_CMODE";
        case (NC_ETYPESIZE):			return "NC_ETYPESIZE";
        case (NC_ETYPE_MISMATCH):		return "NC_ETYPE_MISMATCH";
        case (NC_ETYPESIZE_MISMATCH):		return "NC_ETYPESIZE_MISMATCH";
        case (NC_ESTRICTCDF2):			return "NC_ESTRICTCDF2";
        case (NC_ENOTRECVAR):			return "NC_ENOTRECVAR";
        case (NC_ENOTFILL):			return "NC_ENOTFILL";
        case (NC_EINVAL_OMODE):			return "NC_EINVAL_OMODE";
        case (NC_EPENDING):			return "NC_EPENDING";
        case (NC_EMAX_REQ):			return "NC_EMAX_REQ";
        case (NC_EBADLOG):			return "NC_EBADLOG";
        case (NC_EFLUSHED):			return "NC_EFLUSHED";
        case (NC_EADIOS):			return "NC_EADIOS";

        case (NC_EMULTIDEFINE):			return "NC_EMULTIDEFINE";
        case (NC_EMULTIDEFINE_OMODE):		return "NC_EMULTIDEFINE_OMODE";
        case (NC_EMULTIDEFINE_DIM_NUM):		return "NC_EMULTIDEFINE_DIM_NUM";
        case (NC_EMULTIDEFINE_DIM_SIZE):	return "NC_EMULTIDEFINE_DIM_SIZE";
        case (NC_EMULTIDEFINE_DIM_NAME):	return "NC_EMULTIDEFINE_DIM_NAME";
        case (NC_EMULTIDEFINE_VAR_NUM):		return "NC_EMULTIDEFINE_VAR_NUM";
        case (NC_EMULTIDEFINE_VAR_NAME):	return "NC_EMULTIDEFINE_VAR_NAME";
        case (NC_EMULTIDEFINE_VAR_NDIMS):	return "NC_EMULTIDEFINE_VAR_NDIMS";
        case (NC_EMULTIDEFINE_VAR_DIMIDS):	return "NC_EMULTIDEFINE_VAR_DIMIDS";
        case (NC_EMULTIDEFINE_VAR_TYPE):	return "NC_EMULTIDEFINE_VAR_TYPE";
        case (NC_EMULTIDEFINE_VAR_LEN):		return "NC_EMULTIDEFINE_VAR_LEN";
        case (NC_EMULTIDEFINE_NUMRECS):		return "NC_EMULTIDEFINE_NUMRECS";
        case (NC_EMULTIDEFINE_VAR_BEGIN):	return "NC_EMULTIDEFINE_VAR_BEGIN";
        case (NC_EMULTIDEFINE_ATTR_NUM):	return "NC_EMULTIDEFINE_ATTR_NUM";
        case (NC_EMULTIDEFINE_ATTR_SIZE):	return "NC_EMULTIDEFINE_ATTR_SIZE";
        case (NC_EMULTIDEFINE_ATTR_NAME):	return "NC_EMULTIDEFINE_ATTR_NAME";
        case (NC_EMULTIDEFINE_ATTR_TYPE):	return "NC_EMULTIDEFINE_ATTR_TYPE";
        case (NC_EMULTIDEFINE_ATTR_LEN):	return "NC_EMULTIDEFINE_ATTR_LEN";
        case (NC_EMULTIDEFINE_ATTR_VAL):	return "NC_EMULTIDEFINE_ATTR_VAL";
        case (NC_EMULTIDEFINE_FNC_ARGS):	return "NC_EMULTIDEFINE_FNC_ARGS";
        case (NC_EMULTIDEFINE_FILL_MODE):	return "NC_EMULTIDEFINE_FILL_MODE";
        case (NC_EMULTIDEFINE_VAR_FILL_MODE):	return "NC_EMULTIDEFINE_VAR_FILL_MODE";
        case (NC_EMULTIDEFINE_VAR_FILL_VALUE):	return "NC_EMULTIDEFINE_VAR_FILL_VALUE";
        case (NC_EMULTIDEFINE_CMODE):		return "NC_EMULTIDEFINE_CMODE";
        default:
              sprintf(unknown_str,"Unknown code %d",err);
    }
    return unknown_str;
}

