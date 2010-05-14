/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: error.c,v 1.14 90/02/23 16:08:55 davis Exp */

/*LINTLIBRARY*/

#include "ncconfig.h"

#include <stdio.h>
#include <stdarg.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include "pnetcdf.h"
#include "nc.h"

#ifndef NO_STRERROR
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
#endif /* NO_STRERROR */


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


static char unknown[] = "Unknown Error";


const char *
ncmpi_strerror(int err)
{

#ifdef vms 
	if(err == EVMSERR)
	{
		return vms_strerror(err);
	}	
	/* else */
#endif /* vms */

	if(NC_ISSYSERR(err))
	{
		const char *cp = (const char *) strerror(err);
		if(cp == NULL)
			return unknown;
		/* else */
		return cp;
	}
	/* else */

	switch (err) {
	case NC_NOERR:
	    return "No error";
	case NC_ESMALL:
	    return "Size of off_t too small for requested format ";
	case NC_ENOTINDEP:
	    return "Operation not allowed in collective data mode";
	case NC_EINDEP:
	    return "Operation not allowed in independent data mode";
	case NC_EFILE:
	    return "Unknown error in file operation";
        case NC_EREAD:
            return "Unknow error occurs in reading file";
	case NC_EWRITE:
	    return "Unknow error occurs in writting file";
	case NC_EMULTIDEFINE:
	    return "File header is inconsistent among processes";
	case NC_EOFILE:
	    return "Can not open/create file";
	case NC_EBADID:
	    return "Not a netCDF id";
	case NC_ENFILE:
	    return "Too many netCDF files open";
	case NC_EEXIST:
	    return "netCDF file exists && NC_NOCLOBBER";
	case NC_EINVAL:
	    return "Invalid argument";
	case NC_EPERM:
	    return "Write to read only";
	case NC_ENOTINDEFINE:
	    return "Operation not allowed in data mode";
	case NC_EINDEFINE:
	    return "Operation not allowed in define mode";
	case NC_EINVALCOORDS:
	    return "Index exceeds dimension bound";
	case NC_EMAXDIMS:
	    return "NC_MAX_DIMS exceeded";
	case NC_ENAMEINUSE:
	    return "String match to name in use";
	case NC_ENOTATT:
	    return "Attribute not found";
	case NC_EMAXATTS:
	    return "NC_MAX_ATTRS exceeded";
	case NC_EBADTYPE:
	    return "Not a netCDF data type or _FillValue type mismatch";
	case NC_EBADDIM:
	    return "Invalid dimension id or name";
	case NC_EUNLIMPOS:
	    return "NC_UNLIMITED in the wrong index";
	case NC_EMAXVARS:
	    return "NC_MAX_VARS exceeded";
	case NC_ENOTVAR:
	    return "Variable not found";
	case NC_EGLOBAL:
	    return "Action prohibited on NC_GLOBAL varid";
	case NC_ENOTNC:
	    return "Not a netCDF file";
	case NC_ESTS:
	    return "In Fortran, string too short";
	case NC_EMAXNAME:
	    return "NC_MAX_NAME exceeded";
	case NC_EUNLIMIT:
	    return "NC_UNLIMITED size already in use";
	case NC_ENORECVARS:
	    return "nc_rec op when there are no record vars";
	case NC_ECHAR:
	    return "Attempt to convert between text & numbers";
	case NC_EEDGE:
	    return "Edge+start exceeds dimension bound";
	case NC_ESTRIDE:
	    return "Illegal stride";
	case NC_EBADNAME:
	    return "Attribute or variable name contains illegal characters";
	case NC_ERANGE:
	    return "Numeric conversion not representable";
	case NC_ENOMEM:
	    return "Memory allocation (malloc) failure";
	case NC_EVARSIZE:
	    return "One or more variable sizes violate format constraints";
	case NC_EDIMSIZE:
	    return "Invalid dimension size";
	case NC_EMULTITYPES:
	    return "Multiple types used in memory data";
	case NC_EIOMISMATCH:
	    return "Input/Output data amount mismatch";
	case NC_ENEGATIVECNT:
	    return "Negative count is prohibited";
	case NC_EUNSPTETYPE:
	    return "Unsupported etype is used in MPI datatype for memory data";
	case NC_EDIMS_NELEMS_MULTIDEFINE:
	    return "Number of dimensions is defined inconsistently among processes.";
	case NC_EDIMS_SIZE_MULTIDEFINE:
	    return "Dimension size is defined inconsistently among processes.";
	case NC_EVARS_NELEMS_MULTIDEFINE:
	    return "Number of variables is defined inconsistently among processes.";
	case NC_EVARS_NDIMS_MULTIDEFINE:
	    return "Dimensionality of this variable is defined inconsistently among processes.";
	case NC_EVARS_DIMIDS_MULTIDEFINE:
	    return "Dimension IDs used to define this variable is inconsistent among processes.";
	case NC_EVARS_TYPE_MULTIDEFINE:
	    return "Data type of this variable is defined inconsistently among processes.";
	case NC_EVARS_LEN_MULTIDEFINE:
	    return "Total number of elements of this variable is defined inconsistently among processes.";
	case NC_EVARS_BEGIN_MULTIDEFINE:
	    return "(Internal error) beginning file offset of this variable is inconsistent among processes.";
	case NC_ENUMRECS_MULTIDEFINE:
	    return "Number of records is inconsistent among processes.";
	case NC_EINVAL_REQUEST:
	    return "Invalid nonblocking request ID.";
	
	default:
	    return unknown;
	}
	/* default */
	return unknown;
}

void ncmpii_handle_error(int rank, int mpi_status, char * msg)
{
	char errorString[MPI_MAX_ERROR_STRING];
	int errorStringLen;
	MPI_Error_string(mpi_status, errorString, &errorStringLen);
	printf("%2d: %s : %s\n", rank, msg, errorString);
}
