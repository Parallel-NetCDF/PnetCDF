/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs
 *
 * ncmpi_create()           : dispatcher->create()
 * ncmpi_open()             : dispatcher->open()
 * ncmpi_close()            : dispatcher->close()
 * ncmpi_enddef()           : dispatcher->enddef()
 * ncmpi__enddef()          : dispatcher->_enddef()
 * ncmpi_redef()            : dispatcher->redef()
 * ncmpi_begin_indep_data() : dispatcher->begin_indep_data()
 * ncmpi_end_indep_data()   : dispatcher->end_indep_data()
 * ncmpi_abort()            : dispatcher->abort()
 * ncmpi_inq()              : dispatcher->inq()
 * ncmpi_inq_misc()         : dispatcher->inq_misc()
 * ncmpi_wait()             : dispatcher->wait()
 * ncmpi_wait_all()         : dispatcher->wait()
 * ncmpi_cancel()           : dispatcher->cancel()
 *
 * ncmpi_set_fill()         : dispatcher->set_fill()
 * ncmpi_fill_var_rec()     : dispatcher->fill_rec()
 * ncmpi_def_var_fill()     : dispatcher->def_var_fill()
 * ncmpi_inq_var_fill()     : dispatcher->inq()
 *
 * ncmpi_sync()             : dispatcher->sync()
 * ncmpi_sync_numrecs()     : dispatcher->sync_numrecs()
 *
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>
#include <netcdf_par.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strlen() */

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

size_t nc4io_nc_type_size(nc_type type){
    if (type == NC_BYTE){
        return 1;
    }
    else if (type == NC_UBYTE){
        return 1;
    }
    else if (type == NC_CHAR){
        return 1;
    }
    else if (type == NC_STRING){
        return 1;
    }
    else if (type == NC_SHORT){
        return 2;
    }
    else if (type == NC_USHORT){
        return 2;
    }
    else if (type == NC_INT){
        return 4;
    }
    else if (type == NC_LONG){
        return 4;
    }
    else if (type == NC_UINT){
        return 4;
    }
    else if (type == NC_FLOAT){
        return 4;
    }
    else if (type == NC_INT64){
        return 8;
    }
    else if (type == NC_UINT64){
        return 8;
    }
    else if (type == NC_DOUBLE){
        return 8;
    }
    return 0;
}
