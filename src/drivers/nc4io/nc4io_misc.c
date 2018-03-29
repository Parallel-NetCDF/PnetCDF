#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>
#include <netcdf_par.h>

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int nc4io_nc_type_size(nc_type type){
    
    if (type == NC_NAT){
        return 0;
    }
    else if (type == NC_BYTE){
        return 1;
    }
    else if (type == NC_CHAR){
        return 1;
    }
    else if (type == NC_SHORT){
        return 2;
    }
    else if (type == NC_INT){
        return 4;
    }
    else if (type == NC_LONG){
        return 4;
    }
    else if (type == NC_FLOAT){
        return 4;
    }
    else if (type == NC_DOUBLE){
        return 8;
    }
    else if (type == NC_UBYTE){
        return 1;
    }
    else if (type == NC_USHORT){
        return 2;
    }
    else if (type == NC_UINT){
        return 4;
    }
    else if (type == NC_INT64){
        return 8;
    }
    else if (type == NC_UINT64){
        return 8;
    }
    else if (type == NC_STRING){
        return 1;
    }
    else{
        return 0;
    }

    return 0;
}