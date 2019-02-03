#ifndef _nczipio_INTERNAL_H
#define _nczipio_INTERNAL_H

#include "nczipio_driver.h"

extern int 
nczipioi_init(NC_zip*);

extern int 
nczipioi_var_list_init(NC_zip_var_list*);

extern int 
nczipioi_var_list_free(NC_zip_var_list*);

extern int 
nczipioi_var_list_add(NC_zip_var_list*, NC_zip_var); 

extern int 
nczipioi_var_init(NC_zip*, NC_zip_var*);

extern int 
nczipioi_extract_hint(NC_zip*, MPI_Info);

extern int 
nczipioi_export_hint(NC_zip *nczipp, MPI_Info info);

extern int
nczipioi_get_var(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*, MPI_Offset, MPI_Datatype, int);

extern int
nczipioi_get_varn(NC_zip*, NC_zip_var*, int, MPI_Offset* const*, MPI_Offset* const*, void*, MPI_Offset, MPI_Datatype, int);

extern int
nczipioi_put_var(NC_zip*, NC_zip_var*, const MPI_Offset*, const MPI_Offset*, const MPI_Offset*, void*);

extern MPI_Offset 
NC_Type_size(nc_type);

#endif