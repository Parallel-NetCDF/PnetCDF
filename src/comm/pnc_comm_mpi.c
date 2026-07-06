/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pnetcdf_comm.h>

const char*
ncmpix_comm_backend(void)
{
    return "mpi";
}

int
ncmpix_create(PNC_Comm    comm,
              const char *path,
              int         cmode,
              PNC_Info    info,
              int        *ncidp)
{
    return ncmpi_create(comm, path, cmode, info, ncidp);
}

int
ncmpix_open(PNC_Comm    comm,
            const char *path,
            int         omode,
            PNC_Info    info,
            int        *ncidp)
{
    return ncmpi_open(comm, path, omode, info, ncidp);
}
