#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pnetcdf_comm.h>

const char*
ncmpix_comm_backend(void)
{
    return "commux";
}
