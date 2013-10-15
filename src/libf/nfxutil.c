
#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include "nc.h"
#include "mpinetcdf_impl.h"
int ncmpixVardim( int ncid, int varid ) 
{
  NC_var *varp;
  NC     *ncp;
  int    status;
#ifdef ENABLE_SUBFILING
  int    ndims_org;
#endif

  status = ncmpii_NC_check_id(ncid, &ncp);
  if (status != NC_NOERR) return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if (varp == NULL) return NC_ENOTVAR;
#ifndef ENABLE_SUBFILING
  return varp->ndims;
#else
  /* NOTE: it should check varp->num_subfiles 
     because some variables can't be subfiled */ 
  if (varp->num_subfiles > 1) { 
      status = ncmpi_get_att_int (ncid, varid, 
				  "ndims_org", &ndims_org); 
      return ndims_org; 
  }
  else
      return varp->ndims;
#endif
}
