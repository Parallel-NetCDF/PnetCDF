
#include "nc.h"
#include "mpinetcdf_impl.h"
int ncmpixVardim( int ncid, int varid ) 
{
  NC_var *varp;
  NC     *ncp;
  int    status;

  status = ncmpii_NC_check_id(ncid, &ncp);
  if (status != NC_NOERR) return status;
  varp = ncmpii_NC_lookupvar(ncp, varid);
  if (varp == NULL) return NC_ENOTVAR;
  return varp->ndims;
}
