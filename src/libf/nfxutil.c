
#include "nc.h"
#include "mpinetcdf_impl.h"
int ncmpixVardim( int ncid, int varid ) 
{
  NC_var *varp;
  NC     *ncp;

  ncmpii_NC_check_id(ncid, &ncp);
  varp = ncmpii_NC_lookupvar(ncp, varid);
  return varp->ndims;
}
