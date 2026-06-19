#include "ncmpiType.h"

#ifndef NcmpiUint64Class
#define NcmpiUint64Class

namespace PnetCDF
{

  /*! Class represents a netCDF atomic Uint64 type.*/
  class PUBLIC_API NcmpiUint64 : public NcmpiType
  {
  public:

    /*! equivalence operator */
    bool operator==(const NcmpiUint64 & rhs);

    /*! destructor */
    ~NcmpiUint64();

    /*! Constructor */
    NcmpiUint64();
  };

  /*! A global instance  of the NcmpiUint64 class within the netCDF namespace. */
  extern PUBLIC_API NcmpiUint64 ncmpiUint64;

}
#endif
