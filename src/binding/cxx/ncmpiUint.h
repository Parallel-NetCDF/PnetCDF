#include "ncmpiType.h"

#ifndef NcmpiUintClass
#define NcmpiUintClass

namespace PnetCDF
{

  /*! Class represents a netCDF atomic Uint type. */
  class PUBLIC_API NcmpiUint : public NcmpiType
  {
  public:

    /*! equivalence operator */
    bool operator==(const NcmpiUint & rhs);

    /*! destructor */
    ~NcmpiUint();

    /*! Constructor */
    NcmpiUint();
  };

  /*! A global instance  of the NcmpiUint class within the netCDF namespace. */
  extern PUBLIC_API NcmpiUint ncmpiUint;

}
#endif
