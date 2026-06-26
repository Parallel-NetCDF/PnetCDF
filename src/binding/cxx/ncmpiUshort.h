#include "ncmpiType.h"

#ifndef NcmpiUshortClass
#define NcmpiUshortClass

namespace PnetCDF
{

  /*! Class represents a netCDF atomic Ushort type. */
  class PNETCDF_PUBLIC_API NcmpiUshort : public NcmpiType
  {
  public:

    /*! equivalence operator */
    bool operator==(const NcmpiUshort & rhs);

    /*! destructor */
    ~NcmpiUshort();

    /*! Constructor */
    NcmpiUshort();
  };

  // declare that the class instance ncmpiUshort is known by all....
  extern PNETCDF_PUBLIC_API NcmpiUshort ncmpiUshort;

}
#endif
