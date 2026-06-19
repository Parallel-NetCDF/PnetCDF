#include "ncmpiType.h"

#ifndef NcmpiDoubleClass
#define NcmpiDoubleClass

namespace PnetCDF
{

  /*! Class represents a netCDF atomic Double type. */
  class PUBLIC_API NcmpiDouble : public NcmpiType
  {
  public:

    /*! equivalence operator */
    bool operator==(const NcmpiDouble & rhs);

    /*!  destructor */
    ~NcmpiDouble();

    /*! Constructor */
    NcmpiDouble();
  };

  /*! A global instance  of the NcmpiDouble class within the netCDF namespace. */
  extern PUBLIC_API NcmpiDouble ncmpiDouble;

}
#endif
