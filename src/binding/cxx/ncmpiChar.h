#include "ncmpiType.h"

#ifndef NcmpiCharClass
#define NcmpiCharClass

namespace PnetCDF
{

  /*! Class represents a netCDF atomic Char type. */
  class PNETCDF_PUBLIC_API NcmpiChar : public NcmpiType
  {
  public:

    /*! equivalence operator */
    bool operator==(const NcmpiChar & rhs);

    ~NcmpiChar();

    /*! Constructor */
    NcmpiChar();
  };

  /*! A global instance  of the NcmpiChar class within the netCDF namespace. */
  extern PNETCDF_PUBLIC_API NcmpiChar ncmpiChar;

}
#endif
