#include "ncmpiType.h"

#ifndef NcmpiStringClass
#define NcmpiStringClass

namespace PnetCDF
{
  
  /*! Class represents a netCDF atomic String type. */
  class NcmpiString : public NcmpiType
  {
  public: 
    
    /*! equivalence operator */
    bool operator==(const NcmpiString & rhs);
    
    /*! destructor */
    ~NcmpiString();
    
    /*! Constructor */
    NcmpiString();
  };

  /*! A global instance  of the NcmpiString class within the netCDF namespace. */
  extern NcmpiString ncmpiString;

}
#endif
