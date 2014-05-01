#include "ncmpiString.h"
#include <pnetcdf.h>
using namespace PnetCDF;

// create an instance of NcmpiString  called PnetCDF::ncmpiString
namespace PnetCDF {
  NcmpiString ncmpiString;
}

// constructor
NcmpiString::NcmpiString() : NcmpiType(NC_STRING){
}

NcmpiString::~NcmpiString() {
}


// equivalence operator
bool NcmpiString::operator==(const NcmpiString & rhs)    {
  // simply check the netCDF id.
  return myId == rhs.myId;
}  
