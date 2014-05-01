#include <iostream>
#include <pnetcdf>

using namespace std;
using namespace PnetCDF;
using namespace PnetCDF::exceptions;

int main( int argc, char *argv[] )
{
   MPI_Init(&argc, &argv);

   try
   {
      cout << "Test creation of classic format file" << endl;
      {
	 NcmpiFile ncFile(MPI_COMM_WORLD, "test_classic.nc", NcmpiFile::replace, NcmpiFile::classic);
	 NcmpiDim dim1 = ncFile.addDim("dim1",11);
	 NcmpiDim dim2 = ncFile.addDim("dim2");
	 NcmpiDim dim3 = ncFile.addDim("dim3",13);

	 NcmpiVar var_gw  = ncFile.addVar("George_Washington", ncmpiInt, dim1);
	 // add a 2D record variable
	 vector<NcmpiDim> dimArray(2);
	 dimArray[0]=dim2;
	 dimArray[1]=dim1;
	 NcmpiVar varA1_3  = ncFile.addVar("varA1_3", ncmpiInt, dimArray);

	 // ncFile.enddef(); is no need in C++ program

         // and inserting some data that needs leaving the define mode
         cout << "testing the swith to DATA mode..." << endl;
         int arr[] = {1,2,3,4,5,6,7,8,9,10,11};
         var_gw.putVar_all(arr);
      }

      // Now test reading.
      {
	 NcmpiFile ncFile(MPI_COMM_WORLD, "test_classic.nc", NcmpiFile::read);

	 if (ncFile.getVarCount() != 2)
	    throw NcmpiException( "Holy Mother of Pearl!", __FILE__, __LINE__);
      }

      // and redefinition
      {
        NcmpiFile ncFile(MPI_COMM_WORLD, "test_classic.nc", NcmpiFile::write);
        cout << "testing the swith to DEFINE mode..." << endl;
        ncFile.putAtt(string("name"),string("value"));
      }

      cout << "    -----------   passed\n";
   }
   catch(NcmpiException& e)
   {
      cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
      return 1;
   }

   MPI_Finalize();
   return 0;
}
