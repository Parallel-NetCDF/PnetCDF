#include <stdio.h>
#include <string.h>
#include <iostream>
#include <pnetcdf>

using namespace std;
using namespace PnetCDF;
using namespace PnetCDF::exceptions;

int main( int argc, char *argv[] )
{
   char filename[128];
   int rank, pass=1, verbose=0;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (argc > 2) {
       if (!rank) printf("Usage: %s [filename]\n",argv[0]);
       MPI_Finalize();
       return 0;
   }
   strcpy(filename, "testfile.nc");
   if (argc == 2) strcpy(filename, argv[1]);

   try
   {
      if (verbose) cout << "Test creation of classic format file" << endl;
      {
	 NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::replace,
                          NcmpiFile::classic);
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
         if (verbose) cout << "testing the switch to DATA mode..." << endl;
         int arr[] = {1,2,3,4,5,6,7,8,9,10,11};
         var_gw.putVar_all(arr);
      }

      // Now test reading.
      {
	 NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::read);

	 if (ncFile.getVarCount() != 2)
	    throw NcmpiException( "Holy Mother of Pearl!", __FILE__, __LINE__);
      }

      // and redefinition
      {
        NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::write);
        if (verbose) cout << "testing the switch to DEFINE mode..." << endl;
        ncFile.putAtt(string("name"),string("value"));
      }

      if (verbose) cout << "    -----------   passed\n";
   }
   catch(NcmpiException& e)
   {
      cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
      pass = 0;
   }

    char cmd_str[80];
    sprintf(cmd_str, "*** TESTING C++ %s for creation of classic format file", argv[0]);
    if (rank == 0) {
        if (pass) printf("%-66s ------ pass\n", cmd_str);
        else      printf("%-66s ------ failed\n", cmd_str);
    }

   MPI_Finalize();
   return 0;
}
