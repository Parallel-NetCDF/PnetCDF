/*********************************************************************
 *
 *  Copyright (C) 2012, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This program is adopted from the netCDFexample program, URL:
   http://www.unidata.ucar.edu/software/netcdf/examples/programs/SimpleXyWr.cpp

   This example writes a 2D array of sample data. We create two shared
   dimensions, "x" and "y", and a netCDF variable, called "data".
*/

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <pnetcdf>

using namespace std;
using namespace PnetCDF;
using namespace PnetCDF::exceptions;

// We are writing 2D data, a 6 x 12 grid. 
static const int NX = 6;
static const int NY = 12;

// Return this in event of a problem.
static const int NC_ERR = 2;

int main(int argc, char *argv[])
{
    int err=0;

    MPI_Init(&argc, &argv);

    // This is the data array we will write. It will just be filled
    // with a progression of numbers for this example.
    int dataOut[NX][NY];
  
    // Create some pretend data. If this wasn't an example program, we
    // would have some real data to write, for example, model output.
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            dataOut[i][j] = i * NY + j;
  
    // The default behavior of the C++ API is to throw an exception i
    // an error occurs. A try catch block is necessary.
   
    try {  
        // Create the file. The Replace parameter tells netCDF to overwrite
        // this file, if it already exists.
        NcmpiFile dataFile(MPI_COMM_WORLD, "simple_xy.nc", NcmpiFile::replace);
      
        // Create netCDF dimensions
        NcmpiDim xDim = dataFile.addDim("x", NX);
        NcmpiDim yDim = dataFile.addDim("y", NY);
      
        // Define the variable. The type of the variable in this case is
        // ncInt (32-bit integer).
        vector<NcmpiDim> dims;
        dims.push_back(xDim);
        dims.push_back(yDim);
        NcmpiVar data = dataFile.addVar("data", ncmpiInt, dims);
   
        // Write the data to the file. Although netCDF supports
        // reading and writing subsets of data, in this case we write all
        // the data in one operation.
        data.putVar_all(&dataOut[0][0]);
      
        // The file will be automatically close when the NcmpiFile object goes
        // out of scope. This frees up any internal netCDF resources
        // associated with the file, and flushes any buffers.
      
        cout << "*** SUCCESS writing example file simple_xy.nc!" << endl;
    }
    catch(NcmpiException& e) {
        cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
        err = 1;
    }
    MPI_Finalize();
    return err;
}
