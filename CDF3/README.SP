# John Tannahill <tannahill1@llnl.gov> reported success building
# parallel-netcdf on the 'seaborg' cluster (an IBM-SP at NERSC) and Tyce
# Mclarty <mclarty3@llnl.gov> reported success on LLNL's 'frost' cluster by
# setting these environment variables:

   setenv MPICC  mpcc_r
   setenv MPIF77 mpxlf_r
   setenv F77    xlf
   setenv FC     xlf
   setenv CC     xlc
   setenv CXX    xlC

# after setting these variables, configure/make/make install should "just work"


# We also successfully tested 64-bit mode on 'DataStar' (an IBM-SP at SDSC)
# with environment variables:

   setenv OBJECT_MODE 64
   setenv MPICC  mpcc_r
   setenv MPIF77 mpxlf_r
   setenv F77    xlf
   setenv FC     xlf
   setenv CC     xlc
   setenv CXX    xlC
   setenv CFLAGS -q64
   setenv FFLAGS -q64
   setenv F90FLAGS -q64
   setenv CXXFLAGS -q64
