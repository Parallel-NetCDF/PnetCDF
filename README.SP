# John Tannahill <tannahill1@llnl.gov> reported success building
# parallel-netcdf on the 'seaborg' cluster (an IBM-SP at NERSC) by setting
# these environment variables:

   setenv MPICC  mpcc_r
   setenv MPIF77 mpxlf_r
   setenv F77    xlf
   setenv FC     xlf
   setenv CC     xlc
   setenv CXX    xlC

# after setting these variables, configure/make/make install should "just work"
