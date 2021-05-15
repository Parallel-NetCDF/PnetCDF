## PnetCDF subfiling feature test

In order to use the subfiling module, pnetcdf must be configured with
"--enable-subfiling" option.

This is a test run that will create two subfiles, in addition to the master
file, the file used in the call to ncmpi_creat(). Note that the number of MPI
processes must be equal or larger than the number of subfiles.

Run command:
```
$ mpiexec -n 2 ./test_subfile -f test_subfile.nc -s 2
```

This will create 1 master file (the original file, i.e., test_subfile.nc) and
two subfiles:
```
    test_subfile.nc.subfile_0.nc
    test_subfile.nc.subfile_1.nc
```
The master file contains only metadata about the information of subfiling. The
variable data is actually stored in the two subfiles.

In this example program, a 3D integer array, named var0_0, of size 8 x 16 x 8
is created, which is partitioned along the most significant dimension into two
subarrays of size 8 x 16 x 8 each. The two subarrays are saved in two subfiles
separately. The first half is stored in file test_subfile.nc.subfile_0.nc and
the second half is stored in file test_subfile.nc.subfile_1.nc.

Below shows the header of the master file, test_subfile.nc.
```
% ncmpidump -h test_subfile.nc
netcdf test_subfile {
// file format: CDF-5 (big variables)
dimensions:
	dim0_0 = 16 ;
	dim0_1 = 16 ;
	dim0_2 = 8 ;
variables:
	int var0_0 ;
		var0_0:_PnetCDF_SubFiling.par_dim_index = 0 ;
		var0_0:_PnetCDF_SubFiling.ndims_org = 3 ;
		var0_0:_PnetCDF_SubFiling.dimids_org = 0, 1, 2 ;
		var0_0:_PnetCDF_SubFiling.num_subfiles = 2 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_0).subfile.0 = 0, 7 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_1).subfile.0 = 0, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_2).subfile.0 = 0, 7 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_0).subfile.1 = 8, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_1).subfile.1 = 0, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_2).subfile.1 = 0, 7 ;

// global attributes:
		:_PnetCDF_SubFiling.num_subfiles = 2 ;
}
```

Below shows the header of subfile test_subfile.nc.subfile_0.nc:
```
% ncmpidump -h test_subfile.nc.subfile_0.nc
netcdf test_subfile.nc.subfile_0 {
// file format: CDF-5 (big variables)
dimensions:
	dim0_0.var0_0 = 8 ;
	dim0_1.var0_0 = 16 ;
	dim0_2.var0_0 = 8 ;
variables:
	int var0_0(dim0_0.var0_0, dim0_1.var0_0, dim0_2.var0_0) ;
		var0_0:_PnetCDF_SubFiling.range(dim0_0).subfile.0 = 0, 7 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_1).subfile.0 = 0, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_2).subfile.0 = 0, 7 ;
		var0_0:_PnetCDF_SubFiling.subfile_index = 0 ;
}
```

Below shows the header of subfile test_subfile.nc.subfile_1.nc:
```
% ncmpidump -h test_subfile.nc.subfile_1.nc
netcdf test_subfile.nc.subfile_1 {
// file format: CDF-5 (big variables)
dimensions:
	dim0_0.var0_0 = 8 ;
	dim0_1.var0_0 = 16 ;
	dim0_2.var0_0 = 8 ;
variables:
	int var0_0(dim0_0.var0_0, dim0_1.var0_0, dim0_2.var0_0) ;
		var0_0:_PnetCDF_SubFiling.range(dim0_0).subfile.1 = 8, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_1).subfile.1 = 0, 15 ;
		var0_0:_PnetCDF_SubFiling.range(dim0_2).subfile.1 = 0, 7 ;
		var0_0:_PnetCDF_SubFiling.subfile_index = 1 ;
}
```

* See COPYRIGHT notice in top-level directory.
