export PNETCDF_HINTS="nc_hash_size_dim=1048576;nc_hash_size_var=1048576"
mpiexec -n 4 ./lib_baseline_test_all data/dataset_568k_metadata ./data/dataset_568k_metadata_lib_baseline_test_all.nc

export PNETCDF_HINTS="nc_hash_size_dim=16777216;nc_hash_size_var=16777216"
mpiexec -n 4 ./lib_baseline_test_all data/dataset_5m_metadata ./data/dataset_5m_metadata_lib_baseline_test_all.nc