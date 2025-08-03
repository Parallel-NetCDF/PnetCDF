export PNETCDF_HINTS="nc_hash_size_dim=4096;nc_hash_size_var=4096"
mpiexec -n 4 ./new_format_test_all "data/dataset_568k_metadata" "data/dataset_568k_new_format_test_all.pnc"
export PNETCDF_HINTS="nc_hash_size_dim=16384;nc_hash_size_var=16384"
mpiexec -n 4 ./new_format_test_all "data/dataset_5m_metadata" "data/dataset_5m_new_format_test_all.pnc"
# optinoal: read the created file
mpiexec -n 4 ./new_format_read_test_all "data/dataset_568k_new_format_test_all.pnc"
