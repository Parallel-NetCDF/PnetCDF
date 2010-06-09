mpiexec -n 4 mcoll_perf
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.0.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.1.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.2.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.3.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.4.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.5.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test0.nc nc_2.4.6.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test1.nc nc_2.4.7.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test1.nc nc_2.4.8.nc
mpiexec -n 4 ../../src/utils/ncdiff/ncmpidiff ../data/mcoll_perf_test1.nc nc_2.4.9.nc

