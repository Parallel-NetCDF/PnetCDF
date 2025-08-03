# Data Object Creation Performance with Parallel I/O library

This artifact repository accompanies our study on scalable parallel metadata management in paralllel I/O libraries. It evaluates three approaches for large-volume data object creation: (1) an application-level approach using the default I/O library that gathers metadata on all processes in order to generate a consistent copy before before invoking library APIs, (2) a library-level approach with internal metadata exchange and consistency check which enables independent object creation, and (3) a new file header format with partitioned metadata blocks that enables multiple MPI processes to create data objects independently and write
metadata to the file header in parallel. The latter two approaches require modifications to the I/O library. This program reproduces the results from the paper and compares scalability across all approaches.

### Data Object Creation Benchmark and Performance Measure

The test experiments in this repository are designed to measure the performance and scalability of parallel data object creation using two datasets derived from the Exa.TrkX project, a High Energy Physics application focused on neutrino trajectory reconstruction. The first dataset, referred to as dataset 560k, contains 568,480 data arrays organized into 35,530 groups when stored in HDF5. The second dataset, dataset 5m, contains 5,684,800 data arrays (10× larger) and is used to evaluate performance under increased metadata volume. Both datasets are evenly partitioned among MPI processes to simulate the graph generation workload. In the real-world graph generation application, metadata for each data object is initially present in each process’s memory as the result of parallel computation. To mimic this condition in our tests, the test program first reads the metadata into memory from a prepared binary input file, and then performs the write operation on a shared output file using parallel I/O library, starting from data object creation.

We measure runtime performance by collecting the end-to-end time, starting after metadata is loaded into memory and ending once all data objects are created and metadata is written to file. To identify major cost components, the test program profiles the runtime into metadata exchange time (MPI communication), metadata consistency check time (string and numerical comparisons), and other costs such as file writing and closing. All timings are reported as the maximum across all processes, except for the “others” category, which is calculated as the difference between the total end-to-end time and the sum of the metadata exchange and consistency check times. The detailed timing breakdowns collected are shown below:
* Application-level approach:  Metadata exchange time is collected at the test program level as the time to serialize metadata, perform MPI_Allgather, and deserialize the result. Consistency check time is collected as the time to create data objects using library APIs (e.g., `ncmpi_def_var`, `ncmpi_def_dim` with PnetCDF), measured with a timer around data object creation stage.
* Library-level approach: Metadata exchange time is collected within the library, using a timer wrapped around the relevant code inside the `ncmpi_enddef` function. Consistency check time is also collected within the library, measured as the time taken by `ncmpi_enddef` to reconstruct a complete metadata copy on each process
* New header format approach: Metadata exchange time is collected as the time to synchronize metadata block information (index table), measured using a timer wrapped around the relevant code inside `ncmpi_enddef`. Consistency check time includes two components: (1) collected at the test program level as the time to create independent data objects on each process, using a timer wrapped around the object creation calls where intra-block consistency is checked; and (2) collected within `ncmpi_enddef` as the time to construct the metadata block index, where inter-block consistency and intra-block consistency for shared metadata block are checked. Given the dataset used and the scale of experiments in our study, the second component is small in pratice and can be considered negligible.

### Software Libraries

We use the following software libraries. Here, we provide scripts to install them on Perlmutter. The installation process on other platforms should be similar.

1. HDF5 1.14.4-2 (for application-level baseline approach using HDF5)
    ```shell
    # download source codes
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.14.4-2/src/hdf5-1.14.4-2.tar.gz
    tar -xf hdf5-1.14.4-2.tar.gz
    cd hdf5-1.14.4-2

    # prefix of install dir. One should modify its value. (optional)
    export PREFIX=/path/to/hddf/install

    # configure
    ../configure --prefix=${PREFIX} --enable-parallel --enable-build-mode=production
    
    # compile
    make

    # install
    make install
    ```