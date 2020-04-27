# Note for Developers

### Table of contents
- [Future Work]
- [Internal global attributes]
- [Anchor variable (one per variable with chunking enabled)]
- [Reference table]
- [Chunks]
- [Requirement for compressed variables]

---

## Internal global attributes:
  * Number of chunked variables

## Anchor variable (one per variable with chunking enabled):
  * A scalar variable
  * Data type is the same as user defined
  * Internal attributes
    + Dimension IDs are saved as an attribute of an array of integer type
    + Number of dimensions is saved as an internal attribute
    + An attribute to tell whether it is a fixed-size or record variable
    + An attribute offset pointer to reference table
      * For fixed-size variable, it is a scalar
      * For record variable, it is an array of 8-type integers, one for each record
        * This array can be allocated in multiple of 16 for example
        * Need an integer for allocated size, e.g. multiple of 16
        * Need an integer for size (true number of records written)
    + An attributes for chunk sizes, an integer array
    + An attributes for compression algorithm
    + An attributes for compression level
  * If a variable missing these internal attributes, it is a traditional variable

## Reference table:
  * An array stores offsets of individual chunks
  * Not a NetCDF variable. But we use the CDF5 format specification to define it
    + TODO: give it a formal spec in BNF grammar
  * For a fixed-size variable, it is a 1D array of size equal to the number of chunks
  * This table is loaded into memory when calling ncmpi_inq_varid
  * For blocking API, it is sync-ed and written to file by root
    + TODO: in future, it can be written by multiple ranks in parallel
  * For nonblocking API, multiple tables are written by multiple ranks in parallel

## Chunks:
  * Chunks are not NetCDF variables
    + TODO: give it a formal spec in BNF grammar?
  * Chunks are stored in space between NetCDF variables, i.e. padding areas in files
  * Data is type-converted and byte-swapped before compression
  * In principle, chunks should be stored in file contiguously with each other,
    for all variables. But they are not required to be stored contiguously.
  * The storage order of chunks is in row major

## Requirement for compressed variables:
  * Collective I/O only (this is the same required by HDF5)
  * Must be chunked (same as HDF5)


## Future Work
* Reuse metadata accross variables
  - Variable from same simulation space may have same access apttern.
  - Instead of generating variable metadata and indexx table separately, we can 
    share information accross variables.
  - Chunk sizeand chunk ownership info can be reused.
* Data seiving
  - When rewriting to a chunk, we do't need to read the background if it is 
    fully overwritten.
  - Need an efficient way to determine whether a chunk is fully rewrititen. 
  - It may be infesible due to communication and computation cost.
  - HDF5 approximate this by checking if owner fully rewriten the chunk.
* Reuse metadata accross records
  - I/O pattern accross time steps are likely the same.
  - If we detect same I/O pattern as previous record, we can skip sending the metadata. 
  - MPI datatype created for previous timestep can also be reused. 
---
