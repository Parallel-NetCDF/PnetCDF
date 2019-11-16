# Note for Developers

### Table of contents
- [Future Work]

---
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
