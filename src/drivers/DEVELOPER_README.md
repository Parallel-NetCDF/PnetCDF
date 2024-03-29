## PnetCDF source code naming conventions

PnetCDF uses the following naming conventions for folders, files, and functions
under this directory, `src/drivers`.

* Folder `common/`
  * Functions implemented in this folder are shared by all drivers, and
    dispatchers.
  * Files should use names easily recognized for the subroutines contained.
  * All non-static functions should use prefix name "ncmpii_".

* Folder `include/`
  * Header files stored in this folder are shared by all drivers and
    dispatchers.
  * Files should use prefix name "ncmpii_".
  * Definition of variables, constants, and struct should also use prefix name
    "ncmpii_"

* Folder `ncmpio/`
  * Files in this folder are the implementation of PnetCDF using MPI-IO.
  * Files should use prefix name "ncmpio_".
  * Non-static functions should also use prefix name "ncmpio_".

* If a new driver, say "foo", is added, please do the followings.
  1. Create a new folder with name "ncfoo"
  2. Store all files implementing this driver under `ncfoo/`
  3. All files stored under `ncfoo/` should use prefix name "ncfoo_".
  4. Non-static functions should use prefix name "ncfoo_".

