!
!  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

    ! This function gets the executable name and output file name from the
    ! command line.
    integer function get_args(cmd, filename)
#ifdef NAGf90Fortran
        USE F90_UNIX_ENV, only : iargc, getarg
        implicit none
#else
        implicit none
        integer iargc
#endif
        integer argc
        character(len=*) cmd, filename

        get_args = 1
        call getarg(0, cmd)
        argc = IARGC()
        if (argc .GT. 1) then
           print*,'Usage: ',trim(cmd),' [filename]'
           get_args = 0
           return
        endif
        if (argc .EQ. 1) call getarg(1, filename)
    end function get_args

