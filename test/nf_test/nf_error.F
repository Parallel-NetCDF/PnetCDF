#if 0
    Copyright 1996, UCAR/Unidata
    See netcdf/COPYRIGHT file for copying and redistribution conditions.
    $Id$
#endif


C
C Use for logging error messages
C
        subroutine error(msg)
        use pnetcdf
        implicit        none
        character*(*)   msg
#include "tests.inc"

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg
        end


C
C Use for logging error conditions
C
        subroutine errori(msg, i)
        use pnetcdf
        implicit        none
        character*(*)   msg
        integer         i
#include "tests.inc"

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, i
        end


C
C Use for logging error conditions
C
        subroutine errord(msg, d)
        use pnetcdf
        implicit        none
        character*(*)   msg
        doubleprecision d
#include "tests.inc"

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, d
        end


C
C Use for logging error conditions
C
        subroutine errorc(msg, string)
        use pnetcdf
        implicit        none
        character*(*)   msg
        character*(*)   string
#include "tests.inc"

        nfails = nfails + 1
        if (nfails .le. max_nmpt) print *, msg, 
     +      string(1:len_trim(string))
        end


C
C Use for logging error conditions
C
        subroutine errore(msg, err)
        use pnetcdf
        implicit        none
        character*(*)   msg
        integer         err
#include "tests.inc"

        nfails = nfails + 1
        call errorc(msg, nfmpi_strerror(err))
        end
