        module mpinetcdf

        interface

        subroutine ncmpif_create( comm, path, cmode, info, ncid, ierr )
        integer comm, cmode, info, ncid, ierr
        character *(*) path
        end

        subroutine ncmpif_open( comm, path, omode, info, ncid, ierr )
        integer comm, omode, info, ncid, ierr
        character *(*) path
        end

        subroutine ncmpif_enddef( ncid, ierr )
        integer ncid, ierr
        end

        subroutine ncmpif_close( ncid, ierr )
        integer ncid, ierr
        end

! Begin {put,get}_vara

        subroutine ncmpif_put_vara_int_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), op(*), ierr
        end

        
        subroutine ncmpif_put_vara_int( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), op(*), ierr
        end

        subroutine ncmpif_put_vara_float_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        real op(*)
        end

        
        subroutine ncmpif_put_vara_float( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        real op(*)
        end

        subroutine ncmpif_put_vara_double_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        double precision op(*)
        end

        
        subroutine ncmpif_put_vara_double( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        double precision op(*)
        end

        subroutine ncmpif_get_vara_int_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), op(*), ierr
        end

        
        subroutine ncmpif_get_vara_int( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), op(*), ierr
        end

        subroutine ncmpif_get_vara_float_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        real op(*)
        end

        
        subroutine ncmpif_get_vara_float( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        real op(*)
        end

        subroutine ncmpif_get_vara_double_all( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        double precision op(*)
        end

        
        subroutine ncmpif_get_vara_double( ncid, varid, size1, size2, op, ierr )
        integer ncid, varid, size1(*), size2(*), ierr
        double precision op(*)
        end

        end interface
        end module mpinetcdf
