C*********************************************************************
C   Copyright 1996, UCAR/Unidata
C   See netcdf/COPYRIGHT file for copying and redistribution conditions.
C   $Id$
C*********************************************************************

C Test nfmpi_strerror.
C    Try on a bad error status.
C    Test for each defined error status.
C
        subroutine test_nfmpi_strerror()
        implicit        none
#include "tests.inc"
        integer         number_of_messages
        parameter       (number_of_messages = 27)

        integer         i
        integer         status(number_of_messages)
        character*80    message
        character*80    msg(number_of_messages)

        data    status(1)  / NF_NOERR/
        data    status(2)  / NF_EBADID /
        data    status(3)  / NF_EEXIST /
        data    status(4)  / NF_EINVAL /
        data    status(5)  / NF_EPERM /
        data    status(6)  / NF_ENOTINDEFINE /
        data    status(7)  / NF_EINDEFINE /
        data    status(8)  / NF_EINVALCOORDS /
        data    status(9)  / NF_EMAXDIMS /
        data    status(10) / NF_ENAMEINUSE /
        data    status(11) / NF_ENOTATT /
        data    status(12) / NF_EMAXATTS /
        data    status(13) / NF_EBADTYPE /
        data    status(14) / NF_EBADDIM /
        data    status(15) / NF_EUNLIMPOS /
        data    status(16) / NF_EMAXVARS /
        data    status(17) / NF_ENOTVAR /
        data    status(18) / NF_EGLOBAL /
        data    status(19) / NF_ENOTNC /
        data    status(20) / NF_ESTS /
        data    status(21) / NF_EMAXNAME /
        data    status(22) / NF_EUNLIMIT /
        data    status(23) / NF_ENORECVARS /
        data    status(24) / NF_ECHAR /
        data    status(25) / NF_EEDGE /
        data    status(26) / NF_ESTRIDE /
        data    status(27) / NF_EBADNAME /

        data    msg(1)  / 'No error' /
        data    msg(2)  / 'Not a netCDF id' /
        data    msg(3)  / 'netCDF file exists && NC_NOCLOBBER' /
        data    msg(4)  / 'Invalid argument' /
        data    msg(5)  / 'Write to read only' /
        data    msg(6)  / 'Operation not allowed in data mode' /
        data    msg(7)  / 'Operation not allowed in define mode' /
        data    msg(8)  / 'Index exceeds dimension bound' /
        data    msg(9)  / 'NC_MAX_DIMS exceeded' /
        data    msg(10) / 'String match to name in use' /
        data    msg(11) / 'Attribute not found' /
        data    msg(12) / 'NC_MAX_ATTRS exceeded' /
        data    msg(13)
     +      / 'Not a netCDF data type or _FillValue type mismatch' /
        data    msg(14) / 'Invalid dimension id or name' /
        data    msg(15) / 'NC_UNLIMITED in the wrong index' /
        data    msg(16) / 'NC_MAX_VARS exceeded' /
        data    msg(17) / 'Variable not found' /
        data    msg(18) / 'Action prohibited on NC_GLOBAL varid' /
        data    msg(19) / 'Not a netCDF file' /
        data    msg(20) / 'In Fortran, string too short' /
        data    msg(21) / 'NC_MAX_NAME exceeded' /
        data    msg(22) / 'NC_UNLIMITED size already in use' /
        data    msg(23) / 'nc_rec op when there are no record vars' /
        data    msg(24) / 'Attempt to convert between text & numbers' /
        data    msg(25) / 'Edge+start exceeds dimension bound' /
        data    msg(26) / 'Illegal stride' /
        data    msg(27) 
     +      / 'Attribute or variable name contains illegal characters' /

C       /* Try on a bad error status */
        message = nfmpi_strerror(-666)!/* should fail */
        if (message .ne. 'Unknown Error')
     +      call errorc('nfmpi_strerror on bad error status returned: ',
     +          message)

C       /* Try on each legitimate error status */
        do 1, i=1, number_of_messages
            message = nfmpi_strerror(status(i))
            if (message .ne. msg(i))
     +          call error('nfmpi_strerror() should return "' 
     +                     // msg(i) // '"')
1       continue
        end


C Test nfmpi_open.
C If in read-only section of tests,
C    Try to open a non-existent netCDF file, check error return.
C    Open a file that is not a netCDF file, check error return.
C    Open a netCDF file with a bad mode argument, check error return.
C    Open a netCDF file with NFMPI_NOWRITE mode, try to write, check error.
C    Try to open a netcdf twice, check whether returned netcdf ids different.
C If in writable section of tests,
C    Open a netCDF file with NFMPI_WRITE mode, write something, close it.
C On exit, any open netCDF files are closed.
        subroutine test_nfmpi_open()
        implicit        none
#include "tests.inc"
        integer err
        integer ncid
        integer ncid2
        
C       /* Try to open a nonexistent file */
        call nfmpi_open(comm, 'tooth-fairy.nc', NF_NOWRITE,
     +                   MPI_INFO_NULL, ncid, err)!/* should fail */
        if (err .eq. NF_NOERR)
     +      call error
     +      ('nfmpi_open of nonexistent file should have failed')
        if (.not. NFMPI_ISSYSERR(err))
     +      call error(
     +  'nfmpi_open of nonexistent file should 
     +   have returned system error')

C       /* Open a file that is not a netCDF file. */
        call nfmpi_open(comm, 'tests.inc', NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)!/* should fail */
        if (err .ne. NF_ENOTNC)
     +      call errore('nfmpi_open of non-netCDF file: ', err)

C       /* Open a netCDF file in read-only mode, check that write fails */
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_redef(ncid, err)    !/* should fail */
        if (err .ne. NF_EPERM)
     +      call error('nfmpi_redef of read-only file should fail')
C       /* Opened OK, see if can open again and get a different netCDF ID */
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid2, err)
        if (err .ne. 0) then
            call errore('nfmpi_open: ', err)
        else
            call nfmpi_close(ncid2, err)
        end if
        if (ncid2 .eq. ncid)
     +      call error(
     +  'netCDF IDs for first and second
     +   nfmpi_open calls should differ')

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            call nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid2, err)
            if (err .ne. 0) then
                call errore('nfmpi_create: ', err)
            else 
                call nfmpi_close(ncid2, err)
            end if
            call nfmpi_open(comm, scratch, NF_WRITE, MPI_INFO_NULL,
     +                       ncid2, err)
            if (err .ne. 0) then
                call errore('nfmpi_open: ', err)
            else 
                call nfmpi_close(ncid2, err)
            end if
            call nfmpi_delete(scratch, err)
            if (err .ne. 0) 
     +          call errorc('delete of scratch file failed: ', scratch)
        end if

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


C 
C Test nfmpi_close.
C    Try to close a netCDF file twice, check whether second close fails.
C    Try on bad handle, check error return.
C    Try in define mode and data mode.
C
        subroutine test_nfmpi_close()
        implicit        none
#include "tests.inc"
        integer ncid
        integer err

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

C       /* Close a netCDF file twice, second time should fail */
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close failed: ', err)
        call nfmpi_close(ncid, err)
        if (err .ne. NF_EBADID)
     +      call error('nfmpi_close of closed file should have failed')
        
C       /* Try with a bad netCDF ID */
        call nfmpi_close(BAD_ID, err)!/* should fail */
        if (err .ne. NF_EBADID)
     +      call errore(
     +         'nfmpi_close with bad netCDF ID returned wrong error: ', 
     +         err)

C       /* Close in data mode */
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close in data mode failed: ', err)

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            call nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid, err)
            if (err .ne. 0) 
     +          call errore('nfmpi_create: ', err)
            call nfmpi_close(ncid, err)
            if (err .ne. 0)
     +          call errore('nfmpi_close in define mode: ', err)
            call nfmpi_delete(scratch, err)
            if (err .ne. 0)
     +          call errorc('delete of scratch file failed: ', 
     +              scratch)
        end if
        end


C Test nfmpi_inq.
C    Try on bad handle, check error return.
C    Try in data mode, check returned values.
C    Try asking for subsets of info.
C If in writable section of tests,
C    Try in define mode, after adding an unlimited dimension, variable.
C On exit, any open netCDF files are closed.
        subroutine test_nfmpi_inq()
        implicit        none
#include "tests.inc"
        integer ncid
        integer ncid2                   !/* for scratch netCDF dataset */
        integer ndims                   !/* number of dimensions */
        integer nvars                   !/* number of variables */
        integer ngatts                  !/* number of global attributes */
        integer recdim                  !/* id of unlimited dimension */
        integer err
        integer ndims0
        integer nvars0
        integer ngatts0
        integer recdim0
        integer did
        integer vid

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        
C       /* Try on bad handle */
        call nfmpi_inq(BAD_ID, ndims, nvars, ngatts, recdim, err)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        
        call nfmpi_inq(ncid, ndims, nvars, ngatts, recdim, err)
        if (err .ne. 0) then
            call errore('nfmpi_inq: ', err)
        else if (ndims .ne. NDIMS) then
            call errori
     +      ('nfmpi_inq: wrong number of dimensions returned: ', ndims)
        else if (nvars .ne. NVARS) then
            call errori
     +      ('nfmpi_inq: wrong number of variables returned: ', nvars)
        else if (ngatts .ne. NGATTS) then
            call errori(
     +          'nfmpi_inq: wrong number of global atts returned: ',
     +          ngatts)
        else if (recdim .ne. RECDIM) then
            call errori
     +      ('nfmpi_inq: wrong record dimension ID returned: ', recdim)
        end if

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            call nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid2, err)
            if (err .ne. 0) then
                call errore('nfmpi_create: ', err)
            else                !/* add dim, var, gatt, check inq */
                call nfmpi_enddef(ncid2, err) !/* enter data mode */
                call nfmpi_inq(ncid2, ndims0, nvars0, 
     +              ngatts0, recdim0, err)
                if (err .ne. 0)
     +              call errore('nfmpi_inq: ', err)
                call nfmpi_redef(ncid2, err) !/* enter define mode */
C               /* Check that inquire still works in define mode */
                call nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim, err)
                if (err .ne. 0) then
                    call errore('nfmpi_inq in define mode: ', err)
                else if (ndims .ne. ndims0) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0) then
                    call errori(
     +              'nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                else if (recdim .ne. recdim0) then
                    call errori
     +              ('nfmpi_inq in define mode: recdim wrong, ', recdim)
                end if

C               /* Add dim, var, global att */
                call nfmpi_def_dim(ncid2, 'inqd', 1, did, err)
                if (err .ne. 0)
     +              call errore('nfmpi_def_dim: ', err)
                call nfmpi_def_var(ncid2, 'inqv', NF_FLOAT, 
     +                0, 0, vid, err)
                if (err .ne. 0)
     +              call errore('nfmpi_def_var: ', err)

                call nfmpi_put_att_text(ncid2, NF_GLOBAL, 'inqa', 
     +                                len('stuff'), 'stuff', err)
                if (err .ne. 0)
     +              call errore('nfmpi_put_att_text: ', err)

C               /* Make sure nfmpi_inq sees the additions while in define mode */
                call nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim, err)
                if (err .ne. 0) then
                    call errore('nfmpi_inq in define mode: ', err)
                else if (ndims .ne. ndims0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                end if
                call nfmpi_enddef(ncid2, err)
                if (err .ne. 0)
     +              call errore('nfmpi_enddef: ', err)

C               /* Make sure nfmpi_inq stills sees additions in data mode */
                call nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim, err)
                if (err .ne. 0) then
                    call errore('nfmpi_inq failed in data mode: ',err)
                else if (ndims .ne. ndims0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                end if
                call nfmpi_close(ncid2, err)
                call nfmpi_delete(scratch, err)
                if (err .ne. 0)
     +              call errorc('delete of scratch file failed: ', 
     +                  scratch)
            end if
        end if

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_natts()
        implicit        none
#include "tests.inc"
        integer ncid
        integer ngatts                  !/* number of global attributes */
        integer err

        call nfmpi_inq_natts(BAD_ID, ngatts, err)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_inq_natts(ncid, ngatts, err)
        if (err .ne. 0) then
            call errore('nfmpi_inq_natts: ', err)
        else if (ngatts .ne. NGATTS) then
            call errori
     +      ('nfmpi_inq_natts: wrong number of global atts returned, ',
     +        ngatts)
        end if
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_ndims()
        implicit        none
#include "tests.inc"
        integer ncid
        integer ndims
        integer err

        call nfmpi_inq_ndims(BAD_ID, ndims, err)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_inq_ndims(ncid, ndims, err)
        if (err .ne. 0) then
            call errore('nfmpi_inq_ndims: ', err)
        else if (ndims .ne. NDIMS) then
            call errori
     +      ('nfmpi_inq_ndims: wrong number returned, ', ndims)
        end if
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_nvars()
        implicit        none
#include "tests.inc"
        integer ncid
        integer nvars
        integer err

        call nfmpi_inq_nvars(BAD_ID, nvars, err)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_inq_nvars(ncid, nvars, err)
        if (err .ne. 0) then
            call errore('nfmpi_inq_nvars: ', err)
        else if (nvars .ne. NVARS) then
            call errori
     +      ('nfmpi_inq_nvars: wrong number returned, ', nvars)
        end if
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_unlimdim()
        implicit        none
#include "tests.inc"
        integer ncid
        integer unlimdim
        integer err

        call nfmpi_inq_unlimdim(BAD_ID, unlimdim, err)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_inq_unlimdim(ncid, unlimdim, err)
        if (err .ne. 0) then
            call errore('nfmpi_inq_unlimdim: ', err)
        else if (unlimdim .ne. RECDIM) then
            call errori
     +      ('nfmpi_inq_unlimdim: wrong number returned, ', unlimdim)
        end if
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimid()
        implicit        none
#include "tests.inc"
        integer ncid
        integer dimid
        integer i
        integer err

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call nfmpi_inq_dimid(ncid, 'noSuch', dimid, err)
        if (err .ne. NF_EBADDIM)
     +      call errore('bad dim name: ', err)
        do 1, i = 1, NDIMS
            call nfmpi_inq_dimid(BAD_ID, dim_name(i), dimid, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_dimid(ncid, dim_name(i), dimid, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimid: ', err)
            else if (dimid .ne. i) then
                call errori('expected ', i)
                call errori('got ', dimid)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dim()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name
        integer length

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            call nfmpi_inq_dim(BAD_ID, i, name, length, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_dim(ncid, BAD_DIMID, name, length, err)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            call nfmpi_inq_dim(ncid, i, name, length, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dim: ', err)
            else if (dim_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
                print *, ' expected ', dim_name(i),' for the ',i,'entry'
            else if (dim_len(i) .ne. length) then
                call errori('size unexpected: ', length)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimlen()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer length

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            call nfmpi_inq_dimlen(BAD_ID, i, length, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_dimlen(ncid, BAD_DIMID, length, err)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            call nfmpi_inq_dimlen(ncid, i, length, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimlen: ', err)
            else if (dim_len(i) .ne. length) then
                call errori('size unexpected: ', length)
                print *, 'expected ', dim_len(i),' for the ',i,'entry'
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimname()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME)  name

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            call nfmpi_inq_dimname(BAD_ID, i, name, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_dimname(ncid, BAD_DIMID, name, err)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            call nfmpi_inq_dimname(ncid, i, name, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimname: ', err)
            else if (dim_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varid()
        implicit        none
#include "tests.inc"
        integer ncid
        integer vid
        integer i
        integer err

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        call nfmpi_inq_varid(ncid, 'noSuch', vid, err)
        if (err .ne. NF_ENOTVAR)
     +      call errore('bad ncid: ', err)

        do 1, i = 1, NVARS
            call nfmpi_inq_varid(BAD_ID, var_name(i), vid, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_varid(ncid, var_name(i), vid, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varid: ', err)
            else if (vid .ne. i) then
                call errori('varid unexpected: ', vid)
            endif
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_var()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name
        integer datatype
        integer ndims
        integer dimids(MAX_RANK)
        integer na

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            call nfmpi_inq_var(BAD_ID, i, name, datatype, ndims, dimids,
     +                       na, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_var(ncid,BAD_VARID,name,datatype,ndims,
     +                       dimids,na, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_var(ncid, i, name, datatype, ndims, dimids, 
     +                       na, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_var: ', err)
            else if (var_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            else if (var_type(i) .ne. datatype) then
                call errori('type unexpected: ', datatype)
            else if (var_rank(i) .ne. ndims) then
                call errori('ndims expected: ', ndims)
            else if (.not.int_vec_eq(var_dimid(1,i),dimids,ndims)) then
                call error('unexpected dimid')
            else if (var_natts(i) .ne. na) then
                call errori('natts unexpected: ', na)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_vardimid()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer dimids(MAX_RANK)

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            call nfmpi_inq_vardimid(BAD_ID, i, dimids, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_vardimid(ncid, BAD_VARID, dimids, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_vardimid(ncid, i, dimids, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_vardimid: ', err)
            else if (.not.int_vec_eq(var_dimid(1,i), dimids, 
     +               var_rank(i))) then
                call error('unexpected dimid')
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varname()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            call nfmpi_inq_varname(BAD_ID, i, name, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_varname(ncid, BAD_VARID, name, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_varname(ncid, i, name, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varname: ', err)
            else if (var_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varnatts()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer na

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 0, NVARS	! start with global attributes
            call nfmpi_inq_varnatts(BAD_ID, i, na, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_varnatts(ncid, BAD_VARID, na, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_varnatts(ncid, VARID(i), na, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varnatts: ', err)
            else if (NATTS(i) .ne. na) then	! works for global attributes
                call errori('natts unexpected: ', na)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varndims()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer ndims

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            call nfmpi_inq_varndims(BAD_ID, i, ndims, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_varndims(ncid, BAD_VARID, ndims, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_varndims(ncid, i, ndims, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varndims: ', err)
            else if (var_rank(i) .ne. ndims) then
                call errori('ndims unexpected: ', ndims)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_vartype()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer datatype

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            call nfmpi_inq_vartype(BAD_ID, i, datatype, err)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            call nfmpi_inq_vartype(ncid, BAD_VARID, datatype, err)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            call nfmpi_inq_vartype(ncid, i, datatype, err)
            if (err .ne. 0) then
                call errore('nfmpi_inq_vartype: ', err)
            else if (var_type(i) .ne. datatype) then
                call errori('type unexpected: ', datatype)
            end if
1       continue
        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_att()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer t
        integer n

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0) 
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            do 2, j = 1, NATTS(i)
                call nfmpi_inq_att(BAD_ID, i, ATT_NAME(j,i), t, n, err)
                if (err .ne. NF_EBADID) 
     +              call errore('bad ncid: ', err)
                call nfmpi_inq_att
     +                     (ncid, BAD_VARID, ATT_NAME(j,i), t, n, err)
                if (err .ne. NF_ENOTVAR) 
     +              call errore('bad var id: ', err)
                call nfmpi_inq_att(ncid, i, 'noSuch', t, n, err)
                if (err .ne. NF_ENOTATT) 
     +              call errore('Bad attribute name: ', err)
                call nfmpi_inq_att(ncid, i, ATT_NAME(j,i), t, n, err)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (t .ne. ATT_TYPE(j,i))
     +                  call error('type not that expected')
                    if (n .ne. ATT_LEN(j,i)) 
     +                  call error('length not that expected')
                end if
2           continue
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attlen()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer len

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            call nfmpi_inq_attlen(ncid, i, 'noSuch', len, err)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                call nfmpi_inq_attlen(BAD_ID, i, 
     +             ATT_NAME(j,i), len, err)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                call nfmpi_inq_attlen
     +                (ncid, BAD_VARID, ATT_NAME(j,i), len, err)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                call nfmpi_inq_attlen(ncid, i, ATT_NAME(j,i), len, err)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (len .ne. ATT_LEN(j,i))
     +                  call error('len not that expected')
                end if
2           continue
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_atttype()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer datatype

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            call nfmpi_inq_atttype(ncid, i, 'noSuch', datatype, err)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                call nfmpi_inq_atttype
     +                (BAD_ID, i, ATT_NAME(j,i), datatype, err)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                call nfmpi_inq_atttype(ncid, BAD_VARID, ATT_NAME(j,i), 
     +                               datatype, err)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                call nfmpi_inq_atttype
     +                (ncid, i, ATT_NAME(j,i), datatype, err)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (datatype .ne. ATT_TYPE(j,i))
     +                  call error('type not that expected')
                end if
2           continue
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attname()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        character*(NF_MAX_NAME) name

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            call nfmpi_inq_attname(ncid, i, BAD_ATTNUM, name, err)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute number: ', err)
            call nfmpi_inq_attname(ncid, i, NATTS(i)+1, name, err)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute number: ', err)
            do 2, j = 1, NATTS(i)
                call nfmpi_inq_attname(BAD_ID, i, j, name, err)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                call nfmpi_inq_attname(ncid, BAD_VARID, j, name, err)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad var id: ', err)
                call nfmpi_inq_attname(ncid, i, j, name, err)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (ATT_NAME(j,i) .ne. name)
     +                  call error('name not that expected')
                end if
2           continue
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attid()
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer attnum

        call nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            call nfmpi_inq_attid(ncid, i, 'noSuch', attnum, err)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                call nfmpi_inq_attid(BAD_ID, i, 
     +                ATT_NAME(j,i), attnum, err)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                call nfmpi_inq_attid(ncid, BAD_VARID, ATT_NAME(j,i), 
     +                             attnum, err)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                call nfmpi_inq_attid(ncid, i, 
     +               ATT_NAME(j,i), attnum, err)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (attnum .ne. j)
     +                  call error('attnum not that expected')
                end if
2           continue
1       continue

        call nfmpi_close(ncid, err)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end
