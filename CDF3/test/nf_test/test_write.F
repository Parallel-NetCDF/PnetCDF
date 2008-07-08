C********************************************************************
C   Copyright 1996, UCAR/Unidata
C   See netcdf/COPYRIGHT file for copying and redistribution conditions.
C   $Id: test_write.F 559 2007-08-06 21:32:51Z robl $
C********************************************************************


C Test nfmpi_create
C    For mode in NF_NOCLOBBER, NF_CLOBBER do:
C       create netcdf file 'scratch.nc' with no data, close it
C       test that it can be opened, do nfmpi_inq to check nvars = 0, etc.
C    Try again in NF_NOCLOBBER mode, check error return
C On exit, delete this file
        subroutine test_nfmpi_create()
        implicit        none
#include "tests.inc"

        integer clobber         !/* 0 for NF_NOCLOBBER, 1 for NF_CLOBBER */
        integer err
        integer ncid
        integer ndims           !/* number of dimensions */
        integer nvars           !/* number of variables */
        integer ngatts          !/* number of global attributes */
        integer recdim          !/* id of unlimited dimension */
        integer flags

        flags = NF_NOCLOBBER
        do 1, clobber = 0, 1
            err = nfmpi_create(comm, scratch, flags,  MPI_INFO_NULL,
     +                         ncid)
            if (err .ne. 0) then
                call errore('nfmpi_create: ', err)
            end if
            err = nfmpi_close(ncid)
            if (err .ne. 0) then
                call errore('nfmpi_close: ', err)
            end if
            err = nfmpi_open(comm, scratch, NF_NOWRITE, MPI_INFO_NULL, 
     +                       ncid)
            if (err .ne. 0) then
                call errore('nfmpi_open: ', err)
            end if
            err = nfmpi_inq(ncid, ndims, nvars, ngatts, recdim)
            if (err .ne. 0) then
                call errore('nfmpi_inq: ', err)
            else if (ndims .ne. 0) then
                call errori(
     +              'nfmpi_inq: wrong number of dimensions returned, ',
     +              ndims)
            else if (nvars .ne. 0) then
                call errori(
     +              'nfmpi_inq: wrong number of variables returned, ',
     +              nvars)
            else if (ngatts .ne. 0) then
                call errori(
     +              'nfmpi_inq: wrong number of global atts returned, ',
     +              ngatts)
            else if (recdim .ge. 1) then
                call errori(
     +              'nfmpi_inq: wrong record dimension ID returned, ',
     +              recdim)
            end if
            err = nfmpi_close(ncid)
            if (err .ne. 0) then
                call errore('nfmpi_close: ', err)
            end if

            flags = NF_CLOBBER
1       continue

        err = nfmpi_create(comm, scratch, NF_NOCLOBBER,  MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. NF_EEXIST) then
            call errore('attempt to overwrite file: ', err)
        end if
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0) then
            call errori('delete of scratch file failed: ', err)
        end if
        end


C Test nfmpi_redef 
C (In fact also tests nfmpi_enddef - called from test_nfmpi_enddef)
C    BAD_ID
C    attempt redef (error) & enddef on read-only file
C    create file, define dims & vars. 
C    attempt put var (error)
C    attempt redef (error) & enddef.
C    put vars
C    attempt def new dims (error)
C    redef
C    def new dims, vars.
C    put atts
C    enddef
C    put vars
C    close
C    check file: vars & atts
        subroutine test_nfmpi_redef()
        implicit        none
#include "tests.inc"
        integer         title_len
        parameter       (title_len = 9)

        integer                 ncid            !/* netcdf id */
        integer                 dimid           !/* dimension id */
        integer                 vid             !/* variable id */
        integer                 err
        character*(title_len)   title
        doubleprecision         var
        character*(NF_MAX_NAME) name
        NFMPI_OFFSET                 length
        integer                 intlen

        title = 'Not funny'

C           /* BAD_ID tests */
        err = nfmpi_redef(BAD_ID)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_enddef(BAD_ID)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)

C           /* read-only tests */
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_redef(ncid)
        if (err .ne. NF_EPERM)
     +      call errore('nfmpi_redef in NF_NOWRITE mode: ', err)
        err = nfmpi_enddef(ncid)
        if (err .ne. NF_ENOTINDEFINE)
     +      call errore('nfmpi_redef in NF_NOWRITE mode: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. 0) 
     +      call errore('nfmpi_close: ', err)

C           /* tests using scratch file */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)
        err = nfmpi_inq_varid(ncid, 'd', vid)
        if (err .ne. 0) 
     +      call errore('nfmpi_inq_varid: ', err)
        var = 1.0
        err = nfmpi_begin_indep_data(ncid)       
        err = nfmpi_put_var1_double(ncid, vid, 0, var)
        if (err .ne. NF_EINDEFINE)
     +      call errore('nfmpi_put_var... in define mode: ', err)
        err = nfmpi_end_indep_data(ncid)
        err = nfmpi_redef(ncid)
        if (err .ne. NF_EINDEFINE)
     +      call errore('nfmpi_redef in define mode: ', err)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        call put_vars(ncid)
        length = 8
        err = nfmpi_def_dim(ncid, 'abc', length, dimid)
        if (err .ne. NF_ENOTINDEFINE)
     +      call errore('nfmpi_def_dim in define mode: ', err)
        err = nfmpi_redef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_redef: ', err)
        length = 8
        err = nfmpi_def_dim(ncid, 'abc', length, dimid)
        if (err .ne. 0)
     +      call errore('nfmpi_def_dim: ', err)
        err = nfmpi_def_var(ncid, 'abc', NF_INT, 0, 0, vid)
        if (err .ne. 0)
     +      call errore('nfmpi_def_var: ', err)
        length = len(title)
        err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'title', length, 
     +                        title)
        if (err .ne .0)
     +      call errore('nfmpi_put_att_text: ', err)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        var = 1.0
        err = nfmpi_begin_indep_data(ncid)
        err = nfmpi_put_var1_double(ncid, vid, 0, var)
        if (err .ne. 0)
     +      call errore('nfmpi_put_var1_double: ', err)
        err = nfmpi_end_indep_data(ncid)       
        err = nfmpi_close(ncid)
        if (err .ne. 0) 
     +      call errore('nfmpi_close: ', err)

C           /* check scratch file written as expected */
        call check_file(scratch)
        err = nfmpi_open(comm, scratch, NF_NOWRITE, 
     &        MPI_INFO_NULL, ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_dim(ncid, dimid, name, length)
        if (err .ne. 0) 
     +      call errore('nfmpi_inq_dim: ', err)
        if (name .ne. "abc")
     +      call errori('Unexpected dim name in netCDF ', ncid)
        if (length .ne. 8) then
            intlen = length
            call errori('Unexpected dim length: ', intlen)
        end if
        err = nfmpi_begin_indep_data(ncid)
        err = nfmpi_get_var1_double(ncid, vid, 0, var)
        if (err .ne. 0)
     +      call errore('nfmpi_get_var1_double: ', err)
        if (var .ne. 1.0)
     +      call errori(
     +          'nfmpi_get_var1_double: unexpected value in netCDF '
     +          , ncid)
        err = nfmpi_end_indep_data(ncid)
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete failed for netCDF: ', err)
        end

C Test nfmpi_enddef 
C Simply calls test_nfmpi_redef which tests both nfmpi_redef & nfmpi_enddef

        subroutine test_nfmpi_enddef()
        implicit        none
#include "tests.inc"

        call test_nfmpi_redef
        end


C Test nfmpi_sync
C    try with bad handle, check error
C    try in define mode, check error
C    try writing with one handle, reading with another on same netCDF
        subroutine test_nfmpi_sync()
        implicit        none
#include "tests.inc"

        integer ncidw         !/* netcdf id for writing */
        integer ncidr         !/* netcdf id for reading */
        integer err

C           /* BAD_ID test */
        err = nfmpi_sync(BAD_ID)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)

C           /* create scratch file & try nfmpi_sync in define mode */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncidw)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_sync(ncidw)
        if (err .ne. NF_EINDEFINE)
     +      call errore('nfmpi_sync called in define mode: ', err)

C           /* write using same handle */
        call def_dims(ncidw)
        call def_vars(ncidw)
        call put_atts(ncidw)
        err = nfmpi_enddef(ncidw)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        call put_vars(ncidw)
        err = nfmpi_sync(ncidw)
        if (err .ne. 0)
     +      call errore('nfmpi_sync of ncidw failed: ', err)

C           /* open another handle, nfmpi_sync, read (check) */
        err = nfmpi_open(comm, scratch, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncidr)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_sync(ncidr)
        if (err .ne. 0)
     +      call errore('nfmpi_sync of ncidr failed: ', err)
        call check_dims(ncidr)
        call check_atts(ncidr)
        call check_vars(ncidr)

C           /* close both handles */
        err = nfmpi_close(ncidr)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_close(ncidw)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end


C Test nfmpi_abort
C    try with bad handle, check error
C    try in define mode before anything written, check that file was deleted
C    try after nfmpi_enddef, nfmpi_redef, define new dims, vars, atts
C    try after writing variable
        subroutine test_nfmpi_abort()
        implicit        none
#include "tests.inc"

        integer ncid          !/* netcdf id */
        integer err
        integer ndims
        integer nvars
        integer ngatts
        integer recdim

C           /* BAD_ID test */
        err = nfmpi_abort(BAD_ID)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: status = ', err)

C           /* create scratch file & try nfmpi_abort in define mode */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)
        err = nfmpi_abort(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_abort of ncid failed: ', err)
        err = nfmpi_close(ncid)    !/* should already be closed */
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)    !/* should already be deleted */
        if (err .eq. 0)
     +      call errori('scratch file should not exist: ', err)

C            create scratch file
C            do nfmpi_enddef & nfmpi_redef
C            define new dims, vars, atts
C            try nfmpi_abort: should restore previous state (no dims, vars, atts)
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        err = nfmpi_redef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_redef: ', err)
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)
        err = nfmpi_abort(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_abort of ncid failed: ', err)
        err = nfmpi_close(ncid)    !/* should already be closed */
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_open(comm, scratch, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq (ncid, ndims, nvars, ngatts, recdim)
        if (err .ne. 0)
     +      call errore('nfmpi_inq: ', err)
        if (ndims .ne. 0)
     +      call errori('ndims should be ', 0)
        if (nvars .ne. 0)
     +      call errori('nvars should be ', 0)
        if (ngatts .ne. 0)
     +      call errori('ngatts should be ', 0)
        err = nfmpi_close (ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

C           /* try nfmpi_abort in data mode - should just close */
        err = nfmpi_create(comm, scratch, NF_CLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        call put_vars(ncid)
        err = nfmpi_abort(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_abort of ncid failed: ', err)
        err = nfmpi_close(ncid)       !/* should already be closed */
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        call check_file(scratch)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end


C Test nfmpi_def_dim
C    try with bad netCDF handle, check error
C    try in data mode, check error
C    check that returned id is one more than previous id
C    try adding same dimension twice, check error
C    try with illegal sizes, check error
C    make sure unlimited size works, shows up in nfmpi_inq_unlimdim
C    try to define a second unlimited dimension, check error
        subroutine test_nfmpi_def_dim()
        implicit        none
#include "tests.inc"

        integer ncid
        integer err             !/* status */
        integer i
        integer dimid         !/* dimension id */
        NFMPI_OFFSET length

C           /* BAD_ID test */
        length = 8
        err = nfmpi_def_dim(BAD_ID, 'abc', length, dimid)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)

C           /* data mode test */
        err = nfmpi_create(comm, scratch, NF_CLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        length = 8
        err = nfmpi_def_dim(ncid, 'abc', length, dimid)
        if (err .ne. NF_ENOTINDEFINE)
     +      call errore('bad ncid: ', err)

C           /* define-mode tests: unlimited dim */
        err = nfmpi_redef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_redef: ', err)
        length = NFMPI_UNLIMITED
        err = nfmpi_def_dim(ncid, dim_name(1), length, dimid)
        if (err .ne. 0) 
     +      call errore('nfmpi_def_dim: ', err)
        if (dimid .ne. 1) 
     +      call errori('Unexpected dimid: ', dimid)
        err = nfmpi_inq_unlimdim(ncid, dimid)
        if (err .ne. 0) 
     +      call errore('nfmpi_inq_unlimdim: ', err)
        if (dimid .ne. RECDIM) 
     +      call error('Unexpected recdim: ')
        err = nfmpi_inq_dimlen(ncid, dimid, length)
        if (length .ne. 0) 
     +      call errori('Unexpected length: ', 0)
        length = NFMPI_UNLIMITED
        err = nfmpi_def_dim(ncid, 'abc', length, dimid)
        if (err .ne. NF_EUNLIMIT)
     +      call errore('2nd unlimited dimension: ', err)

C           /* define-mode tests: remaining dims */
        do 1, i = 2, NDIMS
            err = nfmpi_def_dim(ncid, dim_name(i-1), dim_len(i), 
     +                       dimid)
            if (err .ne. NF_ENAMEINUSE)
     +          call errore('duplicate name: ', err)
            err = nfmpi_def_dim(ncid, BAD_NAME, dim_len(i), dimid)
            if (err .ne. NF_EBADNAME)
     +          call errore('bad name: ', err)
            length = NFMPI_UNLIMITED - 1
            err = nfmpi_def_dim(ncid, dim_name(i), length,
     +                       dimid)
            if (err .ne. NF_EDIMSIZE)
     +          call errore('bad size: ', err)
            err = nfmpi_def_dim(ncid, dim_name(i), 
     &               dim_len(i), dimid)
            if (err .ne. 0) 
     +          call errore('nfmpi_def_dim: ', err)
            if (dimid .ne. i) 
     +          call errori('Unexpected dimid: ', 0)
1       continue

C           /* Following just to expand unlimited dim */
        call def_vars(ncid)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        call put_vars(ncid)

C           /* Check all dims */
        call check_dims(ncid)

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end


C Test nfmpi_rename_dim
C    try with bad netCDF handle, check error
C    check that proper rename worked with nfmpi_inq_dim
C    try renaming to existing dimension name, check error
C    try with bad dimension handle, check error
        subroutine test_nfmpi_rename_dim()
        implicit        none
#include "tests.inc"

        integer ncid
        integer err             !/* status */
        character*(NF_MAX_NAME) name

C           /* BAD_ID test */
        err = nfmpi_rename_dim(BAD_ID, 1, 'abc')
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)

C           /* main tests */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        err = nfmpi_rename_dim(ncid, BAD_DIMID, 'abc')
        if (err .ne. NF_EBADDIM)
     +      call errore('bad dimid: ', err)
        err = nfmpi_rename_dim(ncid, 3, 'abc')
        if (err .ne. 0)
     +      call errore('nfmpi_rename_dim: ', err)
        err = nfmpi_inq_dimname(ncid, 3, name)
        if (name .ne. 'abc')
     +      call errorc('Unexpected name: ', name)
        err = nfmpi_rename_dim(ncid, 1, 'abc')
        if (err .ne. NF_ENAMEINUSE)
     +      call errore('duplicate name: ', err)

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end


C Test nfmpi_def_var
C    try with bad netCDF handle, check error
C    try with bad name, check error
C    scalar tests:
C      check that proper define worked with nfmpi_inq_var
C      try redefining an existing variable, check error
C      try with bad datatype, check error
C      try with bad number of dimensions, check error
C      try in data mode, check error
C    check that returned id is one more than previous id
C    try with bad dimension ids, check error
        subroutine test_nfmpi_def_var()
        implicit        none
#include "tests.inc"

        integer ncid
        integer vid
        integer err             !/* status */
        integer i
        integer ndims
        integer na
        character*(NF_MAX_NAME) name
        integer dimids(MAX_RANK)
        integer datatype

C           /* BAD_ID test */
        err = nfmpi_def_var(BAD_ID, 'abc', NF_SHORT, 0, dimids, vid)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: status = ', err)

C           /* scalar tests */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_def_var(ncid, 'abc', NF_SHORT, 0, dimids, vid)
        if (err .ne. 0)
     +      call errore('nfmpi_def_var: ', err)
        err = nfmpi_inq_var(ncid, vid, name, datatype, ndims, dimids, 
     +                   na)
        if (err .ne. 0)
     +      call errore('nfmpi_inq_var: ', err)
        if (name .ne. 'abc')
     +      call errorc('Unexpected name: ', name)
        if (datatype .ne. NF_SHORT)
     +      call error('Unexpected datatype')
        if (ndims .ne. 0)
     +      call error('Unexpected rank')
        err = nfmpi_def_var(ncid, BAD_NAME, NF_SHORT, 
     &         0, dimids, vid)
        if (err .ne. NF_EBADNAME)
     +      call errore('bad name: ', err)
        err = nfmpi_def_var(ncid, 'abc', NF_SHORT, 0, dimids, vid)
        if (err .ne. NF_ENAMEINUSE)
     +      call errore('duplicate name: ', err)
        err = nfmpi_def_var(ncid, 'ABC', BAD_TYPE, -1, dimids, vid)
        if (err .ne. NF_EBADTYPE)
     +      call errore('bad type: ', err)
        err = nfmpi_def_var(ncid, 'ABC', NF_SHORT, -1, dimids, vid)
        if (err .ne. NF_EINVAL)
     +      call errore('bad rank: ', err)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        err = nfmpi_def_var(ncid, 'ABC', NF_SHORT, 0, dimids, vid)
        if (err .ne. NF_ENOTINDEFINE)
     +      call errore('nfmpi_def_var called in data mode: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errorc('delete of scratch file failed: ', scratch)

C           /* general tests using global vars */
        err = nfmpi_create(comm, scratch, NF_CLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        do 1, i = 1, NVARS
            err = nfmpi_def_var(ncid, var_name(i), var_type(i), 
     +                       var_rank(i), var_dimid(1,i), vid)
            if (err .ne. 0) 
     +          call errore('nfmpi_def_var: ', err)
            if (vid .ne. i)
     +          call error('Unexpected varid')
1       continue

C           /* try bad dim ids */
        dimids(1) = BAD_DIMID
        err = nfmpi_def_var(ncid, 'abc', NF_SHORT, 1, dimids, vid)
        if (err .ne. NF_EBADDIM)
     +      call errore('bad dim ids: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errorc('delete of scratch file failed: ', scratch)
        end


C Test nfmpi_rename_var
C    try with bad netCDF handle, check error
C    try with bad variable handle, check error
C    try renaming to existing variable name, check error
C    check that proper rename worked with nfmpi_inq_varid
C    try in data mode, check error
        subroutine test_nfmpi_rename_var()
        implicit        none
#include "tests.inc"

        integer ncid
        integer vid
        integer err
        integer i
        character*(NF_MAX_NAME) name

        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_rename_var(ncid, BAD_VARID, 'newName')
        if (err .ne. NF_ENOTVAR)
     +      call errore('bad var id: ', err)
        call def_dims(ncid)
        call def_vars(ncid)

C           /* Prefix "new_" to each name */
        do 1, i = 1, NVARS
            err = nfmpi_rename_var(BAD_ID, i, 'newName')
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_rename_var(ncid, i, var_name(NVARS))
            if (err .ne. NF_ENAMEINUSE)
     +          call errore('duplicate name: ', err)
            name = 'new_' // var_name(i)
            err = nfmpi_rename_var(ncid, i, name)
            if (err .ne. 0)
     +          call errore('nfmpi_rename_var: ', err)
            err = nfmpi_inq_varid(ncid, name, vid)
            if (err .ne. 0)
     +          call errore('nfmpi_inq_varid: ', err)
            if (vid .ne. i)
     +          call error('Unexpected varid')
1       continue

C           /* Change to data mode */
C           /* Try making names even longer. Then restore original names */
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        do 2, i = 1, NVARS
            name = 'even_longer_' // var_name(i)
            err = nfmpi_rename_var(ncid, i, name)
            if (err .ne. NF_ENOTINDEFINE)
     +          call errore('longer name in data mode: ', err)
            err = nfmpi_rename_var(ncid, i, var_name(i))
            if (err .ne. 0)
     +          call errore('nfmpi_rename_var: ', err)
            err = nfmpi_inq_varid(ncid, var_name(i), vid)
            if (err .ne. 0)
     +          call errore('nfmpi_inq_varid: ', err)
            if (vid .ne. i)
     +          call error('Unexpected varid')
2       continue

        call put_vars(ncid)
        call check_vars(ncid)

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errorc('delete of scratch file failed: ', scratch)
        end


C Test nfmpi_copy_att
C    try with bad source or target netCDF handles, check error
C    try with bad source or target variable handle, check error
C    try with nonexisting attribute, check error
C    check that NF_GLOBAL variable for source or target works
C    check that new attribute put works with target in define mode
C    check that old attribute put works with target in data mode
C    check that changing type and length of an attribute work OK
C    try with same ncid for source and target, different variables
C    try with same ncid for source and target, same variable
        subroutine test_nfmpi_copy_att()
        implicit        none
#include "tests.inc"

        integer ncid_in
        integer ncid_out
        integer vid
        integer err
        integer i
        integer j
        character*(NF_MAX_NAME) name    !/* of att */
        integer datatype                !/* of att */
        NFMPI_OFFSET length                  !/* of att */
        character*1     value

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid_in)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid_out)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid_out)
        call def_vars(ncid_out)

        do 1, i = 0, NVARS
            vid = VARID(i)
            do 2, j = 1, NATTS(i)
                name = ATT_NAME(j,i)
                err = nfmpi_copy_att(ncid_in, BAD_VARID, name,
     +                                ncid_out, vid)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad var id: ', err)
                err = nfmpi_copy_att(ncid_in, vid, name, ncid_out, 
     +                            BAD_VARID)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad var id: ', err)
                err = nfmpi_copy_att(BAD_ID, vid, name, 
     +                ncid_out, vid)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_copy_att(ncid_in, vid, name, 
     +                BAD_ID, vid)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_copy_att(ncid_in, vid, 'noSuch',
     +                                ncid_out, vid)
                if (err .ne. NF_ENOTATT)
     +              call errore('bad attname: ', err)
                err = nfmpi_copy_att(ncid_in, vid, name, 
     +                 ncid_out, vid)
                if (err .ne. 0)
     +              call errore('nfmpi_copy_att: ', err)
                err = nfmpi_copy_att(ncid_out, vid, name,
     +                                ncid_out, vid)
                if (err .ne. 0)
     +              call errore('source = target: ', err)
2           continue
1       continue

        err = nfmpi_close(ncid_in)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

C           /* Close scratch. Reopen & check attributes */
        err = nfmpi_close(ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_open(comm, scratch, NF_WRITE, MPI_INFO_NULL,
     +                   ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        call check_atts(ncid_out)

C           change to define mode
C           define single char. global att. ':a' with value 'A'
C           This will be used as source for following copies
        err = nfmpi_redef(ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_redef: ', err)
        length = 1
        err = nfmpi_put_att_text(ncid_out, NF_GLOBAL, 'a', length, 'A')
        if (err .ne. 0)
     +      call errore('nfmpi_put_att_text: ', err)

C           change to data mode
C           Use scratch as both source & dest.
C           try copy to existing att. change type & decrease length
C           rename 1st existing att of each var (if any) 'a'
C           if this att. exists them copy ':a' to it
        err = nfmpi_enddef(ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        do 3, i = 1, NVARS
            if (NATTS(i) .gt. 0 .and. ATT_LEN(1,i) .gt. 0) then
                err = nfmpi_rename_att(ncid_out, i, 
     +                att_name(1,i), 'a')
                if (err .ne. 0)
     +              call errore('nfmpi_rename_att: ', err)
                err = nfmpi_copy_att(ncid_out, NF_GLOBAL, 'a',
     +                                ncid_out, i)
                if (err .ne. 0)
     +              call errore('nfmpi_copy_att: ', err)
            end if
3       continue
        err = nfmpi_close(ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

C           /* Reopen & check */
        err = nfmpi_open(comm, scratch, NF_WRITE, MPI_INFO_NULL,
     +                   ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 4, i = 1, NVARS
            if (NATTS(i) .gt. 0 .and. ATT_LEN(1,i) .gt. 0) then
                err = nfmpi_inq_att(ncid_out, i, 'a', 
     +                datatype, length)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_att: ', err)
                if (datatype .ne. NF_CHAR)
     +              call error('Unexpected type')
                if (length .ne. 1)
     +              call error('Unexpected length')
                err = nfmpi_get_att_text(ncid_out, i, 'a', value)
                if (err .ne. 0)
     +              call errore('nfmpi_get_att_text: ', err)
                if (value .ne. 'A')
     +              call error('Unexpected value')
            end if                                                   
4       continue                                                   

        err = nfmpi_close(ncid_out)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errorc('delete of scratch file failed', scratch)
        end


C Test nfmpi_rename_att
C    try with bad netCDF handle, check error
C    try with bad variable handle, check error
C    try with nonexisting att name, check error
C    try renaming to existing att name, check error
C    check that proper rename worked with nfmpi_inq_attid
C    try in data mode, check error
        subroutine test_nfmpi_rename_att()
        implicit        none
#include "tests.inc"

        integer ncid
        integer vid
        integer err
        integer i
        integer j
        integer  k
        integer attnum
        character*(NF_MAX_NAME) atnam
        character*(NF_MAX_NAME) name
        character*(NF_MAX_NAME) oldname
        character*(NF_MAX_NAME) newname
        integer nok             !/* count of valid comparisons */
        integer datatype
        integer attyp
        NFMPI_OFFSET length
        NFMPI_OFFSET attlength
        NFMPI_OFFSET ndx(1)
        character*(MAX_NELS)    text
        doubleprecision value(MAX_NELS)
        doubleprecision expect

        nok = 0

        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_rename_att(ncid, BAD_VARID, 'abc', 'newName')
        if (err .ne. NF_ENOTVAR)
     +      call errore('bad var id: ', err)
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)

        do 1, i = 0, NVARS
            vid = VARID(i)
            do 2, j = 1, NATTS(i)
                atnam = ATT_NAME(j,i)
                err = nfmpi_rename_att(BAD_ID, vid, atnam, 
     +                 'newName')
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_rename_att(ncid, vid, 'noSuch', 
     +                  'newName')
                if (err .ne. NF_ENOTATT)
     +              call errore('bad attname: ', err)
                newname = 'new_' // atnam
                err = nfmpi_rename_att(ncid, vid, atnam, newname)
                if (err .ne. 0)
     +              call errore('nfmpi_rename_att: ', err)
                err = nfmpi_inq_attid(ncid, vid, newname, attnum)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_attid: ', err)
                if (attnum .ne. j)
     +              call error('Unexpected attnum')
2           continue
1       continue

C           /* Close. Reopen & check */
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_open(comm, scratch, NF_WRITE, MPI_INFO_NULL, 
     +      ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 3, i = 0, NVARS
            vid = VARID(i)
            do 4, j = 1, NATTS(i)
                atnam = ATT_NAME(j,i)
                attyp = ATT_TYPE(j,i)
                attlength = ATT_LEN(j,i)
                newname = 'new_' // atnam
                err = nfmpi_inq_attname(ncid, vid, j, name)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_attname: ', err)
                if (name .ne. newname)
     +              call error('nfmpi_inq_attname: unexpected name')
                err = nfmpi_inq_att(ncid, vid, name, 
     +              datatype, length)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_att: ', err)
                if (datatype .ne. attyp)
     +              call error('nfmpi_inq_att: unexpected type')
                if (length .ne. attlength)
     +              call error('nfmpi_inq_att: unexpected length')
                if (datatype .eq. NF_CHAR) then
                    err = nfmpi_get_att_text(ncid, vid, name, text)
                    if (err .ne. 0)
     +                  call errore('nfmpi_get_att_text: ', err)
                    do 5, k = 1, attlength
                        ndx(1) = k
                        expect = hash(datatype, -1, ndx)
                        if (ichar(text(k:k)) .ne. expect) then
                            call error(
     +                          'nfmpi_get_att_text: unexpected value')
                        else
                            nok = nok + 1
                        end if
5                   continue
                else
                    err = nfmpi_get_att_double(ncid, vid, name, 
     +                     value)
                    if (err .ne. 0)
     +                  call errore('nfmpi_get_att_double: ', err)
                    do 6, k = 1, attlength
                        ndx(1) = k
                        expect = hash(datatype, -1, ndx)
                        if (inRange(expect, datatype)) then
                            if (.not. equal(value(k),expect,datatype,
     +                                      NF_DOUBLE)) then
                                call error(
     +                        'nfmpi_get_att_double: unexpected value')
                            else
                                nok = nok + 1
                            end if
                        end if
6                   continue
                end if
4           continue
3       continue
        call print_nok(nok)

C           /* Now in data mode */
C           /* Try making names even longer. Then restore original names */

        do 7, i = 0, NVARS
            vid = VARID(i)
            do 8, j = 1, NATTS(i)
                atnam = ATT_NAME(j,i)
                oldname = 'new_' // atnam
                newname = 'even_longer_' // atnam
                err = nfmpi_rename_att(ncid, vid, oldname, newname)
                if (err .ne. NF_ENOTINDEFINE)
     +              call errore('longer name in data mode: ', err)
                err = nfmpi_rename_att(ncid, vid, oldname, atnam)
                if (err .ne. 0)
     +              call errore('nfmpi_rename_att: ', err)
                err = nfmpi_inq_attid(ncid, vid, atnam, attnum)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_attid: ', err)
                if (attnum .ne. j)
     +              call error('Unexpected attnum')
8           continue
7       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end


C Test nfmpi_del_att
C    try with bad netCDF handle, check error
C    try with bad variable handle, check error
C    try with nonexisting att name, check error
C    check that proper delete worked using:
C      nfmpi_inq_attid, nfmpi_inq_natts, nfmpi_inq_varnatts
        subroutine test_nfmpi_del_att()
        implicit        none
#include "tests.inc"

        integer ncid
        integer err
        integer i
        integer j
        integer attnum
        integer na
        integer numatts
        integer vid
        character*(NF_MAX_NAME)  name           !/* of att */

        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        err = nfmpi_del_att(ncid, BAD_VARID, 'abc')
        if (err .ne. NF_ENOTVAR)
     +      call errore('bad var id: ', err)
        call def_dims(ncid)
        call def_vars(ncid)
        call put_atts(ncid)

        do 1, i = 0, NVARS
            vid = VARID(i)
            numatts = NATTS(i)
            do 2, j = 1, numatts
                name = ATT_NAME(j,i)
                err = nfmpi_del_att(BAD_ID, vid, name)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_del_att(ncid, vid, 'noSuch')
                if (err .ne. NF_ENOTATT)
     +              call errore('bad attname: ', err)
                err = nfmpi_del_att(ncid, vid, name)
                if (err .ne. 0)
     +              call errore('nfmpi_del_att: ', err)
                err = nfmpi_inq_attid(ncid, vid, name, attnum)
                if (err .ne. NF_ENOTATT)
     +              call errore('bad attname: ', err)
                if (i .lt. 1) then
                    err = nfmpi_inq_natts(ncid, na)
                    if (err .ne. 0)
     +                  call errore('nfmpi_inq_natts: ', err)
                if (na .ne. numatts-j) then
                        call errori('natts: expected: ', numatts-j)
                        call errori('natts: got:      ', na)
                    end if
                end if
                        err = nfmpi_inq_varnatts(ncid, vid, na)
                if (err .ne. 0)
     +              call errore('nfmpi_inq_natts: ', err)
                if (na .ne. numatts-j) then
                    call errori('natts: expected: ', numatts-j)
                    call errori('natts: got:      ', na)
                end if
2           continue
1       continue

C           /* Close. Reopen & check no attributes left */
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_open(comm, scratch, NF_WRITE, 
     +     MPI_INFO_NULL, ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_natts(ncid, na)
        if (err .ne. 0)
     +      call errore('nfmpi_inq_natts: ', err)
        if (na .ne. 0)
     +      call errori('natts: expected 0, got ', na)
        do 3, i = 0, NVARS
            vid = VARID(i)
            err = nfmpi_inq_varnatts(ncid, vid, na)
            if (err .ne. 0)
     +          call errore('nfmpi_inq_natts: ', err)
            if (na .ne. 0)
     +          call errori('natts: expected 0, got ', na)
3       continue

C           /* restore attributes. change to data mode. try to delete */
        err = nfmpi_redef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_redef: ', err)
        call put_atts(ncid)
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)

        do 4, i = 0, NVARS
            vid = VARID(i)
            numatts = NATTS(i)
            do 5, j = 1, numatts
                name = ATT_NAME(j,i)
                err = nfmpi_del_att(ncid, vid, name)
                if (err .ne. NF_ENOTINDEFINE)
     +              call errore('in data mode: ', err)
5           continue
4       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
        end

C parallel-netcdf doesn't implement set_fill, so i have not
C    parallel-netcdfified this subroutine
C Test nfmpi_set_fill
C    try with bad netCDF handle, check error
C    try in read-only mode, check error
C    try with bad new_fillmode, check error
C    try in data mode, check error
C    check that proper set to NF_FILL works for record & non-record variables
C    (note that it is not possible to test NF_NOFILL mode!)
C    close file & create again for test using attribute _FillValue
        subroutine test_nfmpi_set_fill()
#if 0   
        implicit none
#include "tests.inc"

        integer ncid
        integer vid
        integer err
        integer i
        integer j
        integer old_fillmode
        integer nok             !/* count of valid comparisons */
        character*1 text
        doubleprecision value
        doubleprecision fill
        NFMPI_OFFSET index(MAX_RANK)
        NFMPI_OFFSET length

        nok = 0
        value = 0

C           /* bad ncid */
        err = nfmpi_set_fill(BAD_ID, NF_NOFILL, old_fillmode)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)

C           /* try in read-only mode */
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_set_fill(ncid, NF_NOFILL, old_fillmode)
        if (err .ne. NF_EPERM)
     +      call errore('read-only: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)

C           /* create scratch */
        err = nfmpi_create(comm, scratch, NF_NOCLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if

C           /* BAD_FILLMODE */
        err = nfmpi_set_fill(ncid, BAD_FILLMODE, old_fillmode)
        if (err .ne. NF_EINVAL)
     +      call errore('bad fillmode: ', err)

C           /* proper calls */
        err = nfmpi_set_fill(ncid, NF_NOFILL, old_fillmode)
        if (err .ne. 0)
     +      call errore('nfmpi_set_fill: ', err)
        if (old_fillmode .ne. NF_FILL)
     +      call errori('Unexpected old fill mode: ', old_fillmode)
        err = nfmpi_set_fill(ncid, NF_FILL, old_fillmode)
        if (err .ne. 0)
     +      call errore('nfmpi_set_fill: ', err)
        if (old_fillmode .ne. NF_NOFILL)
     +      call errori('Unexpected old fill mode: ', old_fillmode)

C           /* define dims & vars */
        call def_dims(ncid)
        call def_vars(ncid)

C           /* Change to data mode. Set fillmode again */
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        err = nfmpi_set_fill(ncid, NF_FILL, old_fillmode)
        if (err .ne. 0)
     +      call errore('nfmpi_set_fill: ', err)
        if (old_fillmode .ne. NF_FILL)
     +      call errori('Unexpected old fill mode: ', old_fillmode)

C       /* Write record number NRECS to force writing of preceding records */
C       /* Assumes variable cr is char vector with UNLIMITED dimension */
        err = nfmpi_inq_varid(ncid, 'cr', vid)
        if (err .ne. 0)
     +      call errore('nfmpi_inq_varid: ', err)
        index(1) = NRECS
        text = char(NF_FILL_CHAR)
        err = nfmpi_put_var1_text(ncid, vid, index, text)
        if (err .ne. 0)
     +      call errore('nfmpi_put_var1_text: ', err)

C           /* get all variables & check all values equal default fill */
        do 1, i = 1, NVARS
            if (var_type(i) .eq. NF_CHAR) then
                fill = NF_FILL_CHAR
            else if (var_type(i) .eq. NF_BYTE) then
                fill = NF_FILL_BYTE
            else if (var_type(i) .eq. NF_SHORT) then
                fill = NF_FILL_SHORT
            else if (var_type(i) .eq. NF_INT) then
                fill = NF_FILL_INT
            else if (var_type(i) .eq. NF_FLOAT) then
                fill = NF_FILL_FLOAT
            else if (var_type(i) .eq. NF_DOUBLE) then
                fill = NF_FILL_DOUBLE
            else
                stop 'test_nfmpi_set_fill(): impossible var_type(i)'
            end if

            do 2, j = 1, var_nels(i)
                err = index2indexes(j, var_rank(i), var_shape(1,i), 
     +                              index)
                if (err .ne. 0)
     +              call error('error in index2indexes()')
                if (var_type(i) .eq. NF_CHAR) then
                    err = nfmpi_get_var1_text(ncid, i, index, text)
                    if (err .ne. 0)
     +                  call errore('nfmpi_get_var1_text failed: ',err)
                    value = ichar(text)
                else
                    err = nfmpi_get_var1_double(ncid, i, index, 
     +                 value)
                    if (err .ne. 0)
     +                  call errore
     +                           ('nfmpi_get_var1_double failed: ',err)
                end if
                if (value .ne. fill .and. 
     +              abs((fill - value)/fill) .gt. 1.0e-9) then
                    call errord('Unexpected fill value: ', value)
                else
                    nok = nok + 1
                end if
2           continue
1       continue

C       /* close scratch & create again for test using attribute _FillValue */
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_create(comm, scratch, NF_CLOBBER, MPI_INFO_NULL,
     +                     ncid)
        if (err .ne. 0) then
            call errore('nfmpi_create: ', err)
            return
        end if
        call def_dims(ncid)
        call def_vars(ncid)

C           /* set _FillValue = 42 for all vars */
        fill = 42
        text = char(int(fill))
        length = 1
        do 3, i = 1, NVARS
            if (var_type(i) .eq. NF_CHAR) then
                err = nfmpi_put_att_text(ncid, i, '_FillValue', length, 
     +             text)
                if (err .ne. 0)
     +              call errore('nfmpi_put_att_text: ', err)
            else
                err = nfmpi_put_att_double(ncid, i, '_FillValue',
     +                                  var_type(i),length,fill)
                if (err .ne. 0)
     +              call errore('nfmpi_put_att_double: ', err)
            end if
3       continue

C           /* data mode. write records */
        err = nfmpi_enddef(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_enddef: ', err)
        index(1) = NRECS
        err = nfmpi_put_var1_text(ncid, vid, index, text)
        if (err .ne. 0)
     +      call errore('nfmpi_put_var1_text: ', err)

C           /* get all variables & check all values equal 42 */
        do 4, i = 1, NVARS
            do 5, j = 1, var_nels(i)
                err = index2indexes(j, var_rank(i), var_shape(1,i), 
     +                              index)
                if (err .ne. 0)
     +              call error('error in index2indexes')
                if (var_type(i) .eq. NF_CHAR) then
                    err = nfmpi_get_var1_text(ncid, i, index, text)
                    if (err .ne. 0)
     +                  call errore('nfmpi_get_var1_text failed: ',err)
                    value = ichar(text)
                else
                    err = nfmpi_get_var1_double(ncid, i, 
     +                     index, value)
                    if (err .ne. 0)
     +                  call errore
     +                          ('nfmpi_get_var1_double failed: ', err)
                end if
                if (value .ne. fill) then
                    call errord(' Value expected: ', fill)
                    call errord(' Value read:     ', value)
                else
                    nok = nok + 1
                end if
5           continue
4       continue
        call print_nok(nok)

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        err = nfmpi_delete(scratch, MPI_INFO_NULL)
        if (err .ne. 0)
     +      call errori('delete of scratch file failed: ', err)
#endif
        end

#if 0
C * Test nc_set_default_format
C *    try with bad default format
C *    try with NULL old_formatp
C *    try in data mode, check error
C *    check that proper set to NC_FILL works for record & non-record variables
C *    (note that it is not possible to test NC_NOFILL mode!)
C *    close file & create again for test using attribute _FillValue
      subroutine test_nfmpi_set_default_format()
      implicit none
#include "tests.inc"
      
      integer ncid
      integer err
      integer i
      integer version
      integer old_format
      integer nfmpi_get_file_version
      
C     /* bad format */
      err = nfmpi_set_default_format(3, old_format)
      IF (err .ne. NF_EINVAL)
     +     call errore("bad default format: status = %d", err)
     
C    /* Cycle through available formats. */
      do 1 i=1, 2
         err = nfmpi_set_default_format(i, old_format)
         if (err .ne. 0) 
     +         call errore("setting classic format: status = %d", err)
         err = nfmpi_create(comm, scratch, NF_CLOBBER, 
     +                MPI_INFO_NULL, ncid)
         if (err .ne. 0) 
     +        call errore("bad nfmpi_create: status = %d", err)
         err = nfmpi_put_att_text(ncid, NF_GLOBAL, "testatt", 
     +        4, "blah")
         if (err .ne. 0) call errore("bad put_att: status = %d", err)
         err = nfmpi_close(ncid)
         if (err .ne. 0) call errore("bad close: status = %d", err)
         err = nfmpi_get_file_version(scratch, version)
         if (err .ne. 0) call errore("bad file version = %d", err)
         if (version .ne. i)
     +        call errore("bad file version = %d", err)
 1    continue

C    /* Remove the left-over file. */
      err = nfmpi_delete(scratch)
      if (err .ne. 0) call errore("remove failed", err)
      end

#endif

C     This function looks in a file for the netCDF magic number.
      integer function nfmpi_get_file_version(path, version)
      implicit none
#include "tests.inc"
      
      character*(*) path
      integer version, iosnum
      character magic*4
      integer ver
      integer f
      parameter (f = 10)

      open(f, file=path, status='OLD', form='UNFORMATTED',
     +     access='DIRECT', recl=4)

C     Assume this is not a netcdf file.
      nfmpi_get_file_version = NF_ENOTNC
      version = 0

C     Read the magic number, the first 4 bytes of the file.
      read(f, rec=1, err = 1) magic

C     If the first three characters are not "CDF" we're done.
      if (index(magic, 'CDF') .eq. 1) then
         ver = ichar(magic(4:4))
         if (ver .eq. 1) then
            version = 1
            nfmpi_get_file_version = NF_NOERR
         elseif (ver .eq. 2) then
            version = 2
            nfmpi_get_file_version = NF_NOERR
         endif
      endif

 1    close(f)
      return
      end


