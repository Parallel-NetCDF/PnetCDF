!
!  Copyright (C) 2012, Northwestern University and Argonne National
!  Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id$

!
!    prints all MPI-IO hints used
!
!    Example standard output:
!
!   MPI File Info: nkeys =          18
!   MPI File Info: [ 0] key =            cb_buffer_size, value =16777216
!   MPI File Info: [ 1] key =             romio_cb_read, value =automatic
!   MPI File Info: [ 2] key =            romio_cb_write, value =automatic
!   MPI File Info: [ 3] key =                  cb_nodes, value =1
!   MPI File Info: [ 4] key =         romio_no_indep_rw, value =false
!   MPI File Info: [ 5] key =              romio_cb_pfr, value =disable
!   MPI File Info: [ 6] key =         romio_cb_fr_types, value =aar
!   MPI File Info: [ 7] key =     romio_cb_fr_alignment, value =1
!   MPI File Info: [ 8] key =     romio_cb_ds_threshold, value =0
!   MPI File Info: [ 9] key =         romio_cb_alltoall, value =automatic
!   MPI File Info: [10] key =        ind_rd_buffer_size, value =4194304
!   MPI File Info: [11] key =        ind_wr_buffer_size, value =524288
!   MPI File Info: [12] key =             romio_ds_read, value =automatic
!   MPI File Info: [13] key =            romio_ds_write, value =automatic
!   MPI File Info: [14] key =            cb_config_list, value =*:1
!   MPI File Info: [15] key =      nc_header_align_size, value =0
!   MPI File Info: [16] key =         nc_var_align_size, value =0
!   MPI File Info: [17] key = nc_header_read_chunk_size, value =0


        program main
        implicit none
        include "mpif.h"
        include "pnetcdf.inc"

        integer argc, IARGC, ncid, rank, info, omode, ierr
        character(len = 256) :: filename

        call MPI_Init(ierr)
        call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)

        argc = IARGC()
        if (argc .NE. 1) then
            print *, 'Usage: get_info filename'
            goto 999
        endif
        call getarg(1, filename)

        omode = NF_NOWRITE + NF_64BIT_OFFSET
        ierr = nfmpi_open(MPI_COMM_WORLD, trim(filename), omode,
     +                    MPI_INFO_NULL, ncid)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_open',ierr)


        ierr = nfmpi_get_file_info(ncid, info)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_get_file_info',
     +                                          ierr)

        ierr = nfmpi_close(ncid)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_close',ierr)

        if (rank .EQ. 0) call print_info(info)
        call MPI_Info_free(info, ierr)

 999    call MPI_Finalize(ierr)

        end program main

        subroutine print_info(info_used)
            implicit none
            include "mpif.h"
            include "pnetcdf.inc"

            integer info_used

            ! local variables
            character*(MPI_MAX_INFO_VAL) key, value
            integer nkeys, i, ierr
            logical flag

            call MPI_Info_get_nkeys(info_used, nkeys, ierr)
            print *, 'MPI File Info: nkeys =', nkeys
            do i=0, nkeys-1
                call MPI_Info_get_nthkey(info_used, i, key, ierr)
                call MPI_Info_get(info_used, key, MPI_MAX_INFO_VAL,
     +                            value, flag, ierr)
 123            format('MPI File Info: [',I2,'] key = ',A25,
     +                 ', value =',A)
                print 123, i, trim(key), trim(value)
            enddo
            print *, ''

            return
        end subroutine print_info

        subroutine handle_err(err_msg, errcode)
            implicit none
            include "mpif.h"
            include "pnetcdf.inc"

            character*(*), intent(in) :: err_msg
            integer,       intent(in) :: errcode

            ! local variables
            integer ierr 

            print *, 'Error: ',trim(err_msg),' ',nfmpi_strerror(errcode)
            call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
            return
        end subroutine handle_err

