!*********************************************************************
!*
!*  Copyright (C) 2010, Northwestern University and Argonne National
!*  Laboratory
!*  See COPYRIGHT notice in top-level directory.
!*
!*********************************************************************/

!*
!*    prints all MPI-IO hints used
!*

        program main
        use mpi
        use pnetcdf
        implicit none

        integer argc, ncid, rank, info, omode, ierr
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
        ierr = nfmpi_open(MPI_COMM_WORLD, trim(filename), omode, &
                          MPI_INFO_NULL, ncid)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_open',ierr)


        ierr = nfmpi_get_file_info(ncid, info)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_get_file_info',ierr)

        ierr = nfmpi_close(ncid)
        if (ierr .ne. NF_NOERR) call handle_err('nfmpi_close',ierr)

        if (rank .EQ. 0) call print_info(info)
        call MPI_Info_free(info, ierr)

 999    call MPI_Finalize(ierr)

        end program main

        subroutine print_info(info_used)
            use mpi
            use pnetcdf
            implicit none

            integer, intent(in) :: info_used

            ! local variables
            character*(MPI_MAX_INFO_VAL) key, value
            integer nkeys, i, ierr
            logical flag

            call MPI_Info_get_nkeys(info_used, nkeys, ierr)
            print *, 'MPI File Info: nkeys =', nkeys
            do i=0, nkeys-1
                call MPI_Info_get_nthkey(info_used, i, key, ierr)
                call MPI_Info_get(info_used, key, MPI_MAX_INFO_VAL, &
                                  value, flag, ierr)
                123 format('MPI File Info: [',I2,'] key = ',A24, &
                           ', value =',A)
                print 123, i, trim(key), trim(value)
            enddo
            print *, ''

            return
        end subroutine print_info

        subroutine handle_err(err_msg, errcode)
            use mpi
            use pnetcdf
            implicit none

            character*(*), intent(in) :: err_msg
            integer,       intent(in) :: errcode

            ! local variables
            integer ierr 

            print *, 'Error: ',trim(err_msg),' ',nfmpi_strerror(errcode)
            call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
            return
        end subroutine handle_err

