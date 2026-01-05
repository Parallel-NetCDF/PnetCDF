!
!  Copyright (C) 2015, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id$

      subroutine fusage(cmd)
          implicit none
          character(len=*) cmd

          print*,'Usage: ',trim(cmd),' [OPTIONS]'
          print*,'       [-h] Print help'
          print*,'       [-q] quiet mode'
          print*,'       [-k] Keep output files (default: no)'
          print*,'       [-i  in_path]: input file path (default: NULL)'
          print*,'       [-o out_path]: output netCDF file name (default: %s.nc)'
      end subroutine fusage

      ! This function gets the executable name and output file name from the
      ! command line.
      integer function get_args(cmd, out_path, in_path, keep_files)
          implicit none
          character(len=*) cmd, out_path, in_path
          character(len=256) :: full_cmd, arg
          logical :: keep_files, skip_next
          integer :: i, j, n_args

          keep_files = .false.

          get_args = 1
          call getarg(0, full_cmd)

          ! remove basename from executable name
          i = INDEX(full_cmd, "/", .TRUE.)
          if (i .EQ. 0) then
              cmd(:) = full_cmd(:)
          else
              cmd(:) = full_cmd(i+1:)
          endif

          n_args = command_argument_count()

          skip_next = .false.
          do j = 1, n_args
              if (skip_next) then
                  skip_next = .false.
                  cycle
              end if

              call get_command_argument(j, arg)
              arg = trim(arg) ! Remove trailing spaces

              if (arg == "-k") then
                  keep_files = .true.
              else if (arg == "-i") then
                  if (j < n_args) then
                      call get_command_argument(j+1, arg)
                      in_path = trim(arg)
                      skip_next = .true.
                  end if
              else if (arg == "-o") then
                  if (j < n_args) then
                      call get_command_argument(j+1, arg)
                      out_path = trim(arg)
                      skip_next = .true.
                  end if
              else if (arg == "-h") then
                  call fusage(cmd)
                  return
              end if
          end do

      end function get_args

      ! This function prints the pass/fail message on screen
      subroutine pass_fail(nerrs, msg, timing)
          implicit none
          integer nerrs
          character(len=*) msg
          double precision timing

          ! local variables
          CHARACTER ESC
          PARAMETER (ESC=char(27))

#ifdef PNETCDF_DEBUG
          CHARACTER (LEN=20) PASS_STR, FAIL_STR
          PARAMETER (PASS_STR='-- '//ESC//'[32mpass'//ESC//'[0m (')
          PARAMETER (FAIL_STR='-- '//ESC//'[31mfail'//ESC//'[0m')
#else
          CHARACTER (LEN=11) PASS_STR, FAIL_STR
          PARAMETER (PASS_STR='-- pass (')
          PARAMETER (FAIL_STR='-- fail')
#endif

          if (nerrs .EQ. 0) then
              write(*,"(A64,A,F4.1,A)") msg, trim(PASS_STR), timing, 's)'
          else
              write(*,"(A64,A)") msg, trim(FAIL_STR)
          endif
      end subroutine pass_fail

      subroutine get_env(hint_str, value)
          character(len=*) hint_str, value
#ifdef HAS_GET_ENVIRONMENT_VARIABLE
          call Get_Environment_Variable(hint_str, Value=value)
#else
          call getenv(hint_str, value)
#endif
      end subroutine get_env

      LOGICAL FUNCTION relax_coord_bound_f()
          character(len=256) :: env_str, env_val
          integer :: ierr

#ifdef RELAX_COORD_BOUND
          relax_coord_bound_f = .TRUE.
#else
          relax_coord_bound_f = .FALSE.
#endif
          env_str = "PNETCDF_RELAX_COORD_BOUND"
          call get_environment_variable(env_str, value=env_val, status=ierr)

          if (ierr == 0) THEN
              ! Environment variable is set
              if (env_val(1:1) == '1') then
                  relax_coord_bound_f = .TRUE.
              else
                  relax_coord_bound_f = .FALSE.
              endif
          endif
      END FUNCTION relax_coord_bound_f

