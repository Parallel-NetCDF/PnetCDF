/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include <mpi.h>

extern
MPI_Offset ncmpidiff_core(const char  *file1,
                          const char  *file2,
                          MPI_Comm     comm,
                          MPI_Info     info,
                          int          verbose,
                          int          quiet,
                          int          check_header,
                          int          check_variable_list,
                          int          check_entire_file,
                          int          num_vars,
                          char       **var_names,
                          int          check_tolerance,
                          int          first_diff,
                          char        *cmd_opts,
                          double       tolerance_difference,
                          double       tolerance_ratio);
