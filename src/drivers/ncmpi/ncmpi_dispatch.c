/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <dispatch.h>
#include <ncmpi_dispatch.h>
#include <nc.h>

static PNC_Dispatch ncmpi_dispatcher = {
    ncmpii_create,
    ncmpii_open,
    ncmpii_close,
    ncmpii_enddef,
    ncmpii__enddef,
    ncmpii_redef,
    ncmpii_sync,
    ncmpii_abort,
    ncmpii_set_fill,
    ncmpii_inq,
    ncmpii_begin_indep_data,
    ncmpii_end_indep_data
};

PNC_Dispatch* ncmpii_inq_dispatcher(void) {
    return &ncmpi_dispatcher;
}

