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
    ncmpii_inq_misc,
    ncmpii_sync_numrecs,
    ncmpii_begin_indep_data,
    ncmpii_end_indep_data,

    ncmpii_def_dim,
    ncmpii_inq_dimid,
    ncmpii_inq_dim,
    ncmpii_rename_dim,

    ncmpii_inq_att,
    ncmpii_inq_attid,
    ncmpii_inq_attname,
    ncmpii_copy_att,
    ncmpii_rename_att,
    ncmpii_del_att,
    ncmpii_get_att,
    ncmpii_put_att,

    ncmpii_def_var,
    ncmpii_def_var_fill,
    ncmpii_fill_var_rec,
    ncmpii_inq_var,
    ncmpii_inq_varid,
    ncmpii_rename_var,
    ncmpii_get_var,
    ncmpii_put_var,
    ncmpii_get_varn,
    ncmpii_put_varn,
    ncmpii_get_vard,
    ncmpii_put_vard,
    ncmpii_iget_var,
    ncmpii_iput_var,
    ncmpii_bput_var,
    ncmpii_iget_varn,
    ncmpii_iput_varn,
    ncmpii_bput_varn,

    ncmpii_buffer_attach,
    ncmpii_buffer_detach,
    ncmpii_wait,
    ncmpii_cancel
};

PNC_Dispatch* ncmpii_inq_dispatcher(void) {
    return &ncmpi_dispatcher;
}

