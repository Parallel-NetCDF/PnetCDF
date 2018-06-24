/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <ncbbio_driver.h>

static PNC_driver ncbbio_driver = {
    /* FILE APIs */
    ncbbio_create,
    ncbbio_open,
    ncbbio_close,
    ncbbio_enddef,
    ncbbio__enddef,
    ncbbio_redef,
    ncbbio_sync,
    ncbbio_flush,
    ncbbio_abort,
    ncbbio_set_fill,
    ncbbio_inq,
    ncbbio_inq_misc,
    ncbbio_sync_numrecs,
    ncbbio_begin_indep_data,
    ncbbio_end_indep_data,

    /* DIMENSION APIs */
    ncbbio_def_dim,
    ncbbio_inq_dimid,
    ncbbio_inq_dim,
    ncbbio_rename_dim,

    /* ATTRIBUTE APIs */
    ncbbio_inq_att,
    ncbbio_inq_attid,
    ncbbio_inq_attname,
    ncbbio_copy_att,
    ncbbio_rename_att,
    ncbbio_del_att,
    ncbbio_get_att,
    ncbbio_put_att,

    /* VARIABLE APIs */
    ncbbio_def_var,
    ncbbio_def_var_fill,
    ncbbio_fill_var_rec,
    ncbbio_inq_var,
    ncbbio_inq_varid,
    ncbbio_rename_var,
    ncbbio_get_var,
    ncbbio_put_var,
    ncbbio_get_varn,
    ncbbio_put_varn,
    ncbbio_get_vard,
    ncbbio_put_vard,
    ncbbio_iget_var,
    ncbbio_iput_var,
    ncbbio_bput_var,
    ncbbio_iget_varn,
    ncbbio_iput_varn,
    ncbbio_bput_varn,

    ncbbio_buffer_attach,
    ncbbio_buffer_detach,
    ncbbio_wait,
    ncbbio_cancel
};

PNC_driver* ncbbio_inq_driver(void) {
    return &ncbbio_driver;
}

