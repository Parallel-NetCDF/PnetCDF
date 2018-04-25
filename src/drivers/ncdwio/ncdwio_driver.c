/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncdwio_driver.h>

static PNC_driver ncdwio_driver = {
    /* FILE APIs */
    ncdwio_create,
    ncdwio_open,
    ncdwio_close,
    ncdwio_enddef,
    ncdwio__enddef,
    ncdwio_redef,
    ncdwio_sync,
    ncdwio_abort,
    ncdwio_set_fill,
    ncdwio_inq,
    ncdwio_inq_misc,
    ncdwio_sync_numrecs,
    ncdwio_begin_indep_data,
    ncdwio_end_indep_data,

    /* DIMENSION APIs */
    ncdwio_def_dim,
    ncdwio_inq_dimid,
    ncdwio_inq_dim,
    ncdwio_rename_dim,

    /* ATTRIBUTE APIs */
    ncdwio_inq_att,
    ncdwio_inq_attid,
    ncdwio_inq_attname,
    ncdwio_copy_att,
    ncdwio_rename_att,
    ncdwio_del_att,
    ncdwio_get_att,
    ncdwio_put_att,

    /* VARIABLE APIs */
    ncdwio_def_var,
    ncdwio_def_var_fill,
    ncdwio_fill_var_rec,
    ncdwio_inq_var,
    ncdwio_inq_varid,
    ncdwio_rename_var,
    ncdwio_get_var,
    ncdwio_put_var,
    ncdwio_get_varn,
    ncdwio_put_varn,
    ncdwio_get_vard,
    ncdwio_put_vard,
    ncdwio_iget_var,
    ncdwio_iput_var,
    ncdwio_bput_var,
    ncdwio_iget_varn,
    ncdwio_iput_varn,
    ncdwio_bput_varn,

    ncdwio_buffer_attach,
    ncdwio_buffer_detach,
    ncdwio_wait,
    ncdwio_cancel
};

PNC_driver* ncdwio_inq_driver(void) {
    return &ncdwio_driver;
}

