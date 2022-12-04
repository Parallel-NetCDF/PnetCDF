/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncchkio_driver.h>

static PNC_driver ncchkio_driver = {
    /* FILE APIs */
    ncchkio_create,
    ncchkio_open,
    ncchkio_close,
    ncchkio_enddef,
    ncchkio__enddef,
    ncchkio_redef,
    ncchkio_sync,
    ncchkio_flush,
    ncchkio_abort,
    ncchkio_set_fill,
    ncchkio_inq,
    ncchkio_inq_misc,
    ncchkio_sync_numrecs,
    ncchkio_begin_indep_data,
    ncchkio_end_indep_data,

    /* DIMENSION APIs */
    ncchkio_def_dim,
    ncchkio_inq_dimid,
    ncchkio_inq_dim,
    ncchkio_rename_dim,

    /* ATTRIBUTE APIs */
    ncchkio_inq_att,
    ncchkio_inq_attid,
    ncchkio_inq_attname,
    ncchkio_copy_att,
    ncchkio_rename_att,
    ncchkio_del_att,
    ncchkio_get_att,
    ncchkio_put_att,

    /* VARIABLE APIs */
    ncchkio_def_var,
    ncchkio_def_var_fill,
    ncchkio_fill_var_rec,
    ncchkio_inq_var,
    ncchkio_inq_varid,
    ncchkio_rename_var,
    ncchkio_get_var,
    ncchkio_put_var,
    ncchkio_get_varn,
    ncchkio_put_varn,
    ncchkio_get_vard,
    ncchkio_put_vard,
    ncchkio_iget_var,
    ncchkio_iput_var,
    ncchkio_bput_var,
    ncchkio_iget_varn,
    ncchkio_iput_varn,
    ncchkio_bput_varn,

    ncchkio_buffer_attach,
    ncchkio_buffer_detach,
    ncchkio_wait,
    ncchkio_cancel
};

PNC_driver* ncchkio_inq_driver(void) {
    return &ncchkio_driver;
}

