/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncncio_driver.h>

static PNC_driver ncncio_driver = {
    /* FILE APIs */
    ncncio_create,
    ncncio_open,
    ncncio_close,
    ncncio_enddef,
    ncncio__enddef,
    ncncio_redef,
    ncncio_sync,
    ncncio_abort,
    ncncio_set_fill,
    ncncio_inq,
    ncncio_inq_misc,
    ncncio_sync_numrecs,
    ncncio_begin_indep_data,
    ncncio_end_indep_data,

    /* DIMENSION APIs */
    ncncio_def_dim,
    ncncio_inq_dimid,
    ncncio_inq_dim,
    ncncio_rename_dim,

    /* ATTRIBUTE APIs */
    ncncio_inq_att,
    ncncio_inq_attid,
    ncncio_inq_attname,
    ncncio_copy_att,
    ncncio_rename_att,
    ncncio_del_att,
    ncncio_get_att,
    ncncio_put_att,

    /* VARIABLE APIs */
    ncncio_def_var,
    ncncio_def_var_fill,
    ncncio_fill_var_rec,
    ncncio_inq_var,
    ncncio_inq_varid,
    ncncio_rename_var,
    ncncio_get_var,
    ncncio_put_var,
    ncncio_get_varn,
    ncncio_put_varn,
    ncncio_get_vard,
    ncncio_put_vard,
    ncncio_iget_var,
    ncncio_iput_var,
    ncncio_bput_var,
    ncncio_iget_varn,
    ncncio_iput_varn,
    ncncio_bput_varn,

    ncncio_buffer_attach,
    ncncio_buffer_detach,
    ncncio_wait,
    ncncio_cancel
};

PNC_driver* ncncio_inq_driver(void) {
    return &ncncio_driver;
}

