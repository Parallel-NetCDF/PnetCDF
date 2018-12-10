/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <nczipio_driver.h>

static PNC_driver nczipio_driver = {
    /* FILE APIs */
    nczipio_create,
    nczipio_open,
    nczipio_close,
    nczipio_enddef,
    nczipio__enddef,
    nczipio_redef,
    nczipio_sync,
    nczipio_flush,
    nczipio_abort,
    nczipio_set_fill,
    nczipio_inq,
    nczipio_inq_misc,
    nczipio_sync_numrecs,
    nczipio_begin_indep_data,
    nczipio_end_indep_data,

    /* DIMENSION APIs */
    nczipio_def_dim,
    nczipio_inq_dimid,
    nczipio_inq_dim,
    nczipio_rename_dim,

    /* ATTRIBUTE APIs */
    nczipio_inq_att,
    nczipio_inq_attid,
    nczipio_inq_attname,
    nczipio_copy_att,
    nczipio_rename_att,
    nczipio_del_att,
    nczipio_get_att,
    nczipio_put_att,

    /* VARIABLE APIs */
    nczipio_def_var,
    nczipio_def_var_fill,
    nczipio_fill_var_rec,
    nczipio_inq_var,
    nczipio_inq_varid,
    nczipio_rename_var,
    nczipio_get_var,
    nczipio_put_var,
    nczipio_get_varn,
    nczipio_put_varn,
    nczipio_get_vard,
    nczipio_put_vard,
    nczipio_iget_var,
    nczipio_iput_var,
    nczipio_bput_var,
    nczipio_iget_varn,
    nczipio_iput_varn,
    nczipio_bput_varn,

    nczipio_buffer_attach,
    nczipio_buffer_detach,
    nczipio_wait,
    nczipio_cancel
};

PNC_driver* nczipio_inq_driver(void) {
    return &nczipio_driver;
}

