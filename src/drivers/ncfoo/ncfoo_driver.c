/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncfoo_driver.h>

static PNC_driver ncfoo_driver = {
    /* FILE APIs */
    ncfoo_create,
    ncfoo_open,
    ncfoo_close,
    ncfoo_enddef,
    ncfoo__enddef,
    ncfoo_redef,
    ncfoo_sync,
    ncfoo_flush,
    ncfoo_abort,
    ncfoo_set_fill,
    ncfoo_inq,
    ncfoo_inq_misc,
    ncfoo_sync_numrecs,
    ncfoo_begin_indep_data,
    ncfoo_end_indep_data,

    /* DIMENSION APIs */
    ncfoo_def_dim,
    ncfoo_inq_dimid,
    ncfoo_inq_dim,
    ncfoo_rename_dim,

    /* ATTRIBUTE APIs */
    ncfoo_inq_att,
    ncfoo_inq_attid,
    ncfoo_inq_attname,
    ncfoo_copy_att,
    ncfoo_rename_att,
    ncfoo_del_att,
    ncfoo_get_att,
    ncfoo_put_att,

    /* VARIABLE APIs */
    ncfoo_def_var,
    ncfoo_def_var_fill,
    ncfoo_fill_var_rec,
    ncfoo_inq_var,
    ncfoo_inq_varid,
    ncfoo_rename_var,
    ncfoo_get_var,
    ncfoo_put_var,
    ncfoo_get_varn,
    ncfoo_put_varn,
    ncfoo_get_vard,
    ncfoo_put_vard,
    ncfoo_iget_var,
    ncfoo_iput_var,
    ncfoo_bput_var,
    ncfoo_iget_varn,
    ncfoo_iput_varn,
    ncfoo_bput_varn,

    ncfoo_buffer_attach,
    ncfoo_buffer_detach,
    ncfoo_wait,
    ncfoo_cancel
};

PNC_driver* ncfoo_inq_driver(void) {
    return &ncfoo_driver;
}

