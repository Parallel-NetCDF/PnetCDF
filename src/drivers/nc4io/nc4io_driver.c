/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <nc4io_driver.h>

static PNC_driver nc4io_driver = {
    /* FILE APIs */
    nc4io_create,
    nc4io_open,
    nc4io_close,
    nc4io_enddef,
    nc4io__enddef,
    nc4io_redef,
    nc4io_sync,
    nc4io_flush,
    nc4io_abort,
    nc4io_set_fill,
    nc4io_inq,
    nc4io_inq_misc,
    nc4io_sync_numrecs,
    nc4io_begin_indep_data,
    nc4io_end_indep_data,

    /* DIMENSION APIs */
    nc4io_def_dim,
    nc4io_inq_dimid,
    nc4io_inq_dim,
    nc4io_rename_dim,

    /* ATTRIBUTE APIs */
    nc4io_inq_att,
    nc4io_inq_attid,
    nc4io_inq_attname,
    nc4io_copy_att,
    nc4io_rename_att,
    nc4io_del_att,
    nc4io_get_att,
    nc4io_put_att,

    /* VARIABLE APIs */
    nc4io_def_var,
    nc4io_def_var_fill,
    nc4io_fill_var_rec,
    nc4io_inq_var,
    nc4io_inq_varid,
    nc4io_rename_var,
    nc4io_get_var,
    nc4io_put_var,
    nc4io_get_varn,
    nc4io_put_varn,
    nc4io_get_vard,
    nc4io_put_vard,
    nc4io_iget_var,
    nc4io_iput_var,
    nc4io_bput_var,
    nc4io_iget_varn,
    nc4io_iput_varn,
    nc4io_bput_varn,

    nc4io_buffer_attach,
    nc4io_buffer_detach,
    nc4io_wait,
    nc4io_cancel
};

PNC_driver* nc4io_inq_driver(void) {
    return &nc4io_driver;
}

