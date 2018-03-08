/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncmpio_driver.h>
#include "ncmpio_NC.h"

static PNC_driver ncmpio_driver = {
    /* FILE APIs */
    ncmpio_create,
    ncmpio_open,
    ncmpio_close,
    ncmpio_enddef,
    ncmpio__enddef,
    ncmpio_redef,
    ncmpio_sync,
    ncmpio_flush,
    ncmpio_abort,
    ncmpio_set_fill,
    ncmpio_inq,
    ncmpio_inq_misc,
    ncmpio_sync_numrecs,
    ncmpio_begin_indep_data,
    ncmpio_end_indep_data,

    /* DIMENSION APIs */
    ncmpio_def_dim,
    ncmpio_inq_dimid,
    ncmpio_inq_dim,
    ncmpio_rename_dim,

    /* ATTRIBUTE APIs */
    ncmpio_inq_att,
    ncmpio_inq_attid,
    ncmpio_inq_attname,
    ncmpio_copy_att,
    ncmpio_rename_att,
    ncmpio_del_att,
    ncmpio_get_att,
    ncmpio_put_att,

    /* VARIABLE APIs */
    ncmpio_def_var,
    ncmpio_def_var_fill,
    ncmpio_fill_var_rec,
    ncmpio_inq_var,
    ncmpio_inq_varid,
    ncmpio_rename_var,
    ncmpio_get_var,
    ncmpio_put_var,
    ncmpio_get_varn,
    ncmpio_put_varn,
    ncmpio_get_vard,
    ncmpio_put_vard,
    ncmpio_iget_var,
    ncmpio_iput_var,
    ncmpio_bput_var,
    ncmpio_iget_varn,
    ncmpio_iput_varn,
    ncmpio_bput_varn,

    ncmpio_buffer_attach,
    ncmpio_buffer_detach,
    ncmpio_wait,
    ncmpio_cancel
};

PNC_driver* ncmpio_inq_driver(void) {
    return &ncmpio_driver;
}

