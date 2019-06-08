/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncadios_driver.h>

static PNC_driver ncadios_driver = {
    /* FILE APIs */
    ncadios_create,
    ncadios_open,
    ncadios_close,
    ncadios_enddef,
    ncadios__enddef,
    ncadios_redef,
    ncadios_sync,
    ncadios_flush,
    ncadios_abort,
    ncadios_set_fill,
    ncadios_inq,
    ncadios_inq_misc,
    ncadios_sync_numrecs,
    ncadios_begin_indep_data,
    ncadios_end_indep_data,

    /* DIMENSION APIs */
    ncadios_def_dim,
    ncadios_inq_dimid,
    ncadios_inq_dim,
    ncadios_rename_dim,

    /* ATTRIBUTE APIs */
    ncadios_inq_att,
    ncadios_inq_attid,
    ncadios_inq_attname,
    ncadios_copy_att,
    ncadios_rename_att,
    ncadios_del_att,
    ncadios_get_att,
    ncadios_put_att,

    /* VARIABLE APIs */
    ncadios_def_var,
    ncadios_def_var_fill,
    ncadios_fill_var_rec,
    ncadios_inq_var,
    ncadios_inq_varid,
    ncadios_rename_var,
    ncadios_get_var,
    ncadios_put_var,
    ncadios_get_varn,
    ncadios_put_varn,
    ncadios_get_vard,
    ncadios_put_vard,
    ncadios_iget_var,
    ncadios_iput_var,
    ncadios_bput_var,
    ncadios_iget_varn,
    ncadios_iput_varn,
    ncadios_bput_varn,

    ncadios_buffer_attach,
    ncadios_buffer_detach,
    ncadios_wait,
    ncadios_cancel
};

PNC_driver* ncadios_inq_driver(void) {
    return &ncadios_driver;
}

