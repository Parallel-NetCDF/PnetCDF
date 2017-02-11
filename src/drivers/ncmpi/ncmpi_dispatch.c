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
    ncmpii_get_att_text,
    ncmpii_get_att_schar,
    ncmpii_get_att_uchar,
    ncmpii_get_att_short,
    ncmpii_get_att_ushort,
    ncmpii_get_att_int,
    ncmpii_get_att_uint,
    ncmpii_get_att_long,
    ncmpii_get_att_float,
    ncmpii_get_att_double,
    ncmpii_get_att_longlong,
    ncmpii_get_att_ulonglong,
    ncmpii_put_att,
    ncmpii_put_att_text,
    ncmpii_put_att_schar,
    ncmpii_put_att_uchar,
    ncmpii_put_att_short,
    ncmpii_put_att_ushort,
    ncmpii_put_att_int,
    ncmpii_put_att_uint,
    ncmpii_put_att_long,
    ncmpii_put_att_float,
    ncmpii_put_att_double,
    ncmpii_put_att_longlong,
    ncmpii_put_att_ulonglong,

    ncmpii_def_var,
    ncmpii_def_var_fill,
    ncmpii_inq_var,
    ncmpii_inq_varid,
    ncmpii_rename_var
};

PNC_Dispatch* ncmpii_inq_dispatcher(void) {
    return &ncmpi_dispatcher;
}

