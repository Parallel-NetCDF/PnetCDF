/*
 *  Copyright (C) 2019, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include <ctype.h>

#define CHECK_NAME(A) (strlen(A) > 0 && (isalpha(A[0]) || A[0] == '_'))

extern int
ncadiosi_parse_header_bp2ncd (NC_ad *ncid);

extern int
ncadiosi_parse_header_readall (NC_ad *ncadp);

extern int
ncadiosi_parse_rec_dim(NC_ad *ncadp);

extern int
ncadiosi_var_list_init(NC_ad_var_list *list);

extern int
ncadiosi_dim_list_init(NC_ad_dim_list *list);

extern int
ncadiosi_att_list_init(NC_ad_att_list *list);

extern int
ncadiosi_var_list_free(NC_ad_var_list *list);

extern int
ncadiosi_dim_list_free(NC_ad_dim_list *list);

extern int
ncadiosi_att_list_free(NC_ad_att_list *list);

extern int
ncadiosi_var_list_add(NC_ad_var_list *list, NC_ad_var data);

extern int
ncadiosi_att_list_add(NC_ad_att_list *list, int data);

extern int
ncadiosi_dim_list_add(NC_ad_dim_list *list, NC_ad_dim data);

extern int
ncadiosi_var_list_find(NC_ad_var_list *list, char *name);

extern int
ncadiosi_att_list_find(NC_ad_att_list *list, int data);

extern int
ncadiosi_dim_list_find(NC_ad_dim_list *list, char *name);

extern int
ncadiosi_inq_varid(NC_ad* ncadp, char* name, int *id);

extern int
ncadiosi_inq_attid(NC_ad* ncadp, int vid, char* name, int *id);

extern int
ncadiosi_inq_dimid(NC_ad* ncadp, char* name, int *id);

extern int
ncadiosi_def_var(NC_ad* ncadp, char* name, nc_type type, int ndim, int *dimids, int *id);

extern int
ncadiosi_def_dim(NC_ad* ncadp, char* name, int len, int *id);

extern int
ncadios_sync_header(NC_ad *ncadp);

extern int
ncadiosi_parse_attrs(NC_ad *ncadp);

extern int
ncadiosiconvert(void *inbuf, void *outbuf, MPI_Datatype intype, MPI_Datatype outtype, int N);

extern nc_type
ncadios_to_nc_type(enum ADIOS_DATATYPES atype);

extern MPI_Datatype
ncadios_to_mpi_type(enum ADIOS_DATATYPES atype);

extern int
ncadiosi_get_list_remove(NC_ad_get_list *lp, int reqid);

extern int
ncadiosi_get_list_add(NC_ad_get_list *lp, int *id);

extern int
ncadiosi_get_list_free(NC_ad_get_list *lp);

extern int
ncadiosi_get_list_init(NC_ad_get_list *lp);

extern int
ncadiosi_perform_read(NC_ad *ncadp);

extern int
ncadiosi_handle_get_req(NC_ad *ncadp, NC_ad_get_req *req);

extern int
ncadiosi_wait_get_req(NC_ad *ncadp, int reqid, int *stat);

extern int
ncadiosi_wait_all_get_req(NC_ad *ncadp);

extern int
ncadiosi_init_get_req(NC_ad *ncadp, NC_ad_get_req *r, ADIOS_VARINFO *v, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype);

extern int
ncadiosi_iget_var(NC_ad *ncadp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid);

