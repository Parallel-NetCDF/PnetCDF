/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef TAU_SSON
#include <TAU.h>
#endif

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_subfile.h"

enum {ONE, BALANCED};
#define SUBFILING_MIN_NDIMS 1

/* set default values for the following values */
static int delegate_scheme = BALANCED; /* default: any proc can be delegate proc */
static int is_partitioned = 0;

#define check_err(fn_name_)                                                 \
    if (mpireturn != MPI_SUCCESS) {                                         \
        errs++;                                                             \
        int _len;                                                           \
        char err_str_[MPI_MAX_ERROR_STRING];                                \
        MPI_Error_string(mpireturn, err_str_, &_len);                       \
        fprintf(stderr, #fn_name_ " failed at line %d, mpireturn=%d: %s\n", \
                __LINE__, mpireturn, err_str_);                             \
    }                                                                       \

#if 0
static int ncmpio_itoa(int val, char* buf)
{
    const unsigned int radix = 10;

    char* p;
    unsigned int a;        /* every digit */
    int len;
    char* b;            /* start of the digit char */
    char temp;
    unsigned int u;

    p = buf;

    if (val < 0) {
        *p++ = '-';
        val = 0 - val;
    }
    u = (unsigned int)val;

    b = p;

    do {
        a = u % radix;
        u /= radix;

        *p++ = a + '0';

    } while (u > 0);

    len = (int)(p - buf);

    *p-- = 0;

    /* swap */
    do {
        temp = *p;
        *p = *b;
        *b = temp;
        --p;
        ++b;
    } while (b < p);

    return len;
}
#endif

/*----< subfile_create() >---------------------------------------------------*/
static int
subfile_create(NC *ncp)
{
    int myrank, nprocs, color, status=NC_NOERR, mpireturn;
    char path_sf[1024];
    double ratio;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Comm_rank(ncp->comm, &myrank);
    MPI_Comm_size(ncp->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): nprocs=%d, ncp->num_subfiles=%d\n",
           __func__, myrank, nprocs, ncp->num_subfiles);
#endif

    /* split the orignial comm to subcomm */
    if (nprocs > ncp->num_subfiles) {
        ratio = (double)nprocs/(double)(ncp->num_subfiles);
        color = (int)((double)myrank/ratio);
    }
    else
        color = myrank%ncp->num_subfiles;

#ifdef SUBFILE_DEBUG
    printf("rank(%d): color=%d\n", myrank, color);
#endif
    /* key = myrank/comm_size; */

    /* TODO: fix error when using generated key value.
     * for now, just passing 0 value. */
    TRACE_COMM(MPI_Comm_split)(ncp->comm, color, myrank, &ncp->comm_sf);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_split");

    sprintf(path_sf, "%s.subfile_%i.%s", ncp->path, color, "nc");

/*
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_lustre_start_iodevice", offset);
    MPI_Info_set(info, "striping_factor", "1");
*/

    void *ncp_sf;
    status = ncmpio_create(ncp->comm_sf, path_sf, ncp->iomode, ncp->ncid, info,
                           &ncp_sf);
    if (status != NC_NOERR && myrank == 0)
        fprintf(stderr, "%s: error in creating file(%s): %s\n",
                __func__, path_sf, ncmpi_strerror(status));

    ncp->ncp_sf = (NC*) ncp_sf;
/*
    MPI_Info_free(&info);
*/

    return status;
}

/*----< ncmpio_subfile_open() >----------------------------------------------*/
int
ncmpio_subfile_open(NC *ncp)
{
    int myrank, nprocs, color, status=NC_NOERR, mpireturn;
    char path_sf[1024];
    double ratio;

    MPI_Comm_rank(ncp->comm, &myrank);
    MPI_Comm_size(ncp->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): nprocs=%d, ncp->num_subfiles=%d\n", __func__,
           myrank, nprocs, ncp->num_subfiles);
#endif

    /* split the original comm to subcomm */
    if (nprocs > ncp->num_subfiles) {
        ratio = (double)nprocs/(double)ncp->num_subfiles;
        color = (int)((double)myrank/ratio);
    }
    else
        color = myrank%ncp->num_subfiles;

#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): color=%d\n", __func__, myrank, color);
#endif
    /* key = myrank/comm_size; */

    /* TODO: fix error when using generated key value.
     * for now, just passing 0 value. */
    TRACE_COMM(MPI_Comm_split)(ncp->comm, color, myrank, &ncp->comm_sf);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_split");

    /* char path[1024], file[1024]; */
    /* find_path_and_fname(ncp->path, path, file); */
    /* sprintf(path_sf, "%s%d/%s", path, color, file); */
    sprintf(path_sf, "%s.subfile_%i.%s", ncp->path, color, "nc");

    void *ncp_sf;
    status = ncmpio_open(ncp->comm_sf, path_sf, ncp->iomode, ncp->ncid,
                         MPI_INFO_NULL, &ncp_sf);

    ncp->ncp_sf = (NC*) ncp_sf;
    return status;
}

/*----< ncmpio_subfile_close() >---------------------------------------------*/
int ncmpio_subfile_close(NC *ncp)
{
    int status = NC_NOERR;

    if (ncp->ncp_sf != NULL) {
        status = ncmpio_close(ncp->ncp_sf);
        if (status != NC_NOERR) return status;
        ncp->ncp_sf = NULL;
        MPI_Comm_free(&ncp->comm_sf);
    }

    /* reset values to 0 */
    is_partitioned = 0;

    return status;
}

/*----< ncmpio_subfile_partition() >-----------------------------------------*/
int ncmpio_subfile_partition(NC *ncp)
{
    int i, j, color, myrank, nprocs, status=NC_NOERR, num_subfiles;
    double ratio;

    MPI_Comm_rank(ncp->comm, &myrank);
    MPI_Comm_size(ncp->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank==0)  /* debug */
    {
        printf("rank(%d): is_partitioned=%d\n", myrank, is_partitioned);
    }
#endif
    if (is_partitioned == 1) return NC_NOERR;

    if (nprocs > ncp->num_subfiles) {
        ratio = (double)nprocs/(double)ncp->num_subfiles;
        color = (int)((double)myrank/ratio);
    }
    else
        color = myrank%ncp->num_subfiles;

#ifdef SUBFILE_DEBUG
    printf("%s: rank(%d): color=%d\n", __func__, myrank, color);
#endif

    /* check whether file is already partitioned */
    /* this is to handle when app has multiple ncmpi_enddef() */
    status = ncmpio_get_att(ncp, NC_GLOBAL, "_PnetCDF_SubFiling.num_subfiles",
                            &num_subfiles, MPI_INT);

    if (status == NC_ENOTATT) { /* if such attr doesn't exist */
        status = subfile_create(ncp);
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

        status = ncmpio_put_att(ncp, NC_GLOBAL, "_PnetCDF_SubFiling.num_subfiles",
                                NC_INT, 1, &ncp->num_subfiles, MPI_INT);
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    }
    else if (status == NC_NOERR) { /* attr is already set */
        assert(num_subfiles == ncp->num_subfiles);
    }
    else /* other error */
        return status;

    /* TODO: ignore UNLIMITED dims */
    /* NOTE: the following "for loop" should be before NC_begins() */

    /* adjust the hints to be used by PnetCDF; use the same value in master */
    ncp->ncp_sf->h_align = ncp->h_align;
    ncp->ncp_sf->v_align = ncp->v_align;
    ncp->ncp_sf->r_align = ncp->r_align;

    for(i=0; i<ncp->vars.ndefined; i++) { /* traverse all variables */
        NC_var **vpp = ncp->vars.value;
        NC_dim **dpp = ncp->dims.value;

        /* check attr for subfiles */
        int par_dim_id = (IS_RECVAR(vpp[i])?1:0); /* default is the most significant dim excluding NC_UNLIMITED */

        /* skip if the var is partitioned already */
        if (vpp[i]->dimids == NULL) continue;

        int par_dim_id_temp;
        status = ncmpio_get_att(ncp, i, "_PnetCDF_SubFiling.par_dim_id",
                                &par_dim_id_temp, MPI_INT);
        if (status == NC_NOERR) { /* if this attr exists */
#ifdef SUBFILE_DEBUG
            printf("par_dim_id_temp=%d\n", par_dim_id_temp);
            printf("par_dim_id=%d\n", par_dim_id);
#endif
            /* convert par_dim_id to dim index defined in a var */
            for (j=0; j<vpp[i]->ndims; j++)
                if (vpp[i]->dimids[j] == par_dim_id_temp) {
#ifdef SUBFILE_DEBUG
                    if (myrank == 0)
                        printf("dimids[%d]=%d\n", j, vpp[i]->dimids[j]);
#endif
                    par_dim_id = j;
                    break;
                }
        }
        else if (status != NC_ENOTATT) /* ignore if this attr does not exist */
            return status; /* other kind of error */

        if (par_dim_id < vpp[i]->ndims) {
            status = ncmpio_put_att(ncp, i, "_PnetCDF_SubFiling.par_dim_index",
                                    NC_INT, 1, &par_dim_id, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
        }

#ifdef SUBFILE_DEBUG
        if (myrank == 0)
            printf ("%s: var(%s): size of partitioning dim (id=%d)=%d\n",
                    __func__, vpp[i]->name, par_dim_id,
                    dpp[par_dim_id]->size);
#endif
        /* divide only when dim is partitionable */
        /* 1. skip sizeof(par_dim_id) is smaller than num_subfiles */
        /* 2. skip if ndims < SUBFILING_MIN_NDIMS */
        if (dpp[vpp[i]->dimids[par_dim_id]]->size/ncp->num_subfiles > 0 &&
            vpp[i]->ndims >= par_dim_id+1 &&
            vpp[i]->ndims >= SUBFILING_MIN_NDIMS) {
            int varid, j, jj, k;
            int var_ndims = vpp[i]->ndims; /* keep org ndims */
            int dimids[var_ndims];
            char *key[ncp->num_subfiles][var_ndims];

            for (jj=0; jj < ncp->num_subfiles; jj++)
                for (k=0; k<var_ndims; k++)
                    key[jj][k] = (char*) NCI_Calloc(NC_MAX_NAME, 1);

            /* save original value ndims to attribute before making to 0 */
            status = ncmpio_put_att(ncp, i, "_PnetCDF_SubFiling.ndims_org", NC_INT, 1, &vpp[i]->ndims, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

            int sf_range[ncp->num_subfiles][var_ndims][3];

            /* j: each dimension */
            /* subfile: create a var with partitioned dim sizes */
            for(j=0; j<var_ndims; j++) {
                MPI_Offset dim_sz0, dim_sz;
                char str[80], *dim_name;
                double x; /* = org dim size / num_subfiles */
                int min, max;

                dim_name = dpp[vpp[i]->dimids[j]]->name;
                dim_sz0 = dpp[vpp[i]->dimids[j]]->size; /* init both to org dim sz */
                dim_sz = dim_sz0;

                /* determine partition ratio */
                x = (double)(dim_sz0)/(double)(ncp->num_subfiles);

                /* don't partition dim if dim size is less than ratio x */
                if ((int)x < 1 && j == par_dim_id)
                    continue;

                /* don't need to check? */
                if (j == par_dim_id)
                {
                    double xx = x*(double)color;
                    double yy = x*(double)(color+1);
                    min = (int)xx+(color==0||(xx-(int)xx==0.0)?0:1);
                    max = (int)yy-(yy-(int)yy==0.0?1:0);
                    if ((MPI_Offset)max >= dim_sz0) max = dim_sz0-1;
                    dim_sz = max-min+1;
                }

                /* name is always NULL terminated */
                sprintf(str, "%s.%s", dim_name, vpp[i]->name);
#ifdef SUBFILE_DEBUG
                printf("rank(%d): new dim name = %s, dim_sz=%d\n", myrank, str, dim_sz);
#endif
                status = ncmpio_def_dim(ncp->ncp_sf, str, dim_sz, &dimids[j]);
                if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

                /* dpp_sf[color][j] = ncp->ncp_sf->dims.value[j]; */

                for (jj=0; jj < ncp->num_subfiles; jj++) {
                    double xx, yy;
                    sprintf(key[jj][j], "_PnetCDF_SubFiling.range(%s).subfile.%d", dim_name, jj); /* dim name*/
                    xx = x*(double)jj;
                    min = (int)xx+(jj==0||(xx-(int)xx==0.0)?0:1);
                    yy = x*(double)(jj+1);
                    max = (int)yy-(yy-(int)yy==0.0?1:0);
                    if ((MPI_Offset)max >= dim_sz0) max = dim_sz0-1;
#ifdef SUBFILE_DEBUG
                    if (myrank == 0) printf("subfile(%d): min=%d, max=%d\n", jj, min, max);
#endif
                    if (j == par_dim_id) { /* partitioning dims? */
                        sf_range[jj][j][0] = min;
                        sf_range[jj][j][1] = max;
                        sf_range[jj][j][2] = (max-min+1);
                    }
                    else {
                        sf_range[jj][j][0] = 0;
                        sf_range[jj][j][2] = dim_sz0;
                        sf_range[jj][j][1] = (dim_sz0!=0) ? dim_sz0-1 : 0;
                    }
                }
            } /* for each dim */

            /* master file: replace the original var with scalar var */
            vpp[i]->dimids_org = (int*) NCI_Malloc((size_t)vpp[i]->ndims * SIZEOF_INT);
            memcpy(vpp[i]->dimids_org, vpp[i]->dimids, (size_t)vpp[i]->ndims*SIZEOF_INT);
            vpp[i]->ndims_org = vpp[i]->ndims;
            vpp[i]->ndims = 0;
            if (IS_RECVAR(vpp[i])) {
                NCI_Free(vpp[i]->dimids);
                vpp[i]->dimids = NULL;
            }
            if (vpp[i]->shape != NULL) {
                NCI_Free(vpp[i]->shape);
                vpp[i]->shape = NULL;
            }
            vpp[i]->len = vpp[i]->xsz; /* size of type  */
            vpp[i]->num_subfiles = ncp->num_subfiles;

            status = ncmpio_put_att(ncp, i, "_PnetCDF_SubFiling.dimids_org",
                                    NC_INT, vpp[i]->ndims_org, vpp[i]->dimids_org, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

            status = ncmpio_put_att(ncp, i, "_PnetCDF_SubFiling.num_subfiles",
                                    NC_INT, 1, &ncp->num_subfiles, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

            for (jj=0; jj < ncp->num_subfiles; jj++)
                for (k=0; k < var_ndims; k++) {
                    status = ncmpio_put_att(ncp, i, key[jj][k], NC_INT,
                                            2, sf_range[jj][k], MPI_INT);
                    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
                }

            /* define a var with new dim */
            status = ncmpio_def_var(ncp->ncp_sf, (*vpp[i]).name,
                                    vpp[i]->xtype, var_ndims, dimids, &varid);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

            /* add an attribute about each dim's range in subfile */
            /* varid: var id in subfile */
            for (k=0; k<var_ndims; k++) {
                status = ncmpio_put_att(ncp->ncp_sf, varid, key[color][k],
                                        NC_INT, 2, sf_range[color][k], MPI_INT);
                if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
            }

            /* deallocate buffers */
            for (jj=0; jj < ncp->num_subfiles; jj++)
                for (k=0; k<var_ndims; k++)
                    NCI_Free(key[jj][k]);

            status = ncmpio_put_att(ncp->ncp_sf, varid,
                                    "_PnetCDF_SubFiling.subfile_index",
                                    NC_INT, 1, &color, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
        } /* end if() */
    } /* for each variable */

    is_partitioned = 1;

    return NC_NOERR;
}

/*----< ncmpio_subfile_getput_vars() >---------------------------------------*/
int
ncmpio_subfile_getput_vars(NC               *ncp,
                           NC_var           *varp,
                           const MPI_Offset  start[],
                           const MPI_Offset  count[],
                           const MPI_Offset  stride[],
                           void             *buf,
                           MPI_Offset        bufcount,
                           MPI_Datatype      buftype,
                           int               reqMode)
{
    int mpireturn, errs=0, status;
    NC_subfile_access *my_req, *others_req;
    int i, j, k, myrank, nprocs;
    int *count_my_req_per_proc, *count_others_req_per_proc;
    int varid, varid_sf;
    MPI_Offset *buf_count_my, *buf_count_others;
    void **xbuf=NULL, *cbuf=NULL;
    MPI_Request *requests=NULL;
    MPI_Status *statuses=NULL;
    int ndims_org = varp->ndims_org;
    int color;
    int nasyncios=0;

    MPI_Comm_rank(ncp->comm, &myrank);
    MPI_Comm_size(ncp->comm, &nprocs);

#ifdef SUBFILE_DEBUG
    for (i=0; i<ndims_org; i++)
        printf("rank(%d): %s: var(%s): start[%d]=%ld, count[%d]=%ld, stride[%d]=%ld, bufcount=%ld\n",
               myrank, __func__, varp->name, i, start[i], i, count[i], i,
               ((stride != NULL)?stride[i]:1), bufcount);
#endif

    /* check attr for subfiles */
    int par_dim_id = 0; /* default is the most significant dim */

    status = ncmpio_inq_varid(ncp, varp->name, &varid);
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

    status = ncmpio_get_att(ncp, varid, "_PnetCDF_SubFiling.par_dim_index",
                            &par_dim_id, MPI_INT);
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
        printf("_PnetCDF_SubFiling.par_dim_index of %s = %d\n", varp->name, par_dim_id);
#endif
    status = ncmpio_inq_varid(ncp->ncp_sf, varp->name, &varid_sf);
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

    status = ncmpio_get_att(ncp->ncp_sf, varid_sf, "_PnetCDF_SubFiling.subfile_index",
                            &color, MPI_INT);
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

    count_my_req_per_proc = (int *)NCI_Calloc((size_t)nprocs, SIZEOF_INT);
    count_others_req_per_proc = (int *)NCI_Calloc((size_t)nprocs, SIZEOF_INT);
    buf_count_my = (MPI_Offset *)NCI_Calloc((size_t)nprocs, SIZEOF_MPI_OFFSET);
    buf_count_others = (MPI_Offset *)NCI_Calloc((size_t)nprocs, SIZEOF_MPI_OFFSET);

    /* TODO: shouldn't it be initialized to 0? */
    for (i=0; i<nprocs; i++) {
        buf_count_my[i] = 1;
        buf_count_others[i] = 1;
    }

    /* allocate space for my_req */
    my_req = (NC_subfile_access*) NCI_Calloc((size_t)nprocs, sizeof(NC_subfile_access));

    /* init allocated my_req */
    for (i=0; i<nprocs; i++) {
        my_req[i].start = (MPI_Offset*)NCI_Malloc((size_t)(3 * ndims_org * SIZEOF_MPI_OFFSET));
        my_req[i].count = my_req[i].start + ndims_org;
        my_req[i].start_org = my_req[i].count + ndims_org;
        /* init start/count to -1 */
        for (j=0; j<ndims_org; j++) {
            my_req[i].start[j] = -1;
            my_req[i].count[j] = -1;
            my_req[i].start_org[j] = -1;
        }
    }

    /* i: for each subfile */
    for (i=0; i<ncp->num_subfiles; i++) {
        int flag = 0; /* set to 1 if par_dim_id is partitioned, initially 0 */
        double ratio = (double)nprocs/(double)ncp->num_subfiles;
        int aproc=-1; /* I/O delegate proc in subfile group */

        if (delegate_scheme == BALANCED) {
            double scaled, xx, yy;
            int min, max;

            xx = (ratio)*(double)i; /* i: each subfile */
            min = (int)xx+(i==0||(xx-(int)xx==0.0)?0:1);
            yy = (ratio)*(double)(i+1);
            max = (int)yy-(yy-(int)yy==0.0?1:0);
            if (max >= nprocs) max = nprocs-1;
            /* scaled = (double)random()/RAND_MAX; */
            scaled = (double)myrank/ratio-(double)color;
            aproc = (i==color)?myrank:((int)(min+(max-min+1)*scaled));
        }
        else if (delegate_scheme == ONE)
        {
            double xx;
            int min;
            xx = (ratio)*(double)i;
            min = (int)xx+(i==0||(xx-(int)xx==0.0)?0:1);
            aproc = (i==color)?myrank:min;
        }

        /* check out of range? */
        if (aproc >= nprocs)
            aproc = nprocs-1;

#ifdef SUBFILE_DEBUG
        printf("rank(%d): color=%d, subfile=%d, aproc=%d\n", myrank, color, i, aproc);
#endif

        /* j: for each dim starting from par_dim_id in round-robin manner */
        for (j=par_dim_id; j<par_dim_id+ndims_org; j++) {
            int jx = j%ndims_org;
            NC_dim *dimp = ncp->ncp_sf->dims.value[ncp->ncp_sf->vars.value[varid_sf]->dimids[jx]];
            int sf_range[2];
            int ii, jj, kk, stride_count;
            char key[256], *org_dim_name;

            org_dim_name = strtok(dimp->name, ".");
            sprintf(key, "_PnetCDF_SubFiling.range(%s).subfile.%d", org_dim_name, i); /* dim name*/
#ifdef SUBFILE_DEBUG
            if (myrank == 0) {
                printf("rank(%d): org_dim_name=%s\n", myrank, org_dim_name);
                printf("rank(%d): key=%s\n", myrank, key);
            }
#endif
            /* should get range info from the master file */
            status = ncmpio_get_att(ncp, varid, key, sf_range, MPI_INT);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

#ifdef SUBFILE_DEBUG
            /* if (myrank == 0) */
                printf("rank(%d): ndims_org=%d, rstart=%d, rend=%d, start[%d]=%d, count[%d]=%d\n",
                       myrank, ndims_org, sf_range[0], sf_range[1], jx, start[jx], jx, count[jx]);
#endif
            /* ii: traverse user's request range, incremented by stride count
               jj: traverse subfile range, incrementing sequentially
               kk: count belong to my subfile range */
            ii=start[jx], jj=sf_range[0], kk=0;
            stride_count = (stride == NULL) ? 1 : stride[jx];
            /* printf("stride_count[%d]=%d\n", jx, stride_count); */

            /* TODO: if sf_range is 1, count[] value is incorrect
               e.g., size of par_dim is 4 and the nproc is also 4 */
            int is_unlimited_dim = sf_range[1]-sf_range[0];

            /* skip the remaining dim if par_dim is not
               intersected with user's request */
            if (jx != par_dim_id && !flag) continue;

            while (ii < (start[jx]+(count[jx]*stride_count)) &&
                   jj <= (is_unlimited_dim==0?(count[jx]*stride_count-1):sf_range[1]))
            {
                if (ii < jj)
                    ii+=stride_count;
                else if (jj < ii)
                    jj++;
                else {
#ifdef SUBFILE_DEBUG
                    printf("rank(%d): var(%s): i=%d, j=%d, ii=%d, jj=%d, kk=%d, jx=%d\n", myrank, varp->name, i, j, ii, jj, kk, jx);
#endif
                    if (kk == 0) {
                        my_req[aproc].start[jx] = ii;
#ifdef SUBFILE_DEBUG
                        printf("rank(%d): var(%s): my_req[%d].start[%d]=%d\n",
                               myrank, varp->name, aproc, jx, ii);
#endif
                    }
                    if (jx == par_dim_id) flag = 1;
                    ii+=stride_count; jj++; kk++;
                }

            }
            if (kk > 0 && flag == 1) my_req[aproc].count[jx] = kk;
            /* adjust start offset based on subfile's range start.
               otherwise, there will be an out of bound error during I/O */
            if (my_req[aproc].start[jx] != -1) {
                my_req[aproc].start_org[jx] = my_req[aproc].start[jx];
                my_req[aproc].start[jx] -= sf_range[0];
            }
#ifdef SUBFILE_DEBUG
            /* if (myrank == 0) */
            {
                printf("rank(%d): my_req[%d].start[%d]=%d\n", myrank, aproc,
                       jx, my_req[aproc].start[jx]);
                printf("rank(%d): my_req[%d].count[%d]=%d\n", myrank, aproc,
                       jx, my_req[aproc].count[jx]);
            }
#endif
        } /* for each dim, j */
        if (my_req[aproc].start[0] == -1)
            count_my_req_per_proc[aproc] = 0;
        else
            count_my_req_per_proc[aproc] = 1;
    } /* for each subfile, i */

#ifdef SUBFILE_DEBUG
    for (i=0; i<ncp->num_subfiles; i++) {
        char str_st[100], str_st_org[100], str_ct[100], str_t1[10];
        sprintf(str_st, ">> rank(%d): subfile(%d): var(%s): start(", myrank, i, varp->name);
        sprintf(str_ct, ">> rank(%d): subfile(%d): count(", myrank, i);
        sprintf(str_st_org, "%d: start_org(", i);
        for (j=0; j<ndims_org; j++) {
            sprintf(str_t1, "%d", my_req[i].start[j]);
            strcat(str_st, str_t1);
            sprintf(str_t1, "%d", my_req[i].count[j]);
            strcat(str_ct, str_t1);
            sprintf(str_t1, "%d", my_req[i].start_org[j]);
            strcat(str_st_org, str_t1);
            if (j != (ndims_org-1)) {
                strcat(str_st, ",");
                strcat(str_ct, ",");
                strcat(str_st_org, ",");
            }

        }
        strcat(str_st, ")");
        strcat(str_ct, ")");
        strcat(str_st_org, ")");
        printf("%s\n", str_st);
        printf("%s\n", str_ct);
        printf("%s\n", str_st_org);
    }
#endif

    /* calculate buf_count based on subfile */
    /* TODO: should do only when count_my_req_per_proc[myrank] != 0??? */
    for (i=0; i<ndims_org; i++)
        buf_count_my[myrank] *= my_req[myrank].count[i];

    /* build other proc's request */
    others_req = (NC_subfile_access *) NCI_Calloc((size_t)nprocs, sizeof(NC_subfile_access));
    /* init allocated others_req */
    for (i=0; i<nprocs; i++) {
        others_req[i].start = (MPI_Offset *)NCI_Malloc((size_t)(3 * ndims_org * SIZEOF_MPI_OFFSET));
        others_req[i].count = others_req[i].start + ndims_org;
        others_req[i].start_org = others_req[i].count + ndims_org;
        /* init start/count to -1 */
        for (j=0; j<ndims_org; j++) {
            others_req[i].start[j] = -1;
            others_req[i].count[j] = -1;
            others_req[i].start_org[j] = -1;
        }
    }

#ifdef SUBFILE_DEBUG
    for (i=0; i<nprocs; i++) {
        printf("%d: count_my_req_per_proc[%d]=%d\n", myrank, i, count_my_req_per_proc[i]);
        printf("%d: count_others_req_per_proc[%d]=%d\n", myrank, i, count_others_req_per_proc[i]);
    }
#endif

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t51, "SSON --- getput_vars MPI_Alltoall", "", TAU_USER);
    TAU_PHASE_START(t51);
#endif

    TRACE_COMM(MPI_Alltoall)(count_my_req_per_proc, 1, MPI_INT, count_others_req_per_proc, 1, MPI_INT, ncp->comm);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t51);
#endif

#ifdef SUBFILE_DEBUG
    for (i=0; i<nprocs; i++) {
        printf("=> %d: count_others_req_per_proc[%d]=%d\n", myrank, i, count_others_req_per_proc[i]);
    }
#endif

    MPI_Offset buf_offset[nprocs];
    for (i=0; i<nprocs; i++)
        buf_offset[i] = 0;

#ifdef SUBFILE_DEBUG
    if (myrank == 0) printf("varname=%s: etype=0x%x, etype_size=%d\n", varp->name, varp->xtype, varp->xsz);
#endif

    /* find the ptype (primitive MPI data type) from buftype
     * el_size is the element size of ptype
     * bnelems is the total number of ptype elements in the I/O buffer, buf
     * fnelems is the number of nc variable elements in nc_type
     * nbytes is the amount of read/write in bytes
     */
    MPI_Datatype ptype;
    int el_size;
    int buftype_is_contig;
    MPI_Offset bnelems;

    if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, bufcount is ignored and will be recalculated to match
         * count[]. Note buf's data type must match the data type of
         * variable defined in the file - no data conversion will be done.
         */
        bufcount = 1;
        for (i=0; i<varp->ndims; i++) {
            if (count[i] < 0)  /* no negative count[] */
                DEBUG_RETURN_ERROR(NC_ENEGATIVECNT)
            bufcount *= count[i];
        }
        /* assign buftype match with the variable's data type */
        buftype = ncmpii_nc2mpitype(varp->xtype);
    }

    status = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems, NULL,
                                 &buftype_is_contig);
    /* bnelems now is the number of ptype in a buftype */
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)

#ifdef SUBFILE_DEBUG
    MPI_Aint lb, extent;
    MPI_Type_get_extent(buftype, &lb, &extent);
    printf("rank(%d): var(%s): ptype=0x%x, el_size=%d, bnelems=%d, buftype_is_contig=%d, lb=%d, extent=%d\n",
           myrank, varp->name, ptype, el_size, bnelems, buftype_is_contig, lb, extent);
#endif

    /* if buftype is non-contiguous, pack to contiguous buffer*/
    /* NOTE: no conversion and byte swap are performed here
       as they are done underneath layer */
    if (!buftype_is_contig && bufcount > 0 && bnelems > 0) {
        int position=0;
        MPI_Offset outsize = bnelems * bufcount * el_size;
        if (outsize  != (int)outsize) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        if (bufcount != (int)bufcount) DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        cbuf = NCI_Malloc((size_t)outsize);
        if (fIsSet(reqMode, NC_REQ_WR))
            MPI_Pack(buf, (int)bufcount, buftype, cbuf, (int)outsize,
                     &position, MPI_COMM_SELF);
    }
    else
        cbuf = (void *)buf;

    MPI_Offset diff[ndims_org];
    for (i=0; i < ndims_org; i++) {
        MPI_Offset stride_count;
        stride_count = (stride == NULL?1:stride[i]);
        /* making diff is necessary?? */
        diff[i] = ABS(my_req[myrank].start_org[i]-start[i])/stride_count;
#ifdef SUBFILE_DEBUG
        if (myrank == 0) printf("rank(%d): my_req[%d].start_org[%d]=%d, start[%d]=%d, diff[%d]=%lld\n", myrank,
               myrank, i, my_req[myrank].start_org[i], i, start[i], i, diff[i]);
#endif
    }

    for (i=0; i<ndims_org; i++) {
        int l;
        MPI_Offset tmp=1;
        for (l=i+1; l < ndims_org ; l++) {
            tmp *= my_req[myrank].count[l];
        }
        buf_offset[myrank] += tmp*diff[i];
#ifdef SUBFILE_DEBUG
        if (myrank == 0) printf("local: rank(%d): buf_offset[%d]=%d\n", myrank, myrank, buf_offset[myrank]);
#endif
    }

    buf_offset[myrank] *= el_size;

#ifdef SUBFILE_DEBUG
    printf("rank(%d): buf_offset[%d]=%d, buf_count_my[%d]=%d\n", myrank,
           myrank, buf_offset[myrank], myrank, buf_count_my[myrank]);

    printf ("rank(%d): %s: var(%s)\n", myrank, __func__, varp->name);
    for (i=0; i<ndims_org; i++)
        printf("my_req[%d].start[%d]=%ld, my_req[%d].count[%d]=%ld, my_req[%d].stride[%d]=%ld, bufcount=%ld\n",
               myrank, i, my_req[myrank].start[i], myrank, i,
               my_req[myrank].count[i], myrank, i,
               ((stride != NULL)?stride[i]:1), bufcount);
#endif

    int *array_of_requests;
    int *array_of_statuses;
    /* TODO: each proc can't get more than nprocs-1?? */
    array_of_requests = (int *)NCI_Malloc((size_t)nprocs * SIZEOF_INT);
    for (i=0; i<nprocs; i++)
        array_of_requests[i] = NC_REQ_NULL;

#ifdef SUBFILE_DEBUG
    printf("buf_count_my[%d]=%d\n", myrank, buf_count_my[myrank]);
#endif
    /* doing my portion of I/O */
    if (count_my_req_per_proc[myrank] != 0) {
        status = ncmpio_igetput_varm(ncp->ncp_sf,
                                     ncp->ncp_sf->vars.value[varid_sf],
                                     my_req[myrank].start,
                                     my_req[myrank].count,
                                     stride,
                                     NULL,
                                     (char*)cbuf+buf_offset[myrank],
                                     buf_count_my[myrank],
                                     (!buftype_is_contig?ptype:buftype),
                                     &array_of_requests[nasyncios++],
                                     reqMode,
                                     0);
        if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    }
#ifdef SUBFILE_DEBUG
    printf("rank(%d): var(%s): pushed local I/O to async calls...\n", myrank, varp->name);
#endif

    /* doing other proc's request to my subfile
       TODO: each proc can't get more than nprocs?? */
    requests = (MPI_Request *)NCI_Malloc((size_t)nprocs*2*sizeof(MPI_Request));

    j = 0;
#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t52, "SSON --- getput_vars: exch count_my_req", "", TAU_USER);
    TAU_PHASE_START(t52);
#endif
    /* TODO: need to exchange count_my_req_per_proc? */
    count_others_req_per_proc[myrank] = count_my_req_per_proc[myrank];

    for (i=0; i<nprocs; i++)
        if (count_my_req_per_proc[i] != 0 && i != myrank)
            TRACE_COMM(MPI_Irecv)(&count_others_req_per_proc[i], 1, MPI_INT, i, i+myrank, ncp->comm, &requests[j++]);

    for (i=0; i<nprocs; i++)
        if (count_others_req_per_proc[i] != 0 && i != myrank)
            TRACE_COMM(MPI_Isend)(&count_my_req_per_proc[i], 1, MPI_INT, i, i+myrank, ncp->comm, &requests[j++]);

    statuses = (MPI_Status *)NCI_Malloc((size_t)j*sizeof(MPI_Status));
    TRACE_COMM(MPI_Waitall)(j, requests, statuses);
#ifdef TAU_SSON
    TAU_PHASE_STOP(t52);
#endif
    check_err(MPI_Waitall);
    NCI_Free(statuses);
    NCI_Free(requests);

    j = 0;
    requests = (MPI_Request *)NCI_Malloc((size_t)nprocs*2*sizeof(MPI_Request));

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t53, "SSON --- getput_vars: exch start,count,", "", TAU_USER);
    TAU_PHASE_START(t53);
#endif
    /* when dest rank == myrank */
    if (count_others_req_per_proc[myrank] > 0)
        memcpy(others_req[myrank].start, my_req[myrank].start,
               (size_t)ndims_org * 3 * SIZEOF_MPI_OFFSET);

    /* exchange my_req's and others_req's start[], count[], and start_org[] */
    for (i=0; i<nprocs; i++) {
        if (count_others_req_per_proc[i] != 0 && i != myrank)
            /* MPI_Offset == MPI_LONG_LONG_INT */
            TRACE_COMM(MPI_Irecv)(others_req[i].start, 3*ndims_org, MPI_LONG_LONG_INT, i,
                      i+myrank, ncp->comm, &requests[j++]);
    }
    for (i=0; i<nprocs; i++) {
        if (count_my_req_per_proc[i] != 0 && i != myrank)
            TRACE_COMM(MPI_Isend)(my_req[i].start, 3 * ndims_org, MPI_LONG_LONG_INT, i,
                      i+myrank, ncp->comm, &requests[j++]);
    }

    statuses = (MPI_Status *)NCI_Malloc((size_t)j*sizeof(MPI_Status));
    TRACE_COMM(MPI_Waitall)(j, requests, statuses);
    NCI_Free(statuses);
    NCI_Free(requests);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t53);
#endif
    check_err(MPI_Waitall);

#ifdef SUBFILE_DEBUG
    /* DEBUG: print out others_req.{start,count} */
    for (i=0; i<nprocs; i++) {
        char str_st[100], str_st_org[100], str_ct[100], str_t1[10];
        sprintf(str_st, "%d: others.start(", i);
        sprintf(str_ct, "%d: others.count(", i);
        sprintf(str_st_org, "%d: others.start_org(", i);

        for (j=0; j<ndims_org; j++) {
            sprintf(str_t1, "%d", others_req[i].start[j]);
            strcat(str_st, str_t1);
            sprintf(str_t1, "%d", others_req[i].count[j]);
            strcat(str_ct, str_t1);
            sprintf(str_t1, "%d", others_req[i].start_org[j]);
            strcat(str_st_org, str_t1);
            if (j != (ndims_org-1)) {
                strcat(str_st, ",");
                strcat(str_ct, ",");
                strcat(str_st_org, ",");
            }

        }
        strcat(str_st, ")");
        strcat(str_ct, ")");
        strcat(str_st_org, ")");
        printf("%s\n", str_st);
        printf("%s\n", str_ct);
        printf("%s\n", str_st_org);
    }
#endif

    j = 0;
    requests = (MPI_Request *)NCI_Malloc((size_t)nprocs*2*sizeof(MPI_Request));

    xbuf = (void**)NCI_Calloc((size_t)nprocs, sizeof(void *));

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t54, "SSON --- getput_vars: exch buf", "", TAU_USER);
    TAU_PHASE_START(t54);
#endif

    for (i=0; i<nprocs; i++) {
        /* xbuf[i] = NULL; */
        buf_count_others[i] = 1;
        if (count_others_req_per_proc[i] != 0 && i != myrank) {
            for (k=0; k<ndims_org; k++) {
                buf_count_others[i] *= others_req[i].count[k];
            }
#ifdef SUBFILE_DEBUG
            printf("rank(%d): recv from rank %d: buf_count_others[%d]=%d\n", myrank, i, i, buf_count_others[i]);
#endif
            xbuf[i] = (void*)NCI_Calloc((size_t)buf_count_others[i], (size_t)el_size);
            TRACE_COMM(MPI_Irecv)(xbuf[i], buf_count_others[i], (!buftype_is_contig?ptype:buftype), i, i+myrank, ncp->comm, &requests[j++]);
        }
    }

    for (i=0; i<nprocs; i++) {
        buf_count_my[i] = 1;
        if (count_my_req_per_proc[i] != 0 && i != myrank) {
            MPI_Offset diff[ndims_org];
            for (k=0; k < ndims_org; k++) {
                int l;
                MPI_Offset stride_count, tmp=1;

                stride_count = (stride == NULL?1:stride[k]);
                diff[k] = ABS(my_req[i].start_org[k]-start[k])/stride_count;

                for (l=k+1; l < ndims_org; l++)
                    tmp *= my_req[i].count[l];

                buf_offset[i] += tmp*diff[k];
#ifdef SUBFILE_DEBUG
                if (myrank == 0) printf("rank(%d): subfile(%d): diff[%d]=%d\n", myrank, i, k, diff[k]);
#endif
            }

            buf_offset[i] *= el_size;

            for (k=0; k<ndims_org; k++) {
                buf_count_my[i] *= my_req[i].count[k];
            }
#ifdef SUBFILE_DEBUG
            printf("rank(%d): send to rank %d: buf_offset[%d]=%d, buf_count_my[%d]=%d\n", myrank, i, i, buf_offset[i], i, buf_count_my[i]);
#endif

            TRACE_COMM(MPI_Isend)((char*)cbuf+buf_offset[i], buf_count_my[i], (!buftype_is_contig?ptype:buftype), i, i+myrank, ncp->comm, &requests[j++]);
        } /* end if() */
    } /* end for() */

    statuses = (MPI_Status *)NCI_Malloc((size_t)j * sizeof(MPI_Status));
    TRACE_COMM(MPI_Waitall)(j, requests, statuses);
    NCI_Free(statuses);
    NCI_Free(requests);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t54);
#endif
    check_err(MPI_Waitall);

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t55, "SSON --- getput_vars igetput", "", TAU_USER);
    TAU_PHASE_START(t55);
#endif

    for (i=0; i<nprocs; i++) {
        if (count_others_req_per_proc[i] != 0 && i != myrank) {
            status = ncmpio_igetput_varm(ncp->ncp_sf,
                                         ncp->ncp_sf->vars.value[varid_sf],
                                         others_req[i].start,
                                         others_req[i].count,
                                         stride,
                                         NULL,
                                         xbuf[i],
                                         buf_count_others[i],
                                         (!buftype_is_contig?ptype:buftype),
                                         &array_of_requests[nasyncios++],
                                         reqMode,
                                         0);
            if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
        }
    }
#ifdef TAU_SSON
    TAU_PHASE_STOP(t55);
#endif

#ifdef SUBFILE_DEBUG
    printf("rank(%d): nasyncios=%d\n", myrank, nasyncios);
#endif

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t56, "SSON --- getput_vars ncmpi_wait_all", "", TAU_USER);
    TAU_PHASE_START(t56);
#endif
    /*
    double stime, wait_time;
    stime = MPI_Wtime();
    */
    array_of_statuses = (int *)NCI_Malloc((size_t)nasyncios * SIZEOF_INT);
    status = ncmpio_wait(ncp->ncp_sf, nasyncios, array_of_requests, array_of_statuses, NC_REQ_COLL);
    if (status != NC_NOERR) DEBUG_RETURN_ERROR(status)
    NCI_Free(array_of_statuses);
    NCI_Free(array_of_requests);

    /*
    int pending_nreqs;
    status = ncmpio_inq_nreqs(ncp, &pending_nreqs);
    printf("myrank(%d): pending_nreqs=%d\n", myrank, pending_nreqs);
    wait_time = MPI_Wtime() - stime;
    printf("myrank(%d): ncmpio_wait time = %f\n", myrank, wait_time);
    */
#ifdef TAU_SSON
    TAU_PHASE_STOP(t56);
#endif

#ifdef SUBFILE_DEBUG
    printf("rank(%d): var(%s): after ncmpi_wait_all()\n", myrank, varp->name);
#endif

    /* MPI_Barrier(ncp->comm); */

    /* free all allocated memories */

    for(i=0; i<nprocs; i++)
        if (xbuf[i] != NULL) NCI_Free(xbuf[i]);
    NCI_Free(xbuf);

    if (cbuf != NULL && cbuf != buf) {
        NCI_Free(cbuf);
    }

    NCI_Free(count_my_req_per_proc);
    NCI_Free(count_others_req_per_proc);
    NCI_Free(buf_count_my);
    NCI_Free(buf_count_others);

    for (i=0; i<nprocs; i++) {
        NCI_Free(my_req[i].start);
        NCI_Free(others_req[i].start);
    }
    NCI_Free(my_req);
    NCI_Free(others_req);

    return status;
}
