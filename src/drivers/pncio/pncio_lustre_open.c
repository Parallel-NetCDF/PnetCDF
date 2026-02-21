/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/errno.h>

#include <fcntl.h>      /* open(), O_CREAT */
#include <sys/types.h>  /* open() */
#include <libgen.h>     /* dirname() */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#ifndef PATH_MAX
#define PATH_MAX 65535
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h> /* open(), fstat() */
#endif

#include <mpi.h>

#include "pncio.h"

#ifdef MIMIC_LUSTRE
#define xstr(s) str(s)
#define str(s) #s
#define STRIPE_SIZE 64
#define STRIPE_COUNT 4
#endif

#ifdef HAVE_LUSTRE
/* /usr/include/lustre/lustreapi.h
 * /usr/include/linux/lustre/lustre_user.h
 */
#include <lustre/lustreapi.h>

#define PNETCDF_LUSTRE_DEBUG
// #define PNETCDF_LUSTRE_DEBUG_VERBOSE

#define PATTERN_STR(pattern, int_str) ( \
    (pattern == LLAPI_LAYOUT_DEFAULT)      ? "LLAPI_LAYOUT_DEFAULT" : \
    (pattern == LLAPI_LAYOUT_RAID0)        ? "LLAPI_LAYOUT_RAID0" : \
    (pattern == LLAPI_LAYOUT_WIDE)         ? "LLAPI_LAYOUT_WIDE" : \
    (pattern == LLAPI_LAYOUT_MDT)          ? "LLAPI_LAYOUT_MDT" : \
    (pattern == LLAPI_LAYOUT_OVERSTRIPING) ? "LLAPI_LAYOUT_OVERSTRIPING" : \
    (pattern == LLAPI_LAYOUT_SPECIFIC)     ? "LLAPI_LAYOUT_SPECIFIC" : \
    int_str)

#define PRINT_LAYOUT(val) { \
    char int_str[32]; \
    snprintf(int_str, 32, "%lu", val); \
    printf("\t%-14s = %-25s (0x%lx)\n",#val,PATTERN_STR(val, int_str),val); \
}

#ifdef HAVE_LLAPI_GET_OBD_COUNT

/*----< get_total_avail_osts() >---------------------------------------------*/
static
int get_total_avail_osts(const char *path)
{
    char *dir_path=NULL, *path_copy=NULL;
    int err, ost_count=0, is_mdt=0;
    struct stat sb;

    path_copy = NCI_Strdup(path);

    err = stat(path_copy, &sb);
    if (errno == ENOENT) { /* file does not exist, try folder */
        /* get the parent folder name */
        dir_path = dirname(path_copy);
        err = stat(dir_path, &sb);
    }
    if (err != 0) {
        printf("Warning at %s (%d): path \"%s\" stat() failed (%s)\n",
               __func__,__LINE__,path,strerror(errno));
        goto err_out;
    }

    /* llapi_get_obd_count() only works for directories */
    if (S_ISDIR(sb.st_mode))
        dir_path = (dir_path == NULL) ? path_copy : dir_path;
    else
        /* get the parent folder name */
        dir_path = dirname(path_copy);

    err = llapi_get_obd_count(dir_path, &ost_count, is_mdt);
    if (err != 0) {
        printf("Warning at %d: path \"%s\" llapi_get_obd_count() failed (%s)\n",
               __LINE__,dir_path,strerror(errno));
        ost_count = 0;
    }

err_out:
    if (path_copy != NULL) NCI_Free(path_copy);

    return ost_count;
}

#else

/*----< get_total_avail_osts() >---------------------------------------------*/
static
int get_total_avail_osts(const char *filename)
{
    char *dirc=NULL, *dname, *tail, **members=NULL, *buffer=NULL;
    char pool_name[64], fsname[64], full_pool_name[128];
    int err, dd, num_members=0;
    int max_members = 2048;    /* Maximum number of members to retrieve */
    int buffer_size = 1048576; /* Buffer size for member names */
    struct llapi_layout *layout=NULL;

    dirc = NCI_Strdup(filename);

    struct stat sb;
    if (stat(filename, &sb) == 0 && S_ISDIR(sb.st_mode))
        dname = dirc;
    else
        /* find the parent folder name */
        dname = dirname(dirc);

    dd = open(dname, O_RDONLY, 0600);
    if (dd < 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) fails to open folder %s (%s)\n",
               __FILE__,__LINE__, dname, strerror(errno));
#endif
        goto err_out;
    }

    /* obtain Lustre layout object */
    layout = llapi_layout_get_by_fd(dd, LLAPI_LAYOUT_GET_COPY);
    if (layout == NULL) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_get_by_fd() fails (%s)\n",
               __FILE__, __LINE__,strerror(errno));
#endif
        goto err_out;
    }

    /* find the pool name */
    err = llapi_layout_pool_name_get(layout, pool_name, sizeof(pool_name)-1);
    if (err < 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_pool_name_get() fails (%s)\n",
               __FILE__, __LINE__,strerror(errno));
#endif
        goto err_out;
    }
    else if (pool_name[0] == '\0') {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("%s at %d: %s has NO Pool Name\n",__FILE__, __LINE__,dname);
#endif
        goto err_out;
    }
    /* For example, Perlmutter @NERSC, pool_name "original" is returned */

    /* Using pool_name returned from llapi_layout_pool_name_get() is not enough
     * when calling  llapi_get_poolmembers(). We need to prepend it with
     * 'fsname', which can be obtained by calling llapi_getname(). Note that
     * console command 'lfs getname -n' returns fsname. For example, on
     * Perlmutter @NERSC:
     *    login39::~/Lustre(12:52) #1165  lfs getname -n $SCRATCH/dummy
     *    scratch
     */
    err = llapi_getname(dname, fsname, 63);
    if (err < 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_getname() fails (%s)\n",
               __FILE__, __LINE__,strerror(errno));
#endif
        goto err_out;
    }

    /* When dname is a folder, fsname returned from llapi_getname() may contain
     * a trailing ID, e.g.  scratch-ffff9ca88d9bd800. Must remove the trailing
     * ID, otherwise llapi_get_poolmembers() is not able to find it.
     */
    tail = strchr(fsname, '-');
    if (tail != NULL) *tail = '\0';

    /* In case either pool_name and fsname are empty. For example, on Polaris
     * @ALCF, the returned pool_name is empty, but fsname is not.
     */
    if (pool_name[0] == '\0' && fsname[0] == '\0')
        goto err_out;
    else if (pool_name[0] == '\0')
        strcpy(full_pool_name, fsname);
    else if (fsname[0] == '\0')
        strcpy(full_pool_name, pool_name);
    else
        sprintf(full_pool_name, "%s.%s", fsname, pool_name);

#ifdef PNETCDF_LUSTRE_DEBUG_VERBOSE
    printf("%s at %d: file=%s dir=%s pool=%s fsname=%s full_pool_name=%s\n",
           __func__,__LINE__, filename,dname,pool_name,fsname,full_pool_name);
#endif

    /* Allocate memory for the members and buffer */
    members = (char **)NCI_Malloc(max_members * sizeof(char *));
    buffer = (char *)NCI_Malloc(buffer_size);

    /* obtain pool's info */
    num_members = llapi_get_poolmembers(full_pool_name, members, max_members,
                                        buffer, buffer_size);
#ifdef PNETCDF_LUSTRE_DEBUG_VERBOSE
    if (num_members > 0) {
        int i, min_nmembers = MIN(num_members, 10);
        printf("%s at %d: Found %d members for pool '%s':\n",
               __func__,__LINE__,num_members, pool_name);
        printf("\tFirst %d OSTs and last are\n",min_nmembers);
        for (i=0; i<min_nmembers; i++)
            printf("\t\tmember[%3d] %s\n",i,members[i]);
        printf("\t ...\tmember[%3d] %s\n",num_members-1,members[num_members-1]);
        printf("------------------------------------\n\n");
    } else {
        printf("%s at %d: EOVERFLOW=%d EINVAL=%d\n",__func__,__LINE__,EOVERFLOW,EINVAL);
        printf("%s at %d: No members found for pool '%s' or an error occurred num_members=%d (%s).\n",
               __func__,__LINE__,pool_name, num_members, strerror(errno));
    }
#endif

err_out:
    if (dd >= 0) close(dd);
    if (layout != NULL) llapi_layout_free(layout);
    if (dirc != NULL) NCI_Free(dirc);
    if (buffer != NULL) NCI_Free(buffer);
    if (members != NULL) NCI_Free(members);

    return num_members;
}
#endif

static
int compare(const void *a, const void *b)
{
     if (*(uint64_t*)a > *(uint64_t*)b) return (1);
     if (*(uint64_t*)a < *(uint64_t*)b) return (-1);
     return (0);
}

static
int sort_ost_ids(struct llapi_layout *layout,
                 uint64_t             stripe_count,
                 uint64_t            *osts)
{
    uint64_t i, numOSTs;

    for (i=0; i<stripe_count; i++) {
        if (llapi_layout_ost_index_get(layout, i, &osts[i]) != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
            printf("Error at %s (%d) llapi_layout_ost_index_get(%lu) (%s)\n",
                   __FILE__,__LINE__,i,strerror(errno));
#endif
            return stripe_count;
        }
    }

    /* count the number of unique OST IDs. When Lustre overstriping is
     * used, the unique OSTs may be less than stripe_count.
     */
    qsort(osts, stripe_count, sizeof(uint64_t), compare);
    numOSTs = 0;
    for (i=1; i<stripe_count; i++)
        if (osts[i] > osts[numOSTs])
            osts[++numOSTs] = osts[i];

    return (numOSTs + 1);
}

/*----< get_striping() >-----------------------------------------------------*/
static
uint64_t get_striping(int         fd,
                      const char *path,
                      uint64_t   *pattern,
                      uint64_t   *stripe_count,
                      uint64_t   *stripe_size,
                      uint64_t   *start_iodevice)
{
    int err;
    struct llapi_layout *layout;
    uint64_t *osts=NULL, numOSTs=0;
#ifdef PNETCDF_LUSTRE_DEBUG
    char int_str[32];
#endif

    *pattern = LLAPI_LAYOUT_RAID0;
    *stripe_count = LLAPI_LAYOUT_DEFAULT;
    *stripe_size = LLAPI_LAYOUT_DEFAULT;
    *start_iodevice = LLAPI_LAYOUT_DEFAULT;

    layout = llapi_layout_get_by_fd(fd, LLAPI_LAYOUT_GET_COPY);
    if (layout == NULL) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_get_by_fd() fails\n",
                __FILE__, __LINE__);
#endif
        goto err_out;
    }

    err = llapi_layout_pattern_get(layout, pattern);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        snprintf(int_str, 32, "%lu", *pattern);
        printf("Error at %s (%d) llapi_layout_pattern_get() fails to get pattern %s\n",
                __FILE__, __LINE__, PATTERN_STR(*pattern, int_str));
#endif
        goto err_out;
    }

    /* obtain file striping count */
    err = llapi_layout_stripe_count_get(layout, stripe_count);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        snprintf(int_str, 32, "%lu", *stripe_count);
        printf("Error at %s (%d) llapi_layout_stripe_count_get() fails to get stripe count %s\n",
            __FILE__, __LINE__, PATTERN_STR(*stripe_count, int_str));
#endif
        goto err_out;
    }

    /* obtain file striping unit size */
    err = llapi_layout_stripe_size_get(layout, stripe_size);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        snprintf(int_str, 32, "%lu", *stripe_size);
        printf("Error at %s (%d) llapi_layout_stripe_size_get() fails to get stripe size %s\n",
            __FILE__,__LINE__, PATTERN_STR(*stripe_size, int_str));
#endif
        goto err_out;
    }

    /* /usr/include/linux/lustre/lustre_user.h
     * The stripe size fields are shared for the extension size storage,
     * however the extension size is stored in KB, not bytes.
     *     #define SEL_UNIT_SIZE 1024llu
     * Therefore, the default stripe_size is (SEL_UNIT_SIZE * 1024)
     */

    if (*stripe_count == LLAPI_LAYOUT_DEFAULT ||  /* not set */
        *stripe_count == LLAPI_LAYOUT_INVALID ||  /* invalid */
        *stripe_count == LLAPI_LAYOUT_WIDE    ||  /* all system's OSTs */
        *stripe_count > 1048576) {                /* abnormally large number */
        return 0;
    }

    /* obtain all OST IDs */
    osts = (uint64_t*) NCI_Malloc(sizeof(uint64_t) * (*stripe_count));
    if (llapi_layout_ost_index_get(layout, 0, &osts[0]) != 0) {
        /* check if is a folder */
        struct stat path_stat;
        fstat(fd, &path_stat);
#ifdef PNETCDF_LUSTRE_DEBUG_VERBOSE
        if (S_ISREG(path_stat.st_mode)) /* not a regular file */
            printf("%s at %d: %s is a regular file\n",__func__,__LINE__,path);
        else if (S_ISDIR(path_stat.st_mode))
            printf("%s at %d: %s is a folder\n",__func__,__LINE__,path);
        else
#endif
        if (!S_ISREG(path_stat.st_mode) && /* not a regular file */
            !S_ISDIR(path_stat.st_mode)) { /* not a folder */
#ifdef PNETCDF_LUSTRE_DEBUG
            printf("Error at %s (%d) calling fstat() file %s (neither a regular file nor a folder)\n", \
                    __FILE__, __LINE__, path);
#endif
            goto err_out;
        }

        *start_iodevice = LLAPI_LAYOUT_DEFAULT;
        numOSTs = *stripe_count;

        goto err_out;
    }
    *start_iodevice = osts[0];

    numOSTs = sort_ost_ids(layout, *stripe_count, osts);
    assert(numOSTs <= *stripe_count);

err_out:
    if (osts != NULL) NCI_Free(osts);
    if (layout != NULL) llapi_layout_free(layout);

    return numOSTs;
}

/*----< set_striping() >-----------------------------------------------------*/
static
int set_striping(const char *path,
                 uint64_t    pattern,
                 uint64_t    numOSTs,
                 uint64_t    stripe_count,
                 uint64_t    stripe_size,
                 uint64_t    start_iodevice)
{
    int fd=-1, err=0;

    struct llapi_layout *layout = llapi_layout_alloc();
    if (layout == NULL) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_alloc() fails (%s)\n",
               __FILE__, __LINE__, strerror(errno));
#endif
        goto err_out;
    }

    /* When an abnormally large stripe_count is set by users, Lustre may just
     * allocate the total number of available OSTs, instead of returning an
     * error.
     */
    err = llapi_layout_stripe_count_set(layout, stripe_count);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_stripe_count_set() fails set stripe count %lu (%s)\n",
               __FILE__, __LINE__, stripe_count, strerror(errno));
#endif
        goto err_out;
    }

    err = llapi_layout_stripe_size_set(layout, stripe_size);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_stripe_size_set() fails to set stripe size %lu (%s)\n",
               __FILE__, __LINE__, stripe_size, strerror(errno));
#endif
        goto err_out;
    }

    if (pattern == LLAPI_LAYOUT_OVERSTRIPING) {
        uint64_t i, ost_id;
        if (start_iodevice == LLAPI_LAYOUT_DEFAULT)
            start_iodevice = 0;
        for (i=0; i<stripe_count; i++) {
            ost_id = (start_iodevice + i) % numOSTs;
            err = llapi_layout_ost_index_set(layout, i, ost_id);
            if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
                printf("Error at %s (%d) llapi_layout_ost_index_set() fails to set OST index %lu to %lu (%s)\n",
                       __FILE__, __LINE__, i, ost_id, strerror(errno));
#endif
                goto err_out;
            }
        }
    }
    else if (start_iodevice != LLAPI_LAYOUT_DEFAULT) {
        /* When an abnormally large start_iodevice is set by users, Lustre may
         * return an error. Instead fail will occur later at calling
         * llapi_layout_file_create().
         */
        err = llapi_layout_ost_index_set(layout, 0, start_iodevice);
        if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
            printf("Error at %s (%d) llapi_layout_ost_index_set() fails to set start iodevice %lu (%s)\n",
                   __FILE__, __LINE__, start_iodevice, strerror(errno));
#endif
            goto err_out;
        }
    }

    err = llapi_layout_pattern_set(layout, pattern);
    if (err != 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        char int_str[32];
        snprintf(int_str, 32, "%lu", pattern);
        printf("Error at %s (%d) llapi_layout_pattern_set() fails to set pattern %s (%s)\n",
               __FILE__, __LINE__, PATTERN_STR(pattern, int_str), strerror(errno));
#endif
        goto err_out;
    }

    /* create a new file with desired striping */
    fd = llapi_layout_file_create(path, O_CREAT|O_RDWR, PNCIO_PERM, layout);
    if (fd < 0) {
#ifdef PNETCDF_LUSTRE_DEBUG
        printf("Error at %s (%d) llapi_layout_file_create() fails (%s)\n",
               __FILE__, __LINE__, strerror(errno));
#endif
        goto err_out;
    }

err_out:
    if (layout != NULL) llapi_layout_free(layout);

#ifdef PNETCDF_LUSTRE_DEBUG
    if (fd < 0)
        printf("Error at %s (%d) fails to create file %s with desired file striping. PnetCDF now tries to inherit it from the parent folder.\n",
                __FILE__,__LINE__, path);
#endif

    return fd;
}
#endif

/*----< Lustre_set_cb_node_list() >------------------------------------------*/
/* Construct the list of I/O aggregators. It sets the followings.
 *   fd->hints->cb_nodes and set file info for hint cb_nodes.
 *   fd->hints->ranklist[], an int array of size fd->hints->cb_nodes.
 *   fd->is_agg: indicating whether this rank is an I/O aggregator
 *   fd->my_cb_nodes_index: index into fd->hints->ranklist[]. -1 if N/A
 */
static
int Lustre_set_cb_node_list(PNCIO_File *fd)
{
    int i, j, k, rank, nprocs, num_aggr, striping_factor;
    int *nprocs_per_node, **ranks_per_node;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    /* number of MPI processes running on each node */
    nprocs_per_node = (int *) NCI_Calloc(fd->node_ids.num_nodes, sizeof(int));

    for (i=0; i<nprocs; i++) nprocs_per_node[fd->node_ids.ids[i]]++;

    /* construct rank IDs of MPI processes running on each node */
    ranks_per_node = (int **) NCI_Malloc(sizeof(int*) * fd->node_ids.num_nodes);
    ranks_per_node[0] = (int *) NCI_Malloc(sizeof(int) * nprocs);
    for (i=1; i<fd->node_ids.num_nodes; i++)
        ranks_per_node[i] = ranks_per_node[i - 1] + nprocs_per_node[i - 1];

    for (i=0; i<fd->node_ids.num_nodes; i++) nprocs_per_node[i] = 0;

    /* Populate ranks_per_node[], list of MPI ranks running on each node.
     * Populate nprocs_per_node[], number of MPI processes on each node.
     */
    for (i=0; i<nprocs; i++) {
        k = fd->node_ids.ids[i];
        ranks_per_node[k][nprocs_per_node[k]] = i;
        nprocs_per_node[k]++;
    }

    /* To save a call to MPI_Bcast(), all processes run the same codes below to
     * calculate num_aggr, the number of aggregators (later becomes cb_nodes).
     *
     * The calculation is based on the number of compute nodes,
     * fd->node_ids.num_nodes, and processes per node, nprocs_per_node. At this
     * moment, all processes should have obtained the Lustre file striping
     * settings.
     */
    striping_factor = fd->hints->striping_factor;

    if (striping_factor > nprocs) {
        /* When number of MPI processes is less than striping_factor, set
         * num_aggr to the max number less than nprocs that divides
         * striping_factor. An naive way is:
         *     num_aggr = nprocs;
         *     while (striping_factor % num_aggr > 0)
         *         num_aggr--;
         * Below is equivalent, but faster.
         */
        int divisor = 2;
        num_aggr = 1;
        /* try to divide */
        while (striping_factor >= divisor * divisor) {
            if ((striping_factor % divisor) == 0) {
                if (striping_factor / divisor <= nprocs) {
                    /* The value is found ! */
                    num_aggr = striping_factor / divisor;
                    break;
                }
                /* if divisor is less than nprocs, divisor is a solution,
                 * but it is not sure that it is the best one
                 */
                else if (divisor <= nprocs)
                    num_aggr = divisor;
            }
            divisor++;
        }
    }
    else { /* striping_factor <= nprocs */
        /* Select striping_factor processes to be I/O aggregators. Note this
         * also applies to collective reads to allow more/less aggregators. In
         * most cases, more aggregators yields better read performance.
         */
        if (fd->hints->cb_nodes == 0) {
            /* User did not set hint "cb_nodes" */
            if (nprocs >= striping_factor * 8 &&
                nprocs/fd->node_ids.num_nodes >= 8)
                num_aggr = striping_factor * 8;
            else if (nprocs >= striping_factor * 4 &&
                     nprocs/fd->node_ids.num_nodes >= 4)
                num_aggr = striping_factor * 4;
            else if (nprocs >= striping_factor * 2 &&
                     nprocs/fd->node_ids.num_nodes >= 2)
                num_aggr = striping_factor * 2;
            else
                num_aggr = striping_factor;
        }
        else if (fd->hints->cb_nodes <= striping_factor) {
            /* User has set hint cb_nodes and cb_nodes <= striping_factor.
             * Ignore user's hint and try to set cb_nodes to be at least
             * striping_factor.
             */
            num_aggr = striping_factor;
        }
        else {
            /* User has set hint cb_nodes and cb_nodes > striping_factor */
            if (nprocs < fd->hints->cb_nodes)
                num_aggr = nprocs; /* BAD cb_nodes set by users */
            else
                num_aggr = fd->hints->cb_nodes;
        }

        /* Number of processes per node may not be enough to be picked as
         * aggregators. If this case, reduce num_aggr (cb_nodes). Consider the
         * following case:
         *   number of nodes = 7,
         *   number of processes = 18,
         *   striping_factor = 8,
         *   cb_nodes = 16.
         * Nodes in this case, nodes 0, 1, 2, 3 run 3 processes each and nodes
         * 4, 5, 6 run 2 processes each. In order to keep each OST only
         * accessed by one or more aggregators running on the same compute
         * node, cb_nodes should be reduced to 8. Thus the ranks of aggregators
         * become 0, 3, 6, 9, 12, 14, 16, 1. The aggregator-OST mapping
         * becomes below.
         *   Aggregator  0, running on node 0, access OST 0.
         *   Aggregator  3, running on node 1, access OST 1.
         *   Aggregator  6, running on node 2, access OST 2.
         *   Aggregator  9, running on node 3, access OST 3.
         *   Aggregator 12, running on node 4, access OST 4.
         *   Aggregator 14, running on node 5, access OST 5.
         *   Aggregator 16, running on node 6, access OST 6.
         *   Aggregator  1, running on node 0, access OST 7.
         *
         * Another case (the total number of processes changes to 25):
         *   number of nodes = 7,
         *   number of processes = 25,
         *   striping_factor = 8,
         *   cb_nodes = 16.
         * In this case, nodes 0, 1, 2, 3 run 4 processes each and nodes 4, 5,
         * 6 run 3 processes each. cb_nodes should remain 16 and the ranks of
         * aggregators become 0, 4, 8, 12, 16, 19, 22, 1, 2, 6, 10, 14, 18, 21,
         * 24, 3. The aggregator-OST mapping becomes below.
         *   Aggregators  0,  2, running on node 0, access OST 0.
         *   Aggregators  4,  6, running on node 1, access OST 1.
         *   Aggregators  8, 10, running on node 2, access OST 2.
         *   Aggregators 12, 14, running on node 3, access OST 3.
         *   Aggregators 16, 18, running on node 4, access OST 4.
         *   Aggregators 19, 21, running on node 5, access OST 5.
         *   Aggregators 22, 24, running on node 6, access OST 6.
         *   Aggregator   3,     running on node 0, access OST 7.
         */
        int max_nprocs_node = 0;
        for (i=0; i<fd->node_ids.num_nodes; i++)
            max_nprocs_node = MAX(max_nprocs_node, nprocs_per_node[i]);
        int max_naggr_node = striping_factor / fd->node_ids.num_nodes;
        if (striping_factor % fd->node_ids.num_nodes) max_naggr_node++;
        /* max_naggr_node is the max number of processes per node to be picked
         * as aggregator in each round.
         */
        int rounds = num_aggr / striping_factor;
        if (num_aggr % striping_factor) rounds++;
        while (max_naggr_node * rounds > max_nprocs_node) rounds--;
        num_aggr = striping_factor * rounds;
    }

    /* TODO: the above setting for num_aggr is for collective writes. Should
     * collective reads use the same?  Or just set cb_nodes to the number of
     * nodes.
     */

    /* Next step is to determine the MPI rank IDs of I/O aggregators and add
     * them into ranklist[]. Note fd->hints->ranklist will be freed in
     * PNCIO_File_close().
     */
    fd->hints->ranklist = (int *) NCI_Malloc(num_aggr * sizeof(int));
    if (fd->hints->ranklist == NULL)
        return NC_ENOMEM;

    int block_assignment=0;
#ifdef TRY_AGGR_BLOCK_ASSIGNMENT
    {
        char *env_str;
        if ((env_str = getenv("PNETCDF_USE_BLOCK_ASSIGN")) != NULL)
            block_assignment = (strcasecmp(env_str, "true") == 0) ? 1 : 0;
        if (rank == 0)
            printf("%s %d: PNETCDF_USE_BLOCK_ASSIGN = %d\n",
            __func__,__LINE__,block_assignment);
    }
#endif

    if (striping_factor <= fd->node_ids.num_nodes) {
        /* When number of OSTs is less than number of compute nodes, first
         * select number of nodes equal to the number of OSTs by spread the
         * selection evenly across all compute nodes (i.e. with a stride
         * between every 2 consecutive nodes).
         * Selection of MPI ranks can be done in 2 ways.
         * 1. block assignment
         *    Select ranks from a node and then move on to the next node.
         * 2. cyclic assignment
         *    Select ranks round-robin across all selected nodes.
         * Note when selecting ranks within a node, the ranks are evenly spread
         * among all processes in the node.
         */
        if (block_assignment) {
            int n=0;
            int remain = num_aggr % striping_factor;
            int node_stride = fd->node_ids.num_nodes / striping_factor;
            /* walk through each node and pick aggregators */
            for (j=0; j<fd->node_ids.num_nodes; j+=node_stride) {
                /* Selecting node IDs with a stride. j is the node ID */
                int nranks_per_node = num_aggr / striping_factor;
                /* front nodes may have 1 more to pick */
                if (remain > 0 && j/node_stride < remain) nranks_per_node++;
                int rank_stride = nprocs_per_node[j] / nranks_per_node;
                for (k=0; k<nranks_per_node; k++) {
                    /* Selecting rank IDs within node j with a stride */
                    fd->hints->ranklist[n] = ranks_per_node[j][k*rank_stride];
                    if (++n == num_aggr) {
                        j = fd->node_ids.num_nodes; /* break loop j */
                        break; /* loop k */
                    }
                }
            }
        }
        else {
            int avg = num_aggr / striping_factor;
            int stride = fd->node_ids.num_nodes / striping_factor;
            if (num_aggr % striping_factor) avg++;
            for (i = 0; i < num_aggr; i++) {
                /* j is the selected node ID. This selection is round-robin
                 * across selected nodes.
                 */
                j = (i % striping_factor) * stride;
                k = (i / striping_factor) * (nprocs_per_node[j] / avg);
                assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
        }
    }
    else { /* striping_factor > fd->node_ids.num_nodes */
        /* When number of OSTs is more than number of compute nodes, I/O
         * aggregators are selected from all nodes. Within each node,
         * aggregators are spread evenly instead of the first few ranks.
         */
        int *naggr_per_node, *idx_per_node, avg;
        idx_per_node = (int*) NCI_Calloc(fd->node_ids.num_nodes, sizeof(int));
        naggr_per_node = (int*) NCI_Malloc(fd->node_ids.num_nodes * sizeof(int));
        for (i = 0; i < striping_factor % fd->node_ids.num_nodes; i++)
            naggr_per_node[i] = striping_factor / fd->node_ids.num_nodes + 1;
        for (; i < fd->node_ids.num_nodes; i++)
            naggr_per_node[i] = striping_factor / fd->node_ids.num_nodes;
        avg = num_aggr / striping_factor;
        if (avg > 0)
            for (i = 0; i < fd->node_ids.num_nodes; i++)
                naggr_per_node[i] *= avg;
        for (i = 0; i < fd->node_ids.num_nodes; i++)
            naggr_per_node[i] = MIN(naggr_per_node[i], nprocs_per_node[i]);
        /* naggr_per_node[] is the number of aggregators that can be
         * selected as I/O aggregators
         */

        if (block_assignment) {
            int n = 0;
            for (j=0; j<fd->node_ids.num_nodes; j++) {
                /* j is the node ID */
                int rank_stride = nprocs_per_node[j] / naggr_per_node[j];
                /* try stride==1 seems no effect, rank_stride = 1; */
                for (k=0; k<naggr_per_node[j]; k++) {
                    fd->hints->ranklist[n] = ranks_per_node[j][k*rank_stride];
                    if (++n == num_aggr) {
                        j = fd->node_ids.num_nodes; /* break loop j */
                        break; /* loop k */
                    }
                }
            }
        }
        else {
            for (i = 0; i < num_aggr; i++) {
                int stripe_i = i % striping_factor;
                j = stripe_i % fd->node_ids.num_nodes; /* select from node j */
                k = nprocs_per_node[j] / naggr_per_node[j];
                k *= idx_per_node[j];
                /* try stride==1 seems no effect, k = idx_per_node[j]; */
                idx_per_node[j]++;
                assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
        }
        NCI_Free(naggr_per_node);
        NCI_Free(idx_per_node);
    }

    /* TODO: we can keep these two arrays in case for dynamic construction
     * of fd->hints->ranklist[], such as in group-cyclic file domain
     * assignment method, used in each collective write call.
     */
    NCI_Free(nprocs_per_node);
    NCI_Free(ranks_per_node[0]);
    NCI_Free(ranks_per_node);

    /* set file striping hints */
    fd->hints->cb_nodes = num_aggr;

    /* check whether this process is selected as an I/O aggregator */
    fd->is_agg = 0;
    fd->my_cb_nodes_index = -1;
    for (i = 0; i < num_aggr; i++) {
        if (rank == fd->hints->ranklist[i]) {
            fd->is_agg = 1;
            fd->my_cb_nodes_index = i;
            break;
        }
    }

    return 0;
}

/*----< PNCIO_Lustre_create() >----------------------------------------------*/
/*   1. root creates the file
 *   2. root sets and obtains striping info
 *   3. root broadcasts striping info
 *   4. non-root processes receive striping info from root
 *   5. non-root processes opens the file
 */
int
PNCIO_Lustre_create(PNCIO_File *fd,
                    int         mpi_io_mode)
{
    char int_str[16];
    int err=NC_NOERR, rank, perm, old_mask;
    int stripin_info[4] = {-1, -1, -1, -1};
#ifdef HAVE_LUSTRE
    int total_num_OSTs;
    uint64_t numOSTs, pattern, stripe_count, stripe_size, start_iodevice;
#endif

#ifdef WKL_DEBUG
extern int first_ost_id;
first_ost_id = -1;
#endif

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
static int wkl=0; if (wkl == 0 && world_rank == 0) { printf("\nxxxx %s at %d: %s ---- %s\n",__func__,__LINE__,(fd->file_system == PNCIO_LUSTRE)?"PNCIO_LUSTRE":"PNCIO_UFS",fd->filename); wkl++; fflush(stdout);}
#endif

#if defined(HAVE_LUSTRE) || defined(MIMIC_LUSTRE)
assert(mpi_io_mode & MPI_MODE_CREATE);

/* Note ncmpi_create always creates a file with readable and writable permission. */
    int amode = O_CREAT;
    if (mpi_io_mode & MPI_MODE_RDWR) amode |= O_RDWR;
#endif

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ PNCIO_PERM;

    /* root process creates the file first, followed by all processes open the
     * file.
     */
    if (rank > 0) goto err_out;

    /* For Lustre, we need to obtain file striping info (striping_factor,
     * striping_unit, and num_osts) in order to select the I/O aggregators
     * in fd->hints->ranklist, no matter its is open or create mode.
     */

// printf("fd->hints->pnc_striping %s\n", (fd->hints->pnc_striping == PNCIO_STRIPING_AUTO)?"AUTO":"INHERIT");

#ifdef HAVE_LUSTRE
    int overstriping_ratio, str_factor, str_unit, start_iodev;

    /* In a call to PNCIO_File_SetInfo() earlier, hints have been validated to
     * be consistent among all processes.
     */

    str_unit           = fd->hints->striping_unit;
    str_factor         = fd->hints->striping_factor;
    start_iodev        = fd->hints->start_iodevice;
    overstriping_ratio = fd->hints->lustre_overstriping_ratio;

    if (overstriping_ratio <= 0) /* hint not set of disabled */
        overstriping_ratio = 1;

    /* obtain the total number of OSTs available */
    total_num_OSTs = get_total_avail_osts(fd->filename);
    if (total_num_OSTs <= 0) /* failed to obtain number of available OSTs */
        total_num_OSTs = PNCIO_LUSTRE_MAX_OSTS;

    /* make sure str_factor <= overstriping_ratio * total_num_OSTs */
    if (str_factor > overstriping_ratio * total_num_OSTs)
        str_factor = overstriping_ratio * total_num_OSTs;

    numOSTs=0;
    pattern = LLAPI_LAYOUT_DEFAULT;
    stripe_count = LLAPI_LAYOUT_DEFAULT;
    stripe_size = LLAPI_LAYOUT_DEFAULT;
    start_iodevice = LLAPI_LAYOUT_DEFAULT;

    fd->fd_sys = -1;

    /* When no file striping hint is set, their default values are:
     * fd->hints->striping_factor = 0;
     * fd->hints->striping_unit = 0;
     * fd->hints->start_iodevice = -1;
     * fd->hints->lustre_overstriping_ratio = 1;
     */

    /* Now select file striping configuration for the new file. In many cases,
     * the Lustre striping configuration of the file to be created is not
     * explicitly set by the users (through I/O hints striping_factor and
     * striping_unit) and the striping configuration of parent folder to store
     * the new file is not explicitly set by the users.
     *
     * Codes below try to set the striping for the new file. Precedences are:
     * 1. When hints striping_factor and striping_unit are explicitly set, they
     *    are used as the top precedence.
     * 2. When hint pnc_striping is set to "inherit", the striping will inherit
     *    from the parent folder. If the parent folder's striping count is not
     *    set, then this hint is ignored.
     * 3. When no hint are set, set the new file's striping count to be equal
     *    to the number of compute nodes allocated to fd->comm and the striping
     *    size to 1 MiB.
     */
    if (fd->hints->striping_factor == 0 &&
        fd->hints->pnc_striping == PNCIO_STRIPING_INHERIT) {
        /* Inherit the file striping settings from the parent folder. */
        int dd;
        char *dirc, *dname;

        dirc = NCI_Strdup(fd->filename);
        dname = dirname(dirc); /* folder name */

        dd = open(dname, O_RDONLY, PNCIO_PERM);

        numOSTs = get_striping(dd, dname, &pattern,
                               &stripe_count,
                               &stripe_size,
                               &start_iodevice);
        close(dd);
        NCI_Free(dirc);

#ifdef PNETCDF_LUSTRE_DEBUG_VERBOSE
        printf("line %d: use parent folder's striping to set file's:\n",__LINE__);
        PRINT_LAYOUT(numOSTs);
        PRINT_LAYOUT(stripe_count);
        PRINT_LAYOUT(stripe_size);
        PRINT_LAYOUT(start_iodevice);
        PRINT_LAYOUT(pattern);
#endif
    }

    /* If hint striping_factor is not set by the user and the new file's folder
     * has not set its striping parameters, then we set the number of unique
     * OSTs, numOSTs, to the number of compute nodes allocated to this job,
     * which sets stripe_count to (numOSTs * overstriping_ratio).
     */
    if (str_factor == 0 && (stripe_count == LLAPI_LAYOUT_DEFAULT ||
                            stripe_count == LLAPI_LAYOUT_WIDE)) {
        stripe_count = MIN(fd->node_ids.num_nodes, total_num_OSTs);
        if (overstriping_ratio > 1) stripe_count *= overstriping_ratio;
    }
    else if (str_factor > 0)
        stripe_count = str_factor;

    /* When overstriping is requested by the user, calculate the number of
     * unique OSTs.
     */
    if (overstriping_ratio > 1) {
        pattern = LLAPI_LAYOUT_OVERSTRIPING;
        if (stripe_count < overstriping_ratio)
            numOSTs = 1;
        else
            numOSTs = stripe_count / overstriping_ratio;
    }
    /* If ill values are detected, fall back to no overstriping */
    if (overstriping_ratio <= 1 || numOSTs == stripe_count) {
        numOSTs = stripe_count;
        pattern = LLAPI_LAYOUT_RAID0;
    }

    /* If user has not set hint striping_unit and the folder's striping size is
     * also not set, then use the default.
     */
    if (str_unit == 0 && stripe_size == LLAPI_LAYOUT_DEFAULT)
        stripe_size = LLAPI_LAYOUT_DEFAULT;
    else if (str_unit > 0)
        stripe_size = str_unit;

    /* If user has not set hint start_iodevice and the folder's start_iodevice
     * is also not set, then use the default.
     */
    if (start_iodev == -1 && start_iodevice == LLAPI_LAYOUT_DEFAULT)
        start_iodevice = LLAPI_LAYOUT_DEFAULT;
    else if (start_iodev > 0)
        start_iodevice = start_iodev;

#ifdef PNETCDF_LUSTRE_DEBUG_VERBOSE
    printf("\n\tAfter adjust striping parameters become:\n");
    PRINT_LAYOUT(numOSTs);
    PRINT_LAYOUT(stripe_count);
    PRINT_LAYOUT(stripe_size);
    PRINT_LAYOUT(start_iodevice);
    PRINT_LAYOUT(pattern);
#endif

    /* create a new file and set striping */
    fd->fd_sys = set_striping(fd->filename, pattern,
                                            numOSTs,
                                            stripe_count,
                                            stripe_size,
                                            start_iodevice);

    if (fd->fd_sys < 0)
        /* If explicitly setting file striping failed, inherit the striping
         * from the folder by simply creating the file.
         */
        fd->fd_sys = open(fd->filename, amode, perm);

    if (fd->fd_sys < 0) {
        fprintf(stderr,"Error at %s (%d) fails to create file %s (%s)\n",
                __FILE__,__LINE__, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("Lustre set striping");
        goto err_out;
    }
    fd->is_open = 1;

    /* Obtain Lustre file striping parameters actually set. */
    numOSTs = get_striping(fd->fd_sys, fd->filename, &pattern,
                                       &stripe_count,
                                       &stripe_size,
                                       &start_iodevice);

    stripin_info[0] = stripe_size;
    stripin_info[1] = stripe_count;
    stripin_info[2] = start_iodevice;
    stripin_info[3] = numOSTs;

#elif defined(MIMIC_LUSTRE)
    fd->fd_sys = open(fd->filename, amode, perm);
    if (fd->fd_sys == -1) {
        printf("%s line %d: rank %d fails to create file %s (%s)\n",
               __FILE__,__LINE__, rank, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("open");
        goto err_out;
    }
    fd->is_open = 1;

    char *env_str = getenv("MIMIC_STRIPE_SIZE");
    if (env_str != NULL)
        stripin_info[0] = atoi(env_str);
    else
        stripin_info[0] = STRIPE_SIZE;
    stripin_info[1] = STRIPE_COUNT;
    stripin_info[2] = 0;
    stripin_info[3] = STRIPE_COUNT;
#endif

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    if (fd->file_system == PNCIO_LUSTRE &&
        (stripin_info[0] == -1 || stripin_info[3] == 0)) {
        fprintf(stderr, "%s line %d: failed to create Lustre file %s\n",
                __FILE__, __LINE__, fd->filename);
        return err;
    }

    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];
    if (fd->file_system == PNCIO_LUSTRE) {
        fd->hints->lustre_num_osts = stripin_info[3];
        fd->hints->lustre_overstriping_ratio = stripin_info[1] / stripin_info[3];
    }

    if (rank > 0) { /* non-root processes */
        fd->fd_sys = open(fd->filename, O_RDWR, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d failure to open file %s (%s)\n",
                    __FILE__,__LINE__, rank, fd->filename, strerror(errno));
            return ncmpii_error_posix2nc("ioctl");
        }
        fd->is_open = 1;
    }

    /* construct cb_nodes rank list */
    Lustre_set_cb_node_list(fd);

    MPI_Info_set(fd->info, "romio_filesystem_type", "LUSTRE:");

    snprintf(int_str, 16, "%d", fd->hints->lustre_num_osts);
    MPI_Info_set(fd->info, "lustre_num_osts", int_str);

    snprintf(int_str, 16, "%d", fd->hints->lustre_overstriping_ratio);
    MPI_Info_set(fd->info, "lustre_overstriping_ratio", int_str);

    return err;
}

/*----< PNCIO_Lustre_open() >------------------------------------------------*/
/*   1. all processes open the file.
 *   2. root obtains striping info and broadcasts to all others
 */
int
PNCIO_Lustre_open(PNCIO_File *fd)
{
    char int_str[16];
    int err=NC_NOERR, rank, perm, old_mask;
    int stripin_info[4] = {1048576, -1, -1, -1};

#ifdef WKL_DEBUG
extern int first_ost_id;
first_ost_id = -1;
#endif

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
static int wkl=0; if (wkl == 0 && world_rank == 0) { printf("\nxxxx %s at %d: %s ---- %s\n",__func__,__LINE__,(fd->file_system == PNCIO_LUSTRE)?"PNCIO_LUSTRE":"PNCIO_UFS",fd->filename); wkl++; fflush(stdout);}
#endif

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ PNCIO_PERM;

    int omode = (fd->access_mode & MPI_MODE_RDWR) ? O_RDWR : O_RDONLY;

    /* All processes open the file. */
    fd->fd_sys = open(fd->filename, omode, perm);
    if (fd->fd_sys == -1) {
        fprintf(stderr,"%s line %d: rank %d fails to open file %s (%s)\n",
                __FILE__,__LINE__, rank, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("open");
        goto err_out;
    }
    fd->is_open = 1;

    /* Only root obtains the striping information and bcast to all other
     * processes.
     */
    if (rank == 0) {
#ifdef HAVE_LUSTRE
        uint64_t numOSTs=0;
        uint64_t pattern = LLAPI_LAYOUT_DEFAULT;
        uint64_t stripe_count = LLAPI_LAYOUT_DEFAULT;
        uint64_t stripe_size = LLAPI_LAYOUT_DEFAULT;
        uint64_t start_iodevice = LLAPI_LAYOUT_DEFAULT;

        numOSTs = get_striping(fd->fd_sys, fd->filename, &pattern,
                                           &stripe_count,
                                           &stripe_size,
                                           &start_iodevice);

        stripin_info[0] = stripe_size;
        stripin_info[1] = stripe_count;
        stripin_info[2] = start_iodevice;
        stripin_info[3] = numOSTs;

#elif defined(MIMIC_LUSTRE)
        char *env_str = getenv("MIMIC_STRIPE_SIZE");
        if (env_str != NULL)
            stripin_info[0] = atoi(env_str);
        else
            stripin_info[0] = STRIPE_SIZE;
        stripin_info[1] = STRIPE_COUNT;
        stripin_info[2] = 0;
        stripin_info[3] = STRIPE_COUNT;
#endif
    }

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];
    fd->hints->lustre_num_osts = stripin_info[3];
    fd->hints->lustre_overstriping_ratio = stripin_info[1] / stripin_info[3];

    /* construct cb_nodes rank list */
    Lustre_set_cb_node_list(fd);

    MPI_Info_set(fd->info, "romio_filesystem_type", "LUSTRE:");

    snprintf(int_str, 16, "%d", fd->hints->lustre_num_osts);
    MPI_Info_set(fd->info, "lustre_num_osts", int_str);

    snprintf(int_str, 16, "%d", fd->hints->lustre_overstriping_ratio);
    MPI_Info_set(fd->info, "lustre_overstriping_ratio", int_str);

    return err;
}

