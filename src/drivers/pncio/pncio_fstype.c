/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>   /* readlink() */
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>
#include <fcntl.h>      /* open(), O_CREAT */
#include <sys/types.h>  /* open() */
#include <libgen.h>     /* basename() */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#ifndef PATH_MAX
#define PATH_MAX 65535
#endif

#ifdef HAVE_SYS_VFS_H
#include <sys/vfs.h>
#endif
#ifdef HAVE_SYS_STATVFS_H
#include <sys/statvfs.h>
#endif
#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* struct statfs */
#endif
#ifdef HAVE_SYS_MOUNT_H
#include <sys/mount.h> /* struct statfs */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h> /* open(), fstat(), lstat(), stat() */
#endif

#include <mpi.h>

#include "pncio.h"

/* In a strict ANSI environment, S_ISLNK may not be defined. Fix that here.
 * We assume that S_ISLNK is *always* defined as a macro. If that is not
 * universally true, then add a test to the configure that tries to link
 * a program that references S_ISLNK
 */
#if !defined(S_ISLNK)
#if defined(S_IFLNK)
/* Check for the link bit */
#define S_ISLNK(mode) ((mode) & S_IFLNK)
#else
/* no way to check if it is a link, so say false */
#define S_ISLNK(mode) 0
#endif
#endif /* !(S_ISLNK) */

/* Returns a string, the parent directory of a given filename.
 * The caller should free the memory located returned by this subroutine.
 */
static
void parentdir(const char *filename, char **dirnamep)
{
    int err;
    char *dir = NULL, *slash;
    struct stat statbuf;

    err = lstat(filename, &statbuf);

    if (err || (!S_ISLNK(statbuf.st_mode))) {
        /* No such file, or file is not a link; these are the "normal" cases
         * where we can just return the parent directory.
         */
        dir = NCI_Strdup(filename);
    } else {
        /* filename is a symlink. We've presumably already tried to stat it
         * and found it to be missing (dangling link), but this code doesn't
         * care if the target is really there or not.
         */
        ssize_t namelen;
        char *linkbuf;

        linkbuf = NCI_Malloc(PATH_MAX + 1);
        namelen = readlink(filename, linkbuf, PATH_MAX + 1);
        if (namelen == -1) {
            /* Something strange has happened between the time that we
             * determined that this was a link and the time that we attempted
             * to read it; punt and use the old name.
             */
            dir = NCI_Strdup(filename);
        } else {
            /* successfully read the link */
            linkbuf[namelen] = '\0';    /* readlink doesn't null terminate */
            dir = NCI_Strdup(linkbuf);
        }
        NCI_Free(linkbuf);
    }

    slash = strrchr(dir, '/');
    if (!slash)
        strncpy(dir, ".", 2);
    else {
        if (slash == dir)
            *(dir + 1) = '\0';
        else
            *slash = '\0';
    }

    *dirnamep = dir;
    return;
}

#define UNKNOWN_SUPER_MAGIC (0xDEADBEEF)
#ifndef LL_SUPER_MAGIC
#define LL_SUPER_MAGIC 0x0BD00BD0
#endif

static int check_statfs(const char *filename, int64_t * file_id)
{
    int err = 0;

#ifdef HAVE_STRUCT_STATVFS_WITH_F_BASETYPE
    /* rare: old solaris machines */
    struct statvfs vfsbuf;
#endif
#if defined(HAVE_STRUCT_STATFS_F_TYPE) || defined(HAVE_STRUCT_STATFS_F_FSTYPENAME)
    /* common fs-detection logic for any modern POSIX-compliant environment,
     * with the one wrinkle that some platforms (Darwin, BSD) give us a file
     * system as a string, not an identifier */
    struct statfs fsbuf;
#endif

    *file_id = UNKNOWN_SUPER_MAGIC;

#ifdef HAVE_STRUCT_STATVFS_WITH_F_BASETYPE
    err = statvfs(filename, &vfsbuf);
    if (err == 0)
        *file_id = vfsbuf.f_basetype;
#endif

    /* remember above how I said 'statfs with f_type' was the common linux-y
     * way to report file system type?  Darwin (and probably the BSDs) *also*
     * uses f_type but it is "reserved" and does not give us anything
     * meaningful.  Fine.  If configure detects f_type we'll use it here and on
     * those "reserved" platforms we'll ignore that result and check the
     * f_fstypename field.
     */

#ifdef HAVE_STRUCT_STATFS_F_TYPE
    err = statfs(filename, &fsbuf);
    if (err == 0) {
        *file_id = fsbuf.f_type;
        return 0;
    }
#endif

#ifdef HAVE_STRUCT_STATFS_F_FSTYPENAME
    /* these stat routines store the file system type in a string */
    err = statfs(filename, &fsbuf);
    if (err == 0 && !strncasecmp(fsbuf.f_fstypename, "lustre", 6)) {
        *file_id = LL_SUPER_MAGIC;
        return 0;
    }
#endif

#ifdef HAVE_STRUCT_STAT_ST_FSTYPE
    struct stat sbuf;
    err = stat(filename, &sbuf);
    if (err == 0) {
        *file_id = sbuf.st_fstype;
        return 0;
    }
#endif
    return err;
}

/* Check if file system type from file name, using a system-dependent function
 * call.
 */
int PNCIO_FileSysType(const char *filename)
{

    int err, retry_cnt;
    int64_t file_id=UNKNOWN_SUPER_MAGIC;

    char *colon = strchr(filename, ':');
    if (colon != NULL) { /* there is a prefix end with : */
        if (!strncmp(filename, "lustre", 6))
            return PNCIO_LUSTRE;
        else if (!strncmp(filename, "ufs", 3))
            return PNCIO_UFS;
        else
            return 0;
    }
#ifdef MIMIC_LUSTRE
    return PNCIO_LUSTRE;
#endif

    /* NFS can get stuck and end up returning ESTALE "forever" */

#define MAX_ESTALE_RETRY 10000

    retry_cnt = 0;
    do {
        err = check_statfs(filename, &file_id);
    } while (err && (errno == ESTALE) && retry_cnt++ < MAX_ESTALE_RETRY);

    if (err) {
        /* ENOENT may be returned in two cases:
         * 1) no directory entry for "filename"
         * 2) "filename" is a dangling symbolic link
         *
         * parentdir() tries to deal with both cases.
         */
        if (errno == ENOENT) {
            char *dir;
            parentdir(filename, &dir);
            err = check_statfs(dir, &file_id);
            NCI_Free(dir);
        } else
            return 0;
    }

    if (file_id == LL_SUPER_MAGIC)
        return PNCIO_LUSTRE;
    else
        return PNCIO_UFS; /* UFS support if we don't know what else to use */
}

