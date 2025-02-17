/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <fcntl.h>

#include <pncio.h>

static
const char *GEN_flock_cmd_to_string(int cmd)
{
    switch (cmd) {
#ifdef F_GETLK64
        case F_GETLK64:
            return "F_GETLK64";
#else
        case F_GETLK:
            return "F_GETLK";
#endif
#ifdef F_SETLK64
        case F_SETLK64:
            return "F_SETLK64";
#else
        case F_SETLK:
            return "F_SETLK";
#endif
#ifdef F_SETLKW64
        case F_SETLKW64:
            return "F_SETLKW64";
#else
        case F_SETLKW:
            return "F_SETLKW";
#endif
        default:
            return "UNEXPECTED";
    }
}

static
const char *GEN_flock_type_to_string(int type)
{
    switch (type) {
        case F_RDLCK:
            return "F_RDLCK";
        case F_WRLCK:
            return "F_WRLCK";
        case F_UNLCK:
            return "F_UNLOCK";
        default:
            return "UNEXPECTED";
    }
}

int PNCIO_GEN_SetLock(PNCIO_File *fd, int cmd, int type, MPI_Offset offset, int whence,
                      MPI_Offset len)
{
    FDTYPE fd_sys = fd->fd_sys;
    int err, error_code, err_count = 0, sav_errno;
    struct flock lock;

    if (len == 0)
        return MPI_SUCCESS;


    /* Depending on the compiler flags and options, struct flock
     * may not be defined with types that are the same size as
     * MPI_Offsets.  */
/* FIXME: This is a temporary hack until we use flock64 where
   available. It also doesn't fix the broken Solaris header sys/types.h
   header file, which declares off_t as a UNION ! Configure tests to
   see if the off64_t is a union if large file support is requested;
   if so, it does not select large file support.
*/
#ifdef NEEDS_INT_CAST_WITH_FLOCK
    lock.l_type = type;
    lock.l_start = (int) offset;
    lock.l_whence = whence;
    lock.l_len = (int) len;
#else
    lock.l_type = type;
    lock.l_whence = whence;
    lock.l_start = offset;
    lock.l_len = len;
#endif

    sav_errno = errno;  /* save previous errno in case we recover from retryable errors */
    errno = 0;
    do {
        err = fcntl(fd_sys, cmd, &lock);
    } while (err && ((errno == EINTR) || ((errno == EINPROGRESS) && (++err_count < 10000))));

    if (err && (errno != EBADF)) {
        /* FIXME: This should use the error message system,
         * especially for MPICH */
        fprintf(stderr,
                "This requires fcntl(2) to be implemented. As of 8/25/2011 it is not. Generic MPICH Message: File locking failed in PNCIO_GEN_SetLock(fd %X,cmd %s/%X,type %s/%X,whence %X) with return value %X and errno %X.\n"
                "- If the file system is NFS, you need to use NFS version 3, ensure that the lockd daemon is running on all the machines, and mount the directory with the 'noac' option (no attribute caching).\n"
                "- If the file system is LUSTRE, ensure that the directory is mounted with the 'flock' option.\n",
                fd_sys, GEN_flock_cmd_to_string(cmd), cmd,
                GEN_flock_type_to_string(type), type, whence, err, errno);
        perror("PNCIO_GEN_SetLock:");
        fprintf(stderr, "PNCIO_GEN_SetLock:offset %llu, length %llu\n", (unsigned long long) offset,
                (unsigned long long) len);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!err)   /* report fcntl failure errno's (EBADF), otherwise */
        errno = sav_errno;      /* restore previous errno in case we recovered from retryable errors */

    error_code = (err == 0) ? MPI_SUCCESS : MPI_ERR_UNKNOWN;
    return error_code;
}

int PNCIO_GEN_SetLock64(PNCIO_File *fd, int cmd, int type, MPI_Offset offset, int whence,
                        MPI_Offset len)
{
    FDTYPE fd_sys = fd->fd_sys;
    int err, error_code;
#ifdef _LARGEFILE64_SOURCE
    struct flock64 lock;
#else
    struct flock lock;
#endif

    if (len == 0)
        return MPI_SUCCESS;

    lock.l_type = type;
    lock.l_start = offset;
    lock.l_whence = whence;
    lock.l_len = len;

    do {
        err = fcntl(fd_sys, cmd, &lock);
    } while (err && (errno == EINTR));

    if (err && (errno != EBADF)) {
        fprintf(stderr,
                "File locking failed in PNCIO_GEN_SetLock64(fd %X,cmd %s/%X,type %s/%X,whence %X) with return value %X and errno %X.\n"
                "If the file system is NFS, you need to use NFS version 3, ensure that the lockd daemon is running on all the machines, and mount the directory with the 'noac' option (no attribute caching).\n",
                fd_sys, GEN_flock_cmd_to_string(cmd), cmd,
                GEN_flock_type_to_string(type), type, whence, err, errno);
        perror("PNCIO_GEN_SetLock64:");
        fprintf(stderr, "PNCIO_GEN_SetLock:offset %llu, length %llu\n", (unsigned long long) offset,
                (unsigned long long) len);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    error_code = (err == 0) ? MPI_SUCCESS : MPI_ERR_UNKNOWN;
    return error_code;
}
