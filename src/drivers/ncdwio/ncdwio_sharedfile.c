/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pnc_debug.h>
#include <common.h>
#include <pnetcdf.h>
#include <ncdwio_driver.h>

#define BLOCK_SIZE 8388608

/*
 * Open shared file
 * IN      comm:    MPI communicator of processes sharing the file
 * IN      path:    Path of the file to open
 * IN      amode:   File open mode (using POSIX constants)
 * IN      info:    File hint for opened file (currently unused)
 * OUT       fd:    File handler
 */
int ncdwio_sharedfile_open(MPI_Comm comm, char *path, int flag, MPI_Info info, NC_dw_sharedfile **fh) {
    int err;
    NC_dw_sharedfile *f;

    // Allocate file handle
    f = (NC_dw_sharedfile*)NCI_Malloc(sizeof(NC_dw_sharedfile));

    /* Initialize metadata associate with the file
     * We assum all processes within the given communicator is sharing the file
     * Due to file sharing, actual file position differs than the logical file position within the file view
     * TODO: Adjustable bsize
     */
    f->bsize = BLOCK_SIZE;
    f->pos = 0;
    f->fsize = 0;
    MPI_Comm_rank(comm, &(f->chanel));
    MPI_Comm_size(comm, &(f->nchanel));

    /* Open file
     * Process using chanel 0 is responsible to create the file if mode is create
     * We do this by letting processes using chanel 0 open the file first with given amode
     * Other process then follows with open flag removed
     */

    // Process using chanel 0 creates the file
    if (f->chanel == 0){
        f->fd = open(path, flag, 0744);
        if (f->fd < 0){
            NCI_Free(f);
            err = ncmpii_error_posix2nc("open");
            DEBUG_RETURN_ERROR(err);
        }
    }

    // Wait for file open if there is more than one process sharing the file
    if (f->nchanel > 1){
        MPI_Barrier(comm);
    }

    // Remaining process opens the file
    if (f->chanel > 0){
        // Remove open flag if it's set, the first process already have the file opened
        f->fd = open(path, flag & (~(O_CREAT)), 0744);
        if (f->fd < 0){
            NCI_Free(f);
            err = ncmpii_error_posix2nc("open");
            DEBUG_RETURN_ERROR(err);
        }
    }

    *fh = f;
    return NC_NOERR;
}

/*
 * Close shared file
 * OUT       fd:    File handler
 */
int ncdwio_sharedfile_close(NC_dw_sharedfile *f) {
    int err;

    /* Close posix file descriptor */
    err = close(f->fd);
    if (err != 0){
        err = ncmpii_error_posix2nc("close");
        DEBUG_RETURN_ERROR(err);
    }

    // Free the file handle
    NCI_Free(f);

    return NC_NOERR;
}

/*
 * This function write <count> bytes of data in buffer <buf> to the file at position
 * specified by <offset>
 * The file position is not changed
 * IN       f:    File handle
 * IN     buf:    Buffer of data to be written
 * IN   count:    Number of bytes to write
 * IN  offset:    Starting write position
 *
 * A continuous region in the logical file view seen by upper layer is not continuous in the shared file
 * As a result, one shared file write operation translates to multiple writes
 * We first divide the writing region in the logical file into blocks
 * Each block in the logical file view is then mapped to corresponding block in the shared file
 * First and last block may be partial, we need to adjust write region of the first block and the last block
 *
 * Figure showing the case using chanel 0 with 2 chanels
 * logical file view: |B0      |B1      |B2      |B3      |...
 * write region:           |B0 |B1      |B2   |
 *                         ^                  ^
 *                       offset        offset + count
 * Block map: B0 -> B0, B1 -> B2, B2 -> B4 ...
 * Physical file:     |B0  ----|B1      |B2------|B3      |B4---      |B5      |B6      |
 *                         ^   ^        ^        ^        ^     ^
 *               offstart(s) offend(e) (s)      (e)      (s)   (e)
 * Dashed line shows write coverage
 */
int ncdwio_sharedfile_pwrite(NC_dw_sharedfile *f, void *buf, size_t count, off_t offset){
    int i, err;
    int sblock, eblock;   // start and end block
    off_t offstart, offend; // Start and end offset to write for current block in physical file
    ssize_t ioret;

    // Write directly if not sharing
    if(f->nchanel == 1){
        ioret = pwrite(f->fd, buf, count, offset);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("write");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != count){
            DEBUG_RETURN_ERROR(NC_EWRITE);
        }
        return NC_NOERR;
    }

    /* Calculate first and last blocks
     * Offset / Block size = Block number
     */
    sblock = offset / f->bsize;
    eblock = (offset + count) / f->bsize;

    /* Write each block to the file
     * We first compute the mapped offset of current block in the shared file
     * The count is set to blocksize
     * For the first block, start offset is increased by the offset within the block to handle partial block
     * For last block, end offset must be adjusted
     * After adjusting start and end offset, we write the block to the correct location
     * Local blocknumber * nchanel = global blocknumber
     * global blocknumber * blocksize = global offset
     * Offset % Block size = Offset within the block
     */
    for(i = sblock; i <= eblock; i++){
        /* Block boundary are inclusive, as a result, final block is always partial (can be empty)
         * In this case, we can assume offend will always be larger than offstart
         */
        // Compute physical offset of th eblock
        offstart = i * f->nchanel * f->bsize;
        // A block can be first and last block at the same time due to short write region
        // Last block must be partial
        // NOTE: offend must be computed before offstart, we reply on unadjusted offstart to mark the start position of the block
        if (i == eblock){
            // The local block offset of the end of writing region is the write size of final block
            offend = offstart + (offset + count) % f->bsize;
        }
        else{
            offend = offstart + f->bsize;
        }
        // First block can be partial
        if (i == sblock){
            // The local block offset of the start of writing region is the amount to skip for the first block
            offstart += offset % f->bsize;
        }

        // Write to file
        ioret = pwrite(f->fd, buf, offend - offstart, offstart);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("write");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != offend - offstart){
            DEBUG_RETURN_ERROR(NC_EWRITE);
        }

        // We increase the buffer pointer by amount writen, so it always points to the data of next block
        buf = (void*)(((char*)buf) + offend - offstart);
    }

    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < offset + count){
        f->fsize = offset + count;
    }

    return NC_NOERR;
}

/*
 * This function write <count> bytes of data in buffer <buf> to the file at it's
 * current file position and increase the file position by <count>
 * IN       f:    File handle
 * IN     buf:    Buffer of data to be written
 * IN   count:    Number of bytes to write
 *
 * We call ncdwio_sharedfile_pwrite and then increase the file position by count
 */
int ncdwio_sharedfile_write(NC_dw_sharedfile *f, void *buf, size_t count){
    int err;

    // Write directly if not sharing
    if(f->nchanel == 1){
        ssize_t ioret;
        ioret = write(f->fd, buf, count);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("write");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != count){
            DEBUG_RETURN_ERROR(NC_EWRITE);
        }
    }
    else{
        // Write at current file position
        err = ncdwio_sharedfile_pwrite(f, buf, count, f->pos);
        if (err != NC_NOERR){
            return err;
        }
    }

    // Increase current file position
    f->pos += count;

    return NC_NOERR;
}

/*
 * This function read <count> bytes of data to buffer <buf> from the file at position
 * specified by <offset>
 * The file position is not changed
 * IN       f:    File handle
 * OUT    buf:    Buffer of data to be written
 * IN   count:    Number of bytes to read
 * IN  offset:    Starting read position
 *
 * A continuous region in the logical file view seen by upper layer is not continuous in the shared file
 * As a result, one shared file read operation translates to multiple reads
 * We first divide the reading region in the logical file into blocks
 * Each block in the logical file view is then mapped to corresponding block in the shared file
 * First and last block may be partial, we need to adjust read region of the first block and the last block
 *
 * Figure showing the case using chanel 0 with 2 chanels
 * logical file view: |B0      |B1      |B2      |B3      |...
 * read region:           |B0 |B1      |B2   |
 *                         ^                  ^
 *                       offset        offset + count
 * Block map: B0 -> B0, B1 -> B2, B2 -> B4 ...
 * Physical file:     |B0  ----|B1      |B2------|B3      |B4---      |B5      |B6      |
 *                         ^   ^        ^        ^        ^     ^
 *               offstart(s) offend(e) (s)      (e)      (s)   (e)
 * Dashed line shows write coverage
 */
int ncdwio_sharedfile_pread(NC_dw_sharedfile *f, void *buf, size_t count, off_t offset){
    int i, err;
    int sblock, eblock;   // start and end block
    off_t offstart, offend; // Start and end offset to write for current block in physical file
    ssize_t ioret;

    // Read directly if not sharing
    if(f->nchanel == 1){
        ioret = pread(f->fd, buf, count, offset);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("read");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != count){
            DEBUG_RETURN_ERROR(NC_EREAD);
        }
        return NC_NOERR;
    }

    /* Calculate first and last blocks
     * Offset / Block size = Block number
     */
    sblock = offset / f->bsize;
    eblock = (offset + count) / f->bsize;

    /* Write each block to the file
     * We first compute the mapped offset of current block in the shared file
     * The count is set to blocksize
     * For the first block, start offset is increased by the offset within the block to handle partial block
     * For last block, end offset must be adjusted
     * After adjusting start and end offset, we write the block to the correct location
     * Local blocknumber * nchanel = global blocknumber
     * global blocknumber * blocksize = global offset
     * Offset % Block size = Offset within the block
     */
    for(i = sblock; i <= eblock; i++){
        /* Block boundary are inclusive, as a result, final block is always partial (can be empty)
         * In this case, we can assume offend will always be larger than offstart
         */
        // Compute physical offset of th eblock
        offstart = i * f->nchanel * f->bsize;
        // A block can be first and last block at the same time due to short write region
        // Last block must be partial
        // NOTE: offend must be computed before offstart, we reply on unadjusted offstart to mark the start position of the block
        if (i == eblock){
            // The local block offset of the end of writing region is the write size of final block
            offend = offstart + (offset + count) % f->bsize;
        }
        else{
            offend = offstart + f->bsize;
        }
        // First block can be partial
        if (i == sblock){
            // The local block offset of the start of writing region is the amount to skip for the first block
            offstart += offset % f->bsize;
        }

        // Read from file
        ioret = pread(f->fd, buf, offend - offstart, offstart);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("read");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != offend - offstart){
            DEBUG_RETURN_ERROR(NC_EREAD);
        }

        // We increase the buffer pointer by amount readn, so it always points to the location for next block
        buf = (void*)(((char*)buf) + offend - offstart);
    }

    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < offset + count){
        f->fsize = offset + count;
    }

    return NC_NOERR;
}

/*
 * This function read <count> bytes of data to buffer <buf> from the file at it's
 * current file position and increase the file position by <count>
 * IN       f:    File handle
 * OUT    buf:    Buffer of data to be written
 * IN   count:    Number of bytes to write
 *
 * We call ncdwio_sharedfile_pread and then increase the file position by count
 */
int ncdwio_sharedfile_read(NC_dw_sharedfile *f, void *buf, size_t count){
    int err;

    // Read directly if not sharing
    if(f->nchanel == 1){
        ssize_t ioret;
        ioret = read(f->fd, buf, count);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("read");
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD);
            DEBUG_RETURN_ERROR(err);
        }
        if (ioret != count){
            DEBUG_RETURN_ERROR(NC_EREAD);
        }
    }
    else{
        // Read from current file position
        err = ncdwio_sharedfile_pread(f, buf, count, f->pos);
        if (err != NC_NOERR){
            return err;
        }
    }

    // Increase current file position
    f->pos += count;

    return NC_NOERR;
}

/*
 * This function change the current file position according to <offset and whence>
 * IN       f:    File handle
 * IN  offset:    New offset
 * IN  whence:    Meaning of new offset
 */
int ncdwio_sharedfile_seek(NC_dw_sharedfile *f, off_t offset, int whence){
    int err;
    off_t ioret;

    /* Move file pointer if not sharing, so write and read can function properly without doing pwrite/read
     * Logical file position is not tracked in this case because it's equal to physical position
     */
    if(f->nchanel == 1){
        ioret = lseek(f->fd, offset, whence);
        if (ioret < 0){
            err = ncmpii_error_posix2nc("lseek");
            DEBUG_RETURN_ERROR(err);
        }
    }
    else{
        // Update file position
        switch (whence){
            case SEEK_SET:  // Offset from begining of the file
                f->pos = offset;
                break;
            case SEEK_CUR:  // Offset related to current position
                f->pos += offset;
                break;
            case SEEK_END: // Offset from end of the file
                f->pos = f->fsize + offset;
                break;
            default:
                DEBUG_RETURN_ERROR(NC_EFILE);
        }
    }

    return NC_NOERR;
}
