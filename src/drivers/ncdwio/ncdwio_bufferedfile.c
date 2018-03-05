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

#define BUFSIZE 8388608

/*
 * Open buffered file
 * IN      comm:    MPI communicator of processes sharing the file
 * IN      path:    Path of the file to open
 * IN      amode:   File open mode (using POSIX constants)
 * IN      info:    File hint for opened file (currently unused)
 * OUT       fd:    File handler
 */
int ncdwio_bufferedfile_open(MPI_Comm comm, char *path, int flag, MPI_Info info, NC_dw_bufferedfile **fh) {
    int err; 
    NC_dw_bufferedfile *f;

    /* Allocate buffer */
    f = (NC_dw_bufferedfile*)NCI_Malloc(sizeof(NC_dw_bufferedfile));

    /* Initialize metadata associate with the file
     * We assum file system block (stripe) size is equal to buffer size
     * TODO: Adjustable bsize
     */
    f->bsize = BUFSIZE;
    f->pos = 0;
    f->bused = 0;
    f->bunused = 0;
    f->fsize = 0;
    // Allocate the buffer
    if (f->bsize > 0){
        f->buffer = NCI_Malloc(BUFSIZE);
        if (f->buffer == NULL){
            DEBUG_RETURN_ERROR(NC_ENOMEM);
        }
    }
    else{
        f->buffer = NULL;
    }
    
    /* Open file */
    err = ncdwio_sharedfile_open(comm, path, flag, info, &(f->fd));
    if (err != NC_NOERR){
        NCI_Free(f->buffer);
        NCI_Free(f);
        return err;
    }

    *fh = f;
    return NC_NOERR;
}

/*
 * Close buffered file
 * OUT       fd:    File handler
 */
int ncdwio_bufferedfile_close(NC_dw_bufferedfile *f) {
    int err;

    /* Close file */
    err = ncdwio_sharedfile_close(f->fd);
    if (err != 0){
        err = ncmpii_error_posix2nc("close");
        DEBUG_RETURN_ERROR(err);
    }

    // Free the file handle
    if (f->buffer != NULL){
        NCI_Free(f->buffer);
    }
    NCI_Free(f);

    return NC_NOERR;
}

/*
 * This function write <count> bytes of data in buffer <buf> to the file at it's 
 * current file position and increase the file position by <count>
 * IN       f:    File handle
 * IN     buf:    Buffer of data to be written
 * IN   count:    Number of bytes to write
 * 
 * We always left the physical file possition on the block boundary after a write operation
 * We keep track of logical file position including data writen to buffer
 * Buffer maps to a continuous region in the file space that is aligned to block boundary
 * Buffer is used to store the partial block from the tail of previous write
 * Whenever there is a seek, we flush the buffer to ensure correctness
 * 
 * Given a write region, we divide it into 3 parts
 * A head section where the start is not on the block boundary
 * A mid section which is perfectly aligned
 * A tail section where end is not on the boundary
 * If the region is content whthin a block, the entire region is considered tail
 * If there are any data in the buffer, we combine them and write to disk
 * Because we flush the buffer after each seek operation, we are sure the data in the buffer is the data right before the head section
 * We then write out mid section as usual
 * Finally, we place the tail section in the buffer
 * 
 * Figure showing the case writing a region
 * | indicate block and write region boundary
 * File space:           |0123|4567|89AB|...
 * 
 *                 midstart (2)    midend (6)
 *                         0  v    v  count
 * Write region:           |23|4567|89|   
 *                         ^          ^
 *                      position   position + count
 *
 * Buffer status:        |01  |
 * Data written to file: |0123|4567|
 * Buffer after write:             |89  |
 */
int ncdwio_bufferedfile_write(NC_dw_bufferedfile *f, void *buf, size_t count){
    int err;
    size_t midstart, midend;    // Start and end offset of the mid section related the file position

    if (f->buffer != NULL){
        /* 
        * The start position of mid section can be calculated as the first position on the block boundary after current file position
        * The end position of mid section can be calculated as the last position on the block boundary within the write region
        * This can be incorrect when write region sits within a block where we will have start > end
        * In this case, we simply set start and end to 0, giving the entire region as tail
        */
        midstart = f->pos + (f->bsize - f->pos % f->bsize) % f->bsize;
        midend = f->pos + count - (f->pos + count) % f->bsize;
        if (midstart > midend) {
            midstart = 0;
            midend = 0;
        }
        else{
            midstart -= f->pos;
            midend -= f->pos;
        }

        /*
        * If there are data in the buffer, we combine the head section with the data and write to the file
        * Because we flush the buffer after each seek operation, we are sure the data in the buffer is the data right before the head section
        * Also, this is guaranteed to fill up the entire buffer since the region of head (if exists) end at block boundary
        * If there are no data in the buffer, we have no chioce but to write the head section directly to the file
        * This is done by merging the head section into the mid section
        * Head section may not exists due to perfectly aligned write region or short region
        */
        if (midstart > 0){
            if (f->bused - f->bunused > 0){
                // Copy data into buffer
                memcpy(f->buffer + f->bused, buf, midstart); 
                // Write combined data to file, skipping unused part
                err = ncdwio_sharedfile_write(f->fd, f->buffer + f->bunused, f->bsize - f->bunused);
                if (err != NC_NOERR){
                    return err;
                }
                // Reset buffer status
                f->bused = 0;
                f->bunused = 0;
            }
            else{
                // Treat head section as mid section if nothing to combine
                midstart = 0;
            }
        }
        
        /*
        * Write mid section as usual
        * Mid section may not exists due to short write region
        */
        if (midend > midstart) {
            err = ncdwio_sharedfile_write(f->fd, (void*)(((char*)buf) + midstart), midend - midstart); 
            if (err != NC_NOERR){
                return err;
            }
        }

        /*
        * Place the tail section to the buffer
        * This will not fill up the buffer since the region of tail section (if exists) never ends at block boundary
        * Tail section may not exists due to perfectly aligned write region
        */
        if (count > midend){
            memcpy((void*)(((char*)f->buffer) + f->bused), (void*)(((char*)buf) + midend), count - midend);
            f->bused += count - midend;
        }
    }
    else{
        // Write to file directly if buffer is disabled
        err = ncdwio_sharedfile_write(f->fd, buf, count);
        if (err != NC_NOERR){
            return err;
        }
    }

    // Increase file position
    f->pos += count;

    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < f->pos){
        f->fsize = f->pos;
    }

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
 * pwrite is not buffered, we write directly to the file
 */
int ncdwio_bufferedfile_pwrite(NC_dw_bufferedfile *f, void *buf, size_t count, off_t offset){
    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < offset + count){
        f->fsize = offset + count;
    }

    // Write directly
    return ncdwio_sharedfile_pwrite(f->fd, buf, count, offset);
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
 * We flush the buffer before read so data in the buffer can be reflected
 * Read is not buffered, we read directly from the file
 */
int ncdwio_bufferedfile_pread(NC_dw_bufferedfile *f, void *buf, size_t count, off_t offset){    
    int err;

    if (f->buffer != NULL){
        // Flush the buffer
        if (f->bused - f->bunused > 0){
            // Write data to file, skipping unused part
            err = ncdwio_sharedfile_write(f->fd, f->buffer + f->bunused, f->bsize - f->bunused);
            if (err != NC_NOERR){
                return err;
            }
            // Reset buffer status
            f->bused = 0;
            f->bunused = 0;
        }
    }

    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < offset + count){
        f->fsize = offset + count;
    }

    // Read directly
    return ncdwio_sharedfile_pread(f->fd, buf, count, offset);
}

/*
 * This function read <count> bytes of data to buffer <buf> from the file at it's 
 * current file position and increase the file position by <count>
 * IN       f:    File handle
 * OUT    buf:    Buffer of data to be written
 * IN   count:    Number of bytes to write
 * 
 * We call ncdwio_bufferedfile_pread and then increase the file position by count
 */
int ncdwio_bufferedfile_read(NC_dw_bufferedfile *f, void *buf, size_t count){
    int err;

    // Read from current file position
    err = ncdwio_bufferedfile_pread(f, buf, count, f->pos);
    if (err != NC_NOERR){
        return err;
    }

    // Increase current file position
    f->pos += count;

    // Record the file size as the largest location ever reach by IO operation
    if (f->fsize < f->pos){
        f->fsize = f->pos;
    }

    return NC_NOERR;
}

/*
 * This function change the current file position according to <offset and whence>
 * IN       f:    File handle
 * IN  offset:    New offset
 * IN  whence:    Meaning of new offset
 */
int ncdwio_bufferedfile_seek(NC_dw_bufferedfile *f, off_t offset, int whence){
    int err;

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

    if (f->buffer != NULL){
        /* Flush the buffer
        * We assume buffer always maps to the block containing current position
        * When current position changes, buffer must be flushed
        */
        if (f->bused - f->bunused > 0){
            // Write data to file, skipping unused part
            err = ncdwio_sharedfile_write(f->fd, f->buffer + f->bunused, f->bsize - f->bunused);
            if (err != NC_NOERR){
                return err;
            }
        }

        /* Adjust buffer status
        * Buffer always maps to aligned block
        * If new position is not on the boundary, we place the buffer at the first aligned position right before current position
        * To maintain continuity of buffer data and current file position, we need to insert padding to make up the gap
        * We do this by pretending they are used
        * We keep track of this padding as (actually) unused space to prevent flushing paddings into the file
        */
        f->bused = f->pos % f->bsize;
        f->bunused = f->bused;
    }

    // Seek physical file position
    err = ncdwio_sharedfile_seek(f->fd, offset, whence);
    if (err != NC_NOERR){
        return err;
    }

    return NC_NOERR;
}
