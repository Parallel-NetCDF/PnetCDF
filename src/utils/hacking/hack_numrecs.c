/*
 *  Copyright (C) 2022, Northwestern University.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strerror() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>   /* open() */
#include <fcntl.h>      /* open() */
#include <unistd.h>     /* read() */
#include <inttypes.h>   /* check for Endianness, uint32_t*/
#include <errno.h>      /* errno */

#define SWAP4B(a) ( ((a) << 24) | \
                   (((a) <<  8) & 0x00ff0000) | \
                   (((a) >>  8) & 0x0000ff00) | \
                   (((a) >> 24) & 0x000000ff) )

#define SWAP8B(a) ( (((a) & 0x00000000000000FFULL) << 56) | \
                    (((a) & 0x000000000000FF00ULL) << 40) | \
                    (((a) & 0x0000000000FF0000ULL) << 24) | \
                    (((a) & 0x00000000FF000000ULL) <<  8) | \
                    (((a) & 0x000000FF00000000ULL) >>  8) | \
                    (((a) & 0x0000FF0000000000ULL) >> 24) | \
                    (((a) & 0x00FF000000000000ULL) >> 40) | \
                    (((a) & 0xFF00000000000000ULL) >> 56) )

static void swap4b(void *val)
{
    uint32_t *op = (uint32_t*)val;
    *op = SWAP4B(*op);
}

static void swap8b(unsigned long long *val)
{
    uint64_t *op = (uint64_t*)val;
    *op = SWAP8B(*op);
}

int main(int argc, char **argv)
{
    char magic[4];
    volatile uint32_t LE=0x01234567;
    int i, fd, is_little_endian = (*((uint8_t*)(&LE))) == 0x67;

    if (argc != 3) {
        printf("Usage: %s filename num_recs\n", argv[0]);
        exit(1);
    }

    fd = open(argv[1], O_RDWR);
    if (fd == -1) {
        fprintf(stderr, "Error on open file %s (%s)\n",argv[1],strerror(errno));
        exit(1);
    }

    /* read file signature */
    read(fd, magic, 4);
    printf("input file format is CDF %d\n",(int)magic[3]);

    if (magic[3] < 5) {
        int old_nrec, new_nrec;
        read(fd, &old_nrec, 4);
        if (is_little_endian) swap4b(&old_nrec);
        lseek(fd, 4, SEEK_SET);
        new_nrec = atoi(argv[2]);
        printf("change number of records from %d to %d\n",old_nrec,new_nrec);
        if (is_little_endian) swap4b(&new_nrec);
        write(fd, &new_nrec, 4);
    } else {
        long long old_nrec, new_nrec;
        read(fd, &old_nrec, 8);
        if (is_little_endian) swap8b(&old_nrec);
        printf("old numrecs=%d\n",old_nrec);
        lseek(fd, 8, SEEK_SET);
        new_nrec = atoi(argv[2]);
        printf("change number of records from %lld to %lld\n",old_nrec,new_nrec);
        if (is_little_endian) swap8b(&new_nrec);
        write(fd, &new_nrec, 8);
    }
    close(fd);
    return 0;
}
