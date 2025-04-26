#ifndef CHKIO_HINCLUDED
#define CHKIO_HINCLUDED

#include <stdio.h>
#include <stdlib.h>

#ifndef N_DIM
# define N_DIM 3
#endif

// for CHKPART
#include "pkd.h"

typedef int ptcl_idx_t;

#ifdef SSIO_USE_MPI

#include "mpi.h"

/*
 *  can write as many particles as we want in parallel when doing ssic etc
 *  - will need more than this amount of memory to run the simulation anyway.
 *  when writing out during the simulation however, still need to buffer this
 */
#define MAX_WRITE 50000

/*
 *  Need to limit amount of particles that are read during simulations so
 *  it doesn't use too much memory on top of the already allocated space
 */
#define MAX_READ 50000

typedef struct chkio {
    MPI_File mfile;
    int fmode;

    size_t max_buf_particles;

    size_t total_to_buffer;
    size_t total_read;
    size_t extra;

    size_t particles_written;
    size_t particles_read;

    char * file_buf;
    char * cur_buf_ptr;
} CHKIO;

void chkioInitialise (void);
void chkioSetGroupFile (void);
void chkioSetSingleFile(void);

int chkioSetPos(CHKIO *chkio, const MPI_Offset pos);

#else

typedef struct chkio {
	FILE *fp;
} CHKIO;

int chkioSetPos(CHKIO *chkio, const size_t pos);
#endif

#define CHKIO_READ 0
#define CHKIO_WRITE 1

int chkioOpen(const char *filename, CHKIO *chkio, uint mode);
int chkioData(CHKIO *chkio, CHKPART *data);
int chkioClose(CHKIO *chkio);

#endif
