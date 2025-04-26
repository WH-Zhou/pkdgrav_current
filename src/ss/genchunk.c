
// if MPI SSIO is not used, this file is never used either
#ifdef SSIO_USE_MPI

#include <ss.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>		/* for getpid(), and getopt() if needed */
#include <sys/types.h>	/* ditto (for IRIX) */
#include <boolean.h>
#include <rdpar.h>
#include <vector.h>

#ifdef sparc
#ifdef sun
#undef sun /* sheesh */
#endif
#endif

double binm1;
double binm2;
#define aul 149597870700.0
#define gconst 6.673e-11
double separation;
double eccentricity;

/*#define SSTEST*/

/*#define TILT_RING*/ /*DEBUG used to tilt a portion of the disk at the start*/

/* Definitions */

typedef struct {

	/* Following data read in from parameter file... */

	double central_mass;
	BOOLEAN heliocentric;
	
	int binaryop;
	double binm1, binm2, separation, eccentricity;

	int n,dst_fnc;
	double total_mass,density,radius,scaling,r_inner,r_outer,surf_den_exp;
	double ecc_dsp,inc_dsp,ecc_max,inc_max;

	double seed_mass,seed_density,seed_radius,seed_scaling;
	double seed_sma,seed_ecc,seed_inc,seed_gap_scale;

	int n_planets;
	BOOLEAN tilt;
	char planet_data[MAXPATHLEN];

	double time;
	BOOLEAN adjust_com,softening;
	char output_file[MAXPATHLEN];

	/* Following derived from supplied parameters... */

	int n_data;
	double mass,red_hill,esc_vel,seed_half_gap;

	} PARAMS;

// prototypes, defined in ssic.c
void
gen_planetesimals(PARAMS *p,SSDATA *data);
void
add_mom(SSDATA *data,VECTOR mom_pos,VECTOR mom_vel);

/*
 *  Generate a fraction of the work to be done
 */
void
generate_chunk
(PARAMS params,
 SSDATA *orig_data)
{
	SSDATA *data = NULL;

    // same data
    PARAMS * p = &params;

    int pltsml_orig = params.n_data;
    int pltsml_left = pltsml_orig;

	SSIO ssio;

    ssioSetGroupFile();

    /* file output */
	if (ssioOpen(params.output_file,&ssio,SSIO_WRITE))
    {
        fprintf(stderr,"Unable to open \"%s\"\n", params.output_file);
        exit(1);
    }

    // written the stars/header
    int initd = 0;

    // index of particle - do not reset across iterations of while loop!
    int idx = 0;

    // serial
    int rank = 0;

    // rank 0 writes header
    int comm_sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    if (!rank)
    {
        fprintf(stdout, "%d initial particles\n", pltsml_orig);

        // rank 0 does the header bits
        SSHEAD head;
        head.time = params.time;
        head.n_data = params.n_data;
        head.iMagicNumber = SSIO_MAGIC_STANDARD;
        ssioHead(&ssio,&head);
    }
    else
    {
        // make sure not to add more stars in children
        p->binaryop = 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    pltsml_left = pltsml_orig / comm_sz;
	// add any extra
	if (!rank)
	{
		int extra = pltsml_orig%pltsml_left;
		pltsml_left += extra;
	}

    // TODO this could probably be more robust
    p->n = (pltsml_left >= MAX_OPEN) ? MAX_OPEN : pltsml_left;
    idx = pltsml_left*rank;

    do
    {
        p->n = (pltsml_left >= MAX_OPEN) ? MAX_OPEN : pltsml_left;
        pltsml_left -= p->n;

        // allocate memory as needed
        if (!data || p->n != MAX_OPEN)
        {
            if (!(data = (SSDATA *) realloc(data,
                    (p->n+1) * sizeof(SSDATA))))
            {
                fprintf(stderr,"Unable to allocate data memory\n");
                exit(1);
            }
        }

        /***********************************/

        SSDATA *sun,*planetesimals;
        double process_mass = 0;
        int i;

        planetesimals = data;

        if (!initd && !rank)
        {
            /* Get handy pointers */

            if (p->heliocentric)
            {
                sun = NULL;
            }
            else
            {
                sun = data;
                ZERO_VEC(sun->pos);
                ZERO_VEC(sun->vel);
                ZERO_VEC(sun->spin);
                sun->mass = 1;
                sun->radius = R_SUN/L_SCALE;
                sun->color = SUN;

                planetesimals = sun + 1;
            }
        }

        gen_planetesimals(p, planetesimals);

// not done in this run
#if 0
        if (p->tilt)
        {
            tilt(p,planets);
    #ifdef SSTEST
            {
                VECTOR z,zp;

                assert(p->seed_mass == 0.0);
                calc_ang_mom(planetesimals,p->n,z);
                calc_ang_mom(planets,p->n_planets,zp);

                NORM_VEC(z,MAG(z));
                NORM_VEC(zp,MAG(zp));

                 printf("net planetesimal ang mom direction = (%e,%e,%e)\n",
                              z[X],z[Y],z[Z]);

                 printf("-- will assume (0,0,1)\n");

                SET_VEC(z,0,0,1);

                SUB_VEC(z,zp,z);

                 printf("invar plane: z - zp = (%e,%e,%e)\n",z[X],z[Y],z[Z]);
            }
    #endif
        }
#endif

        /* Get planetesimal moments and estimate of mass centering error */

        VECTOR com_pos,com_vel;
        double total_mass = 0.0;
        ZERO_VEC(com_pos);
        ZERO_VEC(com_vel);

        // should always be >0
        if (p->n)
        {
            // have to reduce the total mass across all processes
            //get_com(planetesimals,p->n,com_pos,com_vel,&process_mass);

            /* get_com */
            for (i = 0; i < p->n; i++)
            {
                add_mom(&data[i], com_pos, com_vel);
                process_mass += data[i].mass;
            }

            MPI_Allreduce(&process_mass,
                          &total_mass,
                          1,
                          MPI_DOUBLE,
                          MPI_SUM,
                          MPI_COMM_WORLD);

            NORM_VEC(com_pos, total_mass);
            NORM_VEC(com_vel, total_mass);
            /* get_com */

            printf("Planetesimal c-o-m (pos,vel) offset mag = (%e,%e)\n",
                   MAG(com_pos), MAG(com_vel));
            SCALE_VEC(com_pos, total_mass);
            SCALE_VEC(com_vel, total_mass);
        }

// not done in this run
#if 0
        /* Now shift everything so centre-of-mass is at origin */
        SSDATA *planets,*seed,*ptr;

        if (p->adjust_com) {
            assert(p->seed_mass == 0.0);
            if (!initd && !p->heliocentric)
                total_mass += sun->mass;

            for (i=0;i<p->n_planets;i++)
            {
                ptr = &planets[i];
                add_mom(ptr,com_pos,com_vel);
                total_mass += ptr->mass;
            }

            NORM_VEC(com_pos,total_mass);
            NORM_VEC(com_vel,total_mass);

            sub_com(data,p->n_data,com_pos,com_vel);
        }
#endif

        // mark as done so these don't get done again
        initd = 1;
        p->binaryop = 0;

        /***********************************/

        // write out this chunk on the end of the file
        for (i = 0; i < p->n; i++) {
            // update idx at the same time
            data[i].org_idx = idx++;
            ssioData(&ssio, &data[i]);
        }
    }
    while (pltsml_left);

    // free before waiting to close
    free(data);

	ssioClose(&ssio);
    /* file output */

    MPI_Finalize();
}

#endif

