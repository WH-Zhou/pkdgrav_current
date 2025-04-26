/*
 ** ss2vtk.c -- Yun Zhang 16-8-16
 ** ========
 ** Converts *.ss binary data to .vtk format, for visualization in paraview.
 ** unit: m,kg,s
 **
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <vector.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ssio.h"
#include "ss.h"

#define VTK_EXT ".vtk"
#define DEMDIAG_EXT ".demdiag"
#define DEMVTK_EXT "pp.vtk"
#define WALLVTK_EXT "pw.vtk"
#define RED_EXT ".r"

typedef struct force_data {
        int nc;  /* the coordination number of each particle */
        double force[3]; /* the total force acting on each particle */
        } FORCEDATA;

typedef struct contact_data { /* contact information for each particle */
        int *nc; /* the neighborhood list */
        double *forceN; /* the magnitude of the normal force */
        double *forceT; /* the magnitude of the tangetial force */
        } CONTACTDATA;

typedef struct wall_data {
        int bwcontact; /* ==1: if the current particle is in contact with a wall */
        int *nc;
        double *posX,*posY,*posZ; /* the contact position on the wall (allow three contact) */
        double *forceN;  /* the magnitude of the normal force */
        double *forceT;  /* the magnitude of the tangetial force */
        } WALLDATA;

static void
Sphere_paraview(char *outfile, SSDATA *data, int N_n, int bContactP, int bContactPW, int bReduced, FORCEDATA *ff) /* for particle visualization */
{
        FILE *fpo;
        int i;
        assert((fpo = fopen(outfile,"w")) != NULL);
        fprintf(fpo,"# vtk DataFile Version 3.0\n");
        fprintf(fpo,"Particle information visualization\n");
        fprintf(fpo,"ASCII\n");
        fprintf(fpo,"DATASET UNSTRUCTURED_GRID\n");
        fprintf(fpo,"POINTS %d double\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%.16e %.16e %.16e\n", data[i].pos[X]*L_SCALE,data[i].pos[Y]*L_SCALE, data[i].pos[Z]*L_SCALE);  /* the position data*/
        fprintf(fpo,"CELLS %d %d\n",N_n,N_n*2);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"1 %d\n",i);
        fprintf(fpo,"CELL_TYPES %d\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"1\n");
        fprintf(fpo,"POINT_DATA %d\n",N_n);
        if (bReduced == 0) {
                fprintf(fpo,"VECTORS velocity double\n");
                for (i=0;i<N_n;i++)
                        fprintf(fpo,"%.16e %.16e %.16e\n", data[i].vel[X]*V_SCALE,data[i].vel[Y]*V_SCALE, data[i].vel[Z]*V_SCALE);
                fprintf(fpo,"VECTORS spin double\n");
                for (i=0;i<N_n;i++)
                        fprintf(fpo,"%.16e %.16e %.16e\n",data[i].spin[X]/T_SCALE,data[i].spin[Y]/T_SCALE,data[i].spin[Z]/T_SCALE);
                if (bContactP || bContactPW) {
                        fprintf(fpo,"VECTORS contact_force double\n");
                        for (i=0;i<N_n;i++)
                                fprintf(fpo,"%.16e %.16e %.16e\n",ff[i].force[X],ff[i].force[Y],ff[i].force[Z]);
                        }
                if (bContactP || bContactPW)
                        fprintf(fpo,"FIELD FieldData  6\n");
                else
                        fprintf(fpo,"FIELD FieldData  5\n");
                }
        else
                fprintf(fpo,"FIELD FieldData  5\n");
        fprintf(fpo,"mass 1 %d double\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%.16e\n",data[i].mass*M_SCALE);
        fprintf(fpo,"radius 1 %d double\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%.16e\n",data[i].radius*L_SCALE);
        fprintf(fpo,"org_idx 1 %d int\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%d\n",data[i].org_idx);
        fprintf(fpo,"Number 1 %d int\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%d\n",i);
        fprintf(fpo,"color 1 %d int\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%d\n",data[i].color);
        if (bReduced == 0) {
                if (bContactP || bContactPW) {
                        fprintf(fpo,"CN 1 %d int\n",N_n);
                        for (i=0;i<N_n;i++)
                                fprintf(fpo,"%d\n",ff[i].nc);
                        }
                }
        fclose(fpo);
        }

static void
Chain_paraview(char *outfile, CONTACTDATA *cc, SSDATA *data, int N_n, int Par_Contact, int PARCONTACT) /* for force chain between particles */
{
        FILE *fpo;
        int i,j;

        assert((fpo = fopen(outfile,"w")) != NULL);
        fprintf(fpo,"# vtk DataFile Version 3.0\n");
        fprintf(fpo,"Force chain distribution\n");
        fprintf(fpo,"ASCII\n");
        fprintf(fpo,"DATASET POLYDATA\n");
        fprintf(fpo,"POINTS %d double\n",N_n);
        for (i=0;i<N_n;i++)
                fprintf(fpo,"%.16e %.16e %.16e\n", data[i].pos[0]*L_SCALE,data[i].pos[1]*L_SCALE, data[i].pos[2]*L_SCALE);  /* the position data*/
        fprintf(fpo,"LINES %d %d\n",Par_Contact,Par_Contact*3);
        for (i=0;i<N_n;i++)
                for (j=0;j<PARCONTACT;j++)
                        if (cc[i].nc[j] >= 0)
                                fprintf(fpo,"2 %d %d\n",i,cc[i].nc[j]); /* connect line */
        fprintf(fpo,"CELL_DATA  %d\n",Par_Contact);
        fprintf(fpo,"SCALARS fn double 1\n");
        fprintf(fpo,"LOOKUP_TABLE default\n");
        for (i=0;i<N_n;i++)
                for (j=0;j<PARCONTACT;j++)
                        if (cc[i].nc[j] >= 0)
                                fprintf(fpo,"%.16e\n",cc[i].forceN[j]);
        fprintf(fpo,"SCALARS ft double 1\n");
        fprintf(fpo,"LOOKUP_TABLE default\n");
        for (i=0;i<N_n;i++)
                for (j=0;j<PARCONTACT;j++)
                        if (cc[i].nc[j] >= 0)
                                fprintf(fpo,"%.16e\n",cc[i].forceT[j]);
        fclose(fpo);
        }

static void
Wall_contact_paraview(char *outfile, WALLDATA *ww, SSDATA *data, int N_n, int Wall_Contact, int WALLCONTACT) /* for force chain between particles and walls */
{
        FILE *fpo;
        int Particleonwall,i,j,k,l;

        Particleonwall = 0;
        for (i=0;i<N_n;i++)
                if (ww[i].bwcontact > 0)
                        Particleonwall++;
        assert((fpo = fopen(outfile,"w")) != NULL);
        fprintf(fpo,"# vtk DataFile Version 3.0\n");
        fprintf(fpo,"Force Chain on wall\n");
        fprintf(fpo,"ASCII\n");
        fprintf(fpo,"DATASET POLYDATA\n");
        fprintf(fpo,"POINTS %d double\n",Particleonwall + Wall_Contact);
        for (i=0;i<N_n;i++)
                if (ww[i].bwcontact > 0)
                        fprintf(fpo,"%.16e %.16e %.16e\n", data[i].pos[0]*L_SCALE,data[i].pos[1]*L_SCALE, data[i].pos[2]*L_SCALE);  /* the position data*/
        for (i=0;i<N_n;i++)
                if (ww[i].bwcontact > 0)
                        for (j=0;j<WALLCONTACT;j++)
                                if (ww[i].nc[j] >= 0)
                                        fprintf(fpo,"%.16e %.16e %.16e\n", ww[i].posX[j]*L_SCALE,ww[i].posY[j]*L_SCALE,ww[i].posZ[j]*L_SCALE);
        fprintf(fpo,"LINES %d %d\n",Wall_Contact,Wall_Contact*3);
        k = 0;
        l = -1;
        for (i=0;i<N_n;i++)
                if (ww[i].bwcontact > 0) {
                        l++;
                        for (j=0;j<WALLCONTACT;j++)
                                if (ww[i].nc[j] >= 0) {
                                        fprintf(fpo,"2 %d %d\n",l,Particleonwall+k); /* connect line */
                                        k++;
                                        }
                        }
        fprintf(fpo,"CELL_DATA  %d\n",Wall_Contact);
        fprintf(fpo,"SCALARS fn double 1\n");
        fprintf(fpo,"LOOKUP_TABLE default\n");
        for (i=0;i<N_n;i++)
                if (ww[i].bwcontact > 0)
                        for (j=0;j<WALLCONTACT;j++)
                                if (ww[i].nc[j] >= 0)
                                        fprintf(fpo,"%.16e\n",ww[i].forceN[j]);
        fprintf(fpo,"SCALARS ft double 1\n");
        fprintf(fpo,"LOOKUP_TABLE default\n");
        for (i=0;i<N_n;i++)
            if (ww[i].bwcontact > 0)
                    for (j=0;j<WALLCONTACT;j++)
                            if (ww[i].nc[j] >= 0)
                                    fprintf(fpo,"%.16e\n",ww[i].forceT[j*2]);
        fclose(fpo);
        }

static int
convert(char *infile, int bnParticle, int bContactP, int bContactPW, int bQuiet)
{
        SSIO ssio;
        SSHEAD h;
        char outfile[MAXPATHLEN];
        SSDATA *data;
        FORCEDATA *ff = NULL;
        int i,j,bReduced = 0;

        /* Input process */
        if (ssioOpen(infile,&ssio,SSIO_READ)) {
                (void) fprintf(stderr,"Unable to open \"%s\"\n",infile);
                return 1;
                }
        if (ssioHead(&ssio,&h) || h.n_data < 0) {
                (void) fprintf(stderr,"Corrupt header\n");
                (void) ssioClose(&ssio);
                return 1;
                }
        if (h.n_data == 0) {
                (void) fprintf(stderr,"No data found!");
                (void) ssioClose(&ssio);
                return 1;
                }
        if (!bQuiet) printf("%s: Number of particles = %i (time %g)\n",infile,h.n_data,h.time);
        data = (SSDATA *) malloc((size_t) h.n_data*sizeof(SSDATA));
        assert(data != NULL);
        switch(h.iMagicNumber) {
        case SSIO_MAGIC_STANDARD: {
                for (i=0;i<h.n_data;i++)
                        if (ssioData(&ssio,&data[i])) {
                                (void) fprintf(stderr,"Corrupt data\n");
                                (void) ssioClose(&ssio);
                                return 1;
                                }
                break;
                }
        case SSIO_MAGIC_REDUCED: {  /* Reduced ss format */
                SSRDATA data_r;
                bReduced = 1;
                for (i=0;i<h.n_data;i++) {
                        if (ssioDataReduced(&ssio,&data_r)) {
                                (void) fprintf(stderr,"Corrupt data\n");
                                (void) ssioClose(&ssio);
                                return 1;
                                }
                        data[i].org_idx = data_r.iOrgIdx;
                        data[i].mass = data_r.fMass;
                        data[i].radius = data_r.fRadius;
                        data[i].pos[X] = data_r.vPos[X];
                        data[i].pos[Y] = data_r.vPos[Y];
                        data[i].pos[Z] = data_r.vPos[Z];
                        data[i].vel[X] = 0.0;
                        data[i].vel[Y] = 0.0;
                        data[i].vel[Z] = 0.0;
                        data[i].spin[X] = 0.0;
                        data[i].spin[Y] = 0.0;
                        data[i].spin[Z] = 0.0;
                        data[i].color = data_r.iColor;
                        }
                break;
                }
        default:
                fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
                return 1;
                }
        (void) ssioClose(&ssio);

        /* Particle visualization */
        if (ssioNewExt(infile,SS_EXT,outfile,VTK_EXT)) {
                fprintf(stderr,"Unable to generate output filename.\n");
                return 1;
                }

        /* Contact force chain calculation and visualization */
        if ( !bReduced ) {
                if ( bContactP || bContactPW ) {
                        char demfile[MAXPATHLEN],demvtk[MAXPATHLEN];
                        double Fn[3],Ft[3],temp,tempv[3],scale;
                        CONTACTDATA *cc; /* particle contact */
                        WALLDATA *ww = NULL; /* wall contact */
                        int c,pc,pw,Par_Contact,Wall_Contact;
                        FILE *fpi;

                        pc = 0;
                        pw = 0;
                        if ( ssioNewExt(infile,SS_EXT,demfile,DEMDIAG_EXT) ) {
                                fprintf(stderr,"Unable to generate demfile filename.\n");
                                return 1;
                                }
                        assert((fpi = fopen(demfile,"r")) != NULL);
                        fscanf(fpi,"%i",&i); /* The number of particle */
                        fscanf(fpi,"%i",&pc); /* The number of particle-particle contact */
                        if (bContactPW) fscanf(fpi,"%i",&pw); /* The number of particle-wall contact */
                        cc = (CONTACTDATA *) malloc((size_t) h.n_data*sizeof(CONTACTDATA));
                        ff = (FORCEDATA *) malloc((size_t) h.n_data*sizeof(FORCEDATA));
                        for (i=0;i<h.n_data;i++) {
                                ff[i].nc = 0;
                                ff[i].force[X] = 0.0;
                                ff[i].force[Y] = 0.0;
                                ff[i].force[Z] = 0.0;
                                cc[i].nc = (int *) malloc((size_t) pc*sizeof(int));
                                cc[i].forceN = (double *) malloc((size_t) pc*sizeof(double));
                                cc[i].forceT = (double *) malloc((size_t) pc*sizeof(double));
                                for (j=0;j<pc;j++)
                                        cc[i].nc[j] = -1;
                                }
                        if (bContactPW) {
							ww = (WALLDATA *) malloc((size_t) h.n_data*sizeof(WALLDATA));
                                for (i=0;i<h.n_data;i++) {
                                        ww[i].bwcontact = -1;
                                        ww[i].nc = (int *) malloc((size_t) pw*sizeof(int));
                                        ww[i].forceN = (double *) malloc((size_t) pw*sizeof(double));
                                        ww[i].forceT = (double *) malloc((size_t) pw*sizeof(double));
                                        ww[i].posX = (double *) malloc((size_t) pw*sizeof(double));
                                        ww[i].posY = (double *) malloc((size_t) pw*sizeof(double));
                                        ww[i].posZ = (double *) malloc((size_t) pw*sizeof(double));
                                        for (j=0;j<pw;j++)
                                                ww[i].nc[j] = -1;
                                        }
                                }
                        i = j = 0;  /* i is the current particle, j is the order of neighborhood */
                        Par_Contact = 0;
                        Wall_Contact = 0;
                        scale =  M_SCALE*V_SCALE/T_SCALE;
                        while ( fscanf(fpi,"%i",&c) > 0 ) {
                                fscanf(fpi,"%lf%lf%lf%lf%lf%lf",&Fn[X],&Fn[Y],&Fn[Z],&Ft[X],&Ft[Y],&Ft[Z]);
                                if (c >= 0 && j < pc) { /* contact between particle and particle */
                                        if ( !pw && !bContactPW && i >= h.n_data) {
                                                fprintf(stderr,"simulation includes walls: need to use -w mode.\n");
                                                return 1;
                                                }
                                        cc[i].nc[j] = c;
                                        cc[i].forceN[j] = MAG(Fn)*scale;
                                        cc[i].forceT[j] = MAG(Ft)*scale;
                                        ADD_VEC(Fn,Ft,tempv);
                                        ff[i].nc += 1;
                                        ff[i].force[X] += tempv[X]*scale;
                                        ff[i].force[Y] += tempv[Y]*scale;
                                        ff[i].force[Z] += tempv[Z]*scale;
                                        ff[cc[i].nc[j]].nc += 1;
                                        ff[cc[i].nc[j]].force[X] -= tempv[X]*scale;
                                        ff[cc[i].nc[j]].force[Y] -= tempv[Y]*scale;
                                        ff[cc[i].nc[j]].force[Z] -= tempv[Z]*scale;
                                        Par_Contact++;
                                        }
                                else if ( c >= 0 && j >= pc ) { /* contact between particle and wall */
                                        ww[i].nc[j-pc] = c;
                                        ww[i].bwcontact = 1;
                                        temp = - data[i].radius/MAG(Fn);
                                        ww[i].posX[j-pc] = data[i].pos[X] + temp*Fn[X];
                                        ww[i].posY[j-pc] = data[i].pos[Y] + temp*Fn[Y];
                                        ww[i].posZ[j-pc] = data[i].pos[Z] + temp*Fn[Z];
                                        ww[i].forceN[j-pc] = MAG(Fn)*scale;
                                        ww[i].forceT[j-pc] = MAG(Ft)*scale;
                                        ADD_VEC(Fn,Ft,tempv);
                                        ff[i].nc += 1;
                                        ff[i].force[X] += tempv[X]*scale;
                                        ff[i].force[Y] += tempv[Y]*scale;
                                        ff[i].force[Z] += tempv[Z]*scale;
                                        Wall_Contact++;
                                        }
                                j++;
                                if (j >= pc + pw) {
                                        j = 0;
                                        i++;
                                        }
                                }
                        fclose(fpi);
                        if ( ssioNewExt(infile,SS_EXT,demvtk,DEMVTK_EXT) ) {
                                fprintf(stderr,"Unable to generate demvtk filename.\n");
                                return 1;
                                }
                        if ( Par_Contact > 0)
                                Chain_paraview(demvtk, cc, data, h.n_data, Par_Contact, pc); /* Visualization */
                        else
                                if (!bQuiet) printf("No particle contacts with others.\n");
                        free(cc);
                        if ( bContactPW ) {
                                if ( ssioNewExt(infile,SS_EXT,demfile,WALLVTK_EXT) ) {
                                fprintf(stderr,"Unable to generate walldemfile filename.\n");
                                return 1;
                                        }
                                if ( Wall_Contact > 0)
                                        Wall_contact_paraview(demfile, ww, data, h.n_data, Wall_Contact,pw); /* Visualization */
                                else
                                        if (!bQuiet) printf("No particle contacts with walls.\n");
                                }
                        }
                }
        else {
                if ( bContactP || bContactPW )
                        if (!bQuiet) printf("Reduced ss format is not supported for force analyse.\n");
                }

        if (bnParticle) Sphere_paraview(outfile, data, h.n_data, bContactP, bContactPW, bReduced, ff);

        return 0;
        }

static void usage(const char *achProgName)
{
        fprintf(stderr,"Usage: %s [ -n ] [ -c ] [ -w ] [ -q ] file [ file ... ]\n"
                "Where: -n = don't convert .ss data to .vtk format\n"
                "       -c = analyze and output force chain distribution between particles (for simulation without walls) \n"
                "       -w = analyze and output force chain distribution between particles and walls (for simulation with walls)\n"
                "       -q = activate quiet mode\n",achProgName);
        exit(0);
        }

int
main(int argc,char *argv[])
{
        int i;
        int bnParticle = 1; /* default = output particle .vtk data for visualization in paraview*/
        int bQuiet = 0; /* default = verbose */
        int bContactP = 0;    /* default = don't analyze force chain distribution between particles */
        int bContactPW = 0;   /* default = don't analyze force chain distribution between particles and walls */
        int flag = 0;  /* default = get correct output */

        while ( (i = getopt(argc, argv, "n::c::w::q::")) != EOF ) {
                switch (i) {
                case 'n':
                        bnParticle = 0;
                        break;
                case 'c':
                        bContactP = 1;
                        break;
                case 'w':
                        bContactPW = 1;
                        break;
                case 'q':
                        bQuiet = 1;
                        break;
                default:
                        usage(argv[0]);
                        }
                }
        if ( optind >= argc )
                usage(argv[0]);
		/*(void) printf("The number of input file = %i\n",argc-1);*/ /*DEBUG*/
        for ( i=optind;i<argc;++i )
               flag = convert(argv[i],bnParticle,bContactP, bContactPW, bQuiet);
        if (!bQuiet && !flag) printf("Open paraview to check the results and have fun!\n");
        return 0;
	}

/* ss2vtk.c */
