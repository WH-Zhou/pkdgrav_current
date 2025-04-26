/*
 ** carvepoly.c -- Yun Zhang 11-8-16. Modified: Joe DeMartini 7-28-20
 ** ========
 ** Carve input particle data to form a aggregate according to the given Polyhedron.
 ** The result can be visualized in Paraview.
 ** unit: m,kg,s
 **
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "vector.h"
#include "ssio.h"
#include "ss.h"

#define R_scale 1000 /* Normally, the unit in the obj file is km, while in the ssdata file is m */
#define OUTPUTFILE "carvepoly.ss"
#define PARAS "particles.vtk" /* check for paraview: sphere and face */
#define PARAF "faces.vtk"

typedef struct particle_data {  /* particle data */
        SSDATA data;
	int flag; /* the number of intersections between the line O-particle and the surface; even denotes inside the surface, odd denotes outside */
} PARTICLEDATA;

typedef struct v_data {   /* vertex data: rectangular coordinate system & spherical coordinate system */
	double pos[3];
} DATAV;

typedef struct f_data {   /* face data: vertex information & range expressed in spherical coordinate system */
	int v[3];
} DATAF;

typedef struct tetrahedron {
	double ut,vt,wt,Ut,Vt,Wt;  /* the lengths of edges of the tetrahedron */
	double Xt,xt,Yt,yt,Zt,zt;
	double at,bt,ct,dt;  /* the coefficient for calculating the volume of a tetrahedron using Heron-type formula */
} TETRA;

static void face_paraview(const char *outfile, DATAV *datav, DATAF *dataf, int N_f) /* for face paraview */
{
	FILE *fpo;
	int i;
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	fprintf(fpo,"# vtk DataFile Version 3.0\n");
	fprintf(fpo,"Asteroid surface\n");
	fprintf(fpo,"ASCII\n");
	fprintf(fpo,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fpo,"POINTS %d double\n",N_f*3);
	for (i=0;i<N_f;i++)  /* the vertex data */
		fprintf(fpo,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", datav[dataf[i].v[0]].pos[0],
				datav[dataf[i].v[0]].pos[1],datav[dataf[i].v[0]].pos[2],datav[dataf[i].v[1]].pos[0],
				datav[dataf[i].v[1]].pos[1],datav[dataf[i].v[1]].pos[2],datav[dataf[i].v[2]].pos[0],
				datav[dataf[i].v[2]].pos[1],datav[dataf[i].v[2]].pos[2]); 
	fprintf(fpo,"CELLS %d %d\n",N_f,N_f*4);
	for (i=0;i<N_f;i++)  /* the face data */
		fprintf(fpo, "3 %d %d %d\n",i*3,i*3+1,i*3+2 );
	fprintf(fpo,"CELL_TYPES %d\n",N_f);
	for (i=0;i<N_f;i++) 
		fprintf(fpo,"5\n");
	fclose(fpo);
	}

static void Sphere_paraview(const char *outfile, PARTICLEDATA *p, int N_p, int N_n) /* for particle paraview */
{
	FILE *fpo;
	int i;
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	fprintf(fpo,"# vtk DataFile Version 3.0\n");
	fprintf(fpo,"Nbody evolution check\n");
	fprintf(fpo,"ASCII\n");
	fprintf(fpo,"DATASET UNSTRUCTURED_GRID\n");
	if (N_n) { 
		fprintf(fpo,"POINTS %d double\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%.16e %.16e %.16e\n", p[i].data.pos[0],p[i].data.pos[1], p[i].data.pos[2]);  /* the position data*/
		fprintf(fpo,"CELLS %d %d\n",N_n,N_n*2);
		for (i=0;i<N_n;i++) 
			fprintf(fpo,"1 %d\n",i);
		fprintf(fpo,"CELL_TYPES %d\n",N_n);
		for (i=0;i<N_n;i++) 
			fprintf(fpo,"1\n");
		fprintf(fpo,"POINT_DATA %d\n",N_n);
		fprintf(fpo,"VECTORS velocity double\n");
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%.16e %.16e %.16e\n", p[i].data.vel[0],p[i].data.vel[1], p[i].data.vel[2]);
		fprintf(fpo,"VECTORS spin double\n");
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%.16e %.16e %.16e\n",p[i].data.spin[0],p[i].data.spin[1],p[i].data.spin[2]);
		fprintf(fpo,"FIELD FieldData  5\n");
		fprintf(fpo,"mass 1 %d double\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%.16e\n",p[i].data.mass);
		fprintf(fpo,"radius 1 %d double\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%.16e\n",p[i].data.radius);
		fprintf(fpo,"org_idx 1 %d int\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%d\n",p[i].data.org_idx);
		fprintf(fpo,"flag 1 %d int\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%d\n",p[i].flag);
		fprintf(fpo,"color 1 %d int\n",N_n);
                for (i=0;i<N_p;i++)
                        if(p[i].flag%2 == 0)
                                fprintf(fpo,"%d\n",p[i].data.color);
		}
	else{
                fprintf(fpo,"POINTS %d double\n",N_p);
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%.16e %.16e %.16e\n", p[i].data.pos[0],p[i].data.pos[1], p[i].data.pos[2]);  /* the position data*/
                fprintf(fpo,"CELLS %d %d\n",N_p,N_p*2);
                for (i=0;i<N_p;i++)
			fprintf(fpo,"1 %d\n",i);
                fprintf(fpo,"CELL_TYPES %d\n",N_p);
                for (i=0;i<N_p;i++)
			fprintf(fpo,"1\n");
                fprintf(fpo,"POINT_DATA %d\n",N_p);
		fprintf(fpo,"VECTORS velocity double\n");
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%.16e %.16e %.16e\n", p[i].data.vel[0],p[i].data.vel[1], p[i].data.vel[2]);
		fprintf(fpo,"VECTORS spin double\n");
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%.16e %.16e %.16e\n",p[i].data.spin[0],p[i].data.spin[1],p[i].data.spin[2]);
		fprintf(fpo,"FIELD FieldData  4\n");
                fprintf(fpo,"mass 1 %d double\n",N_p);
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%.16e\n",p[i].data.mass);
                fprintf(fpo,"radius 1 %d double\n",N_p);
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%.16e\n",p[i].data.radius);
                fprintf(fpo,"org_idx 1 %d int\n",N_p);
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%d\n",p[i].data.org_idx);
                fprintf(fpo,"color 1 %d int\n",N_p);
                for (i=0;i<N_p;i++)
                        fprintf(fpo,"%d\n",p[i].data.color);
		}
	fclose(fpo);
        }


static int
writess(PARTICLEDATA *p, int N_p_new, int N_p, int bVerbose)
{
        SSIO ssio;
        SSHEAD head;
        SSDATA data;
        int i;

        if ( ssioOpen(OUTPUTFILE,&ssio,SSIO_WRITE) ) {
                fprintf(stderr,"Unable to open \"%s\".\n",OUTPUTFILE);
                return 1;
                }
        head.n_data = N_p_new;
        head.time = 0.0;
        head.iMagicNumber = SSIO_MAGIC_STANDARD;
        if (ssioHead(&ssio,&head)) {
                fprintf(stderr,"Error writing header.\n");
                ssioClose(&ssio);
                }
        for (i=0;i<N_p;i++)
                if(p[i].flag%2 == 0) {
                        data.org_idx = p[i].data.org_idx;
                        data.mass = p[i].data.mass / M_SCALE;
                        data.radius = p[i].data.radius / L_SCALE;
                        data.pos[X] = p[i].data.pos[X] / L_SCALE;
                        data.pos[Y] = p[i].data.pos[Y] / L_SCALE;
                        data.pos[Z] = p[i].data.pos[Z] / L_SCALE;
                        data.vel[X] = p[i].data.vel[X] / V_SCALE;
                        data.vel[Y] = p[i].data.vel[Y] / V_SCALE;
                        data.vel[Z] = p[i].data.vel[Z] / V_SCALE;
                        data.spin[X] = p[i].data.spin[X] * T_SCALE;
                        data.spin[Y] = p[i].data.spin[Y] * T_SCALE;
                        data.spin[Z] = p[i].data.spin[Z] * T_SCALE;
                        data.color = p[i].data.color;
                        if (ssioData(&ssio,&data)) {
                                fprintf(stderr,"Error writing data.\n");
                                ssioClose(&ssio);
                                }
                        }
        if (bVerbose) printf("Output file %s: n_data = %i\n",OUTPUTFILE,head.n_data);

        return 0;
}

static void
process_poly(char *ssfile, char *polyfile, int bVtk, int bScale, int bRetain, int bAid, int bVerbose, int bAdel)
{
        SSIO ssio;
        SSHEAD h;
        FILE *fpi;
        PARTICLEDATA *p;
        DATAV *datav;
        DATAF *dataf;
        double temp1,tempv[3];
        double omf,omf_total,n_r1[3],n_r2[3],n_r3[3];
        double minimum_v,maximum_v; /* the radius of the inscribed sphere and cutting sphere of the surface */
        double VolPar,VolPoly;  /* volume*/
        int i,j,flag,N_p,N_v,N_f,N_p_new,DelIdx,delflag,aggstart;
        char ch1,ch2;

        /* read .ss file */
        if (ssioOpen(ssfile,&ssio,SSIO_READ)) {
                (void) fprintf(stderr,"Unable to open \"%s\"\n",ssfile);
                return;
                }
        if (ssioHead(&ssio,&h) || h.n_data < 0) {
                (void) fprintf(stderr,"Corrupt header\n");
                (void) ssioClose(&ssio);
                return;
                }
        if (h.n_data == 0) {
                (void) fprintf(stderr,"No data found!");
                (void) ssioClose(&ssio);
                return;
                }
        p = (PARTICLEDATA *) malloc((size_t) h.n_data*sizeof(PARTICLEDATA));
        assert(p != NULL);
        switch(h.iMagicNumber) {
        case SSIO_MAGIC_STANDARD: {
                for (i=0;i<h.n_data;i++)
                        if (ssioData(&ssio,&p[i].data)) {
                                (void) fprintf(stderr,"Corrupt data\n");
                                (void) ssioClose(&ssio);
                                return;
                                }
                break;
                }
        case SSIO_MAGIC_REDUCED: {  /* Reduced ss format */
                SSRDATA data_r;
                for (i=0;i<h.n_data;i++) {
                        if (ssioDataReduced(&ssio,&data_r)) {
                                (void) fprintf(stderr,"Corrupt data\n");
                                (void) ssioClose(&ssio);
                                return;
                                }
                        p[i].data.org_idx = data_r.iOrgIdx;
                        p[i].data.mass = data_r.fMass;
                        p[i].data.radius = data_r.fRadius;
                        p[i].data.pos[X] = data_r.vPos[X];
                        p[i].data.pos[Y] = data_r.vPos[Y];
                        p[i].data.pos[Z] = data_r.vPos[Z];
                        p[i].data.vel[X] = 0.0;
                        p[i].data.vel[Y] = 0.0;
                        p[i].data.vel[Z] = 0.0;
                        p[i].data.spin[X] = 0.0;
                        p[i].data.spin[Y] = 0.0;
                        p[i].data.spin[Z] = 0.0;
                        p[i].data.color = data_r.iColor;
                        }
                break;
                }
        default:
                fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
                return;
                }
        (void) ssioClose(&ssio);
        for (i=0;i<h.n_data;i++) { /* convert to unit: kg,m,s */
                p[i].data.mass *= M_SCALE;
                p[i].data.radius *= L_SCALE;
                p[i].data.pos[X] *= L_SCALE;
                p[i].data.pos[Y] *= L_SCALE;
                p[i].data.pos[Z] *= L_SCALE;
                if (bRetain) {
                        p[i].data.vel[X] *= V_SCALE;
                        p[i].data.vel[Y] *= V_SCALE;
                        p[i].data.vel[Z] *= V_SCALE;
                        p[i].data.spin[X] /= T_SCALE;
                        p[i].data.spin[Y] /= T_SCALE;
                        p[i].data.spin[Z] /= T_SCALE;
                        }
                else {
                        p[i].data.vel[X] = 0.0;
                        p[i].data.vel[Y] = 0.0;
                        p[i].data.vel[Z] = 0.0;
                        p[i].data.spin[X] = 0.0;
                        p[i].data.spin[Y] = 0.0;
                        p[i].data.spin[Z] = 0.0;
                        }
                }
        N_p = h.n_data;
        if (bVerbose) printf("%s: Number of particles = %i\n",ssfile,N_p);

        /* read polyhedron file */
        if ( !(fpi = fopen(polyfile,"r"))) {  /* polyhedron model, style: .obj */
                fprintf(stderr,"Unable to open \"%s\".\n",polyfile);
                return;
                }
        N_v = 0;
        N_f = 0;
        while((ch1 = (char) fgetc(fpi))!=EOF){ /* count the number of v and f in the polyhedron file */
                if(ch1 == '#'){
					    ch2 = (char) fgetc(fpi);
                        while(ch2 != EOF && ch2!= '\n') ch2 = (char) fgetc(fpi);
                        }
                else{
                        if(ch1 == 'v') N_v++;
                        else if(ch1 == 'f') N_f++;
                        }
                }
        rewind(fpi);
        datav = (DATAV *) malloc((size_t) N_v *sizeof(DATAV));
        dataf = (DATAF *) malloc((size_t) N_f *sizeof(DATAF));
        i=0; j=0;
        while((ch1 = (char) fgetc(fpi))!=EOF){ /* read v and f information */
                if(ch1 == '#'){
					ch2 = (char) fgetc(fpi);
                        while(ch2 != EOF && ch2!= '\n') ch2 = (char) fgetc(fpi);
                        }
                else{
                        if(ch1 == 'v') {
                                fscanf(fpi,"%lf%lf%lf",&datav[i].pos[0],&datav[i].pos[1],&datav[i].pos[2]);
                                if (bScale) {
                                        datav[i].pos[0] *= R_scale;
                                        datav[i].pos[1] *= R_scale;
                                        datav[i].pos[2] *= R_scale;
                                        }
                                i++;
                                }
                        else if(ch1 == 'f') {
                                fscanf(fpi,"%d%d%d",&dataf[j].v[0],&dataf[j].v[1],&dataf[j].v[2]);
                                dataf[j].v[0] --;
                                dataf[j].v[1] --;
                                dataf[j].v[2] --;
                                j++;
                                }
                        }
                }
        fclose(fpi);
        if (bVerbose) printf("%s: v_data = %i, f_data = %i\n",polyfile,N_v,N_f);

        /* carve shape */
        /* present v in spherical coordinate frame */
        minimum_v = 1.0e15;
        maximum_v = 0.0;
        for (i=0;i<N_v;i++) {
                temp1 = sqrt( datav[i].pos[0]*datav[i].pos[0] + datav[i].pos[1]*datav[i].pos[1] + datav[i].pos[2]*datav[i].pos[2] ); /* r */
                if(temp1 < 1.0e-16) printf("Singular point!");
                if(temp1 < minimum_v) minimum_v = temp1;
                if(temp1 > maximum_v) maximum_v = temp1;
                }
        /* check if the given particle is inside of the polyhedron */
	DelIdx = -N_p-1;	/* A number that the org_idx should never be */
	delflag = 0;
	aggstart = -1;
        for (i=0;i<N_p;i++) {
		if (bAdel && p[i].data.org_idx < 0 && p[i].data.org_idx != p[aggstart].data.org_idx)
			aggstart = i;	/* Store where the agg starts so we can flag preceding members */
		if (bAdel) {	/* If an agg member is flagged for deletion, so are the rest */
			if (p[i].data.org_idx == DelIdx) { /* Handles all following agg members */
				p[i].flag = delflag;
				continue; /* Skip further calculations, if flagged for deletion */
				}
			}
                temp1 = sqrt( p[i].data.pos[0]*p[i].data.pos[0] + p[i].data.pos[1]*p[i].data.pos[1] + p[i].data.pos[2]*p[i].data.pos[2]);
                if(temp1 >= maximum_v) {
                        p[i].flag = -1;  /* outside the largest cutting sphere */
			if (bAdel && p[i].data.org_idx < 0) {	/* If agg flagged for deletion, note org index */
				DelIdx = p[i].data.org_idx;
				delflag = -1;
				for (j=aggstart;j<i;j++) /* Handles all preceding agg members */
					p[j].flag = delflag;
				}
			}
                else if(temp1 <= minimum_v)
                        p[i].flag = 0; /* inside the inscribed sphere */
                else {
                        omf_total = 0; /* initialization */
                        for (j=0;j<N_f;j++) {  /* check the number of intersections between the line O-particle and the surface */
                                SUB_VEC(datav[dataf[j].v[0]].pos, p[i].data.pos, tempv);
                                temp1 =  MAG(tempv);
                                n_r1[0] = tempv[0] / temp1;
                                n_r1[1] = tempv[1] / temp1;
                                n_r1[2] = tempv[2] / temp1;
                                SUB_VEC(datav[dataf[j].v[1]].pos, p[i].data.pos, tempv);
                                temp1 =  MAG(tempv);
                                n_r2[0] = tempv[0] / temp1;
                                n_r2[1] = tempv[1] / temp1;
                                n_r2[2] = tempv[2] / temp1;
                                SUB_VEC(datav[dataf[j].v[2]].pos, p[i].data.pos, tempv);
                                temp1 =  MAG(tempv);
                                n_r3[0] = tempv[0] / temp1;
                                n_r3[1] = tempv[1] / temp1;
                                n_r3[2] = tempv[2] / temp1;
                                CROSS(n_r2, n_r3, tempv);
                                omf = 2.0*atan( DOT(n_r1, tempv) / (1 + DOT(n_r1, n_r2) + DOT(n_r2, n_r3) + DOT(n_r3, n_r1) ) );
                                omf_total += omf;
                                }
                        if(fabs(omf_total)<1.0e-5) {
                                p[i].flag = 1;
				if (bAdel && p[i].data.org_idx < 0) {	/* Flag all members with same org_idx for deletion */
					DelIdx = p[i].data.org_idx;
					delflag = 1;
					for (j=aggstart;j<i;j++) /* Handles all preceding agg members */
						p[j].flag = delflag;
					}
				}
                        else
                                p[i].flag = 0;
                        }
                }
        N_p_new = 0;   /* the number of particles in the output .ss file */
        VolPar = 0.0;  /* the total volume of rubble pile */
        for (i=0;i<N_p;i++)
                if(p[i].flag%2 == 0) {
                        if (bAid)
                                p[i].data.org_idx = N_p_new;
                        N_p_new++;
                        VolPar += 4.0/3.0*PI*p[i].data.radius*p[i].data.radius*p[i].data.radius;
                        }
        if (N_p_new == 0) {
                fprintf(stderr,"No particle found inside of the given polyhedron!\n");
                return;
                }

        /* output */
        flag = writess(p, N_p_new, N_p, bVerbose);
        if (flag) return;

        /* packing efficiency estimation*/
        if (bVerbose) {
                TETRA tetra1;  /* tetraheron */
                VolPoly = 0.0; /* the total volume of given polyhedron*/
                for (j=0;j<N_f;j++) {  /* using Heron-type formula to calculate the volume of a tetrahedron */
                        tetra1.Ut = DIST( datav[dataf[j].v[0]].pos, datav[dataf[j].v[1]].pos);
                        tetra1.Vt = DIST( datav[dataf[j].v[1]].pos, datav[dataf[j].v[2]].pos);
                        tetra1.Wt = DIST( datav[dataf[j].v[2]].pos, datav[dataf[j].v[0]].pos);
                        tetra1.ut = MAG( datav[dataf[j].v[2]].pos);
                        tetra1.vt = MAG( datav[dataf[j].v[0]].pos);
                        tetra1.wt = MAG( datav[dataf[j].v[1]].pos);
                        tetra1.Xt = (tetra1.wt - tetra1.Ut + tetra1.vt) * (tetra1.Ut + tetra1.vt + tetra1.wt);
                        tetra1.xt = (tetra1.Ut - tetra1.vt + tetra1.wt) * (tetra1.vt - tetra1.wt + tetra1.Ut);
                        tetra1.Yt = (tetra1.ut - tetra1.Vt + tetra1.wt) * (tetra1.Vt + tetra1.wt + tetra1.ut);
                        tetra1.yt = (tetra1.Vt - tetra1.wt + tetra1.ut) * (tetra1.wt - tetra1.ut + tetra1.Vt);
                        tetra1.Zt = (tetra1.vt - tetra1.Wt + tetra1.ut) * (tetra1.Wt + tetra1.ut + tetra1.vt);
                        tetra1.zt = (tetra1.Wt - tetra1.ut + tetra1.vt) * (tetra1.ut - tetra1.vt + tetra1.Wt);
                        tetra1.at = sqrt(tetra1.xt * tetra1.Yt * tetra1.Zt);
                        tetra1.bt = sqrt(tetra1.yt * tetra1.Zt * tetra1.Xt);
                        tetra1.ct = sqrt(tetra1.zt * tetra1.Xt * tetra1.Yt);
                        tetra1.dt = sqrt(tetra1.xt * tetra1.yt * tetra1.zt);
                        VolPoly += sqrt( (-tetra1.at + tetra1.bt + tetra1.ct + tetra1.dt)*
                                        (tetra1.at - tetra1.bt + tetra1.ct + tetra1.dt)*
                                        (tetra1.at + tetra1.bt - tetra1.ct + tetra1.dt)*
                                        (tetra1.at + tetra1.bt + tetra1.ct - tetra1.dt) )/
                                (192.0 * tetra1.ut * tetra1.vt * tetra1.wt);
                        }
                printf("The total volume of particles is: %.4e m^3\n",VolPar);
                printf("The volume of polyheron is: %.4e m^3\n",VolPoly);
                printf("The packing efficiency is: %.4f\n",VolPar/VolPoly);
                }

        /* for paraview */
        if (bVtk) {
                face_paraview(PARAF,datav,dataf,N_f);
                Sphere_paraview(PARAS,p,N_p,N_p_new);
                if (bVerbose) printf("Open paraview to check the results and have fun!\n");
                }

        free(datav);
        free(dataf);
        free(p);

        return;
}

static void usage(const char *achProgName)
{
        fprintf(stderr,"Usage: %s [-v] [-s] [-r] [-i] [-q] [-a] ss.file polyhedron.file\n"
                "Where: -v = produce .vtk format files (visualization purpose)\n"
                "       -s = the unit in the obj file is m (km by default)\n"
                "       -r = retain the velocity and spin information (set to 0 by default)\n"
                "       -i = retain the aggregate ID information (set to particle ID by default)\n"
                "       -q = activate quiet mode\n"
                "       -a = remove all members if any aggregate constituent crosses a boundary\n",achProgName);
        exit(0);
        }

int
main(int argc,char *argv[])
{
        int i;
        int bVtk = 0;     /* default = don't output .vtk data for visualization in Paraview*/
        int bScale = 1; /* default = km */
        int bRetain = 0;  /* default = don't retain the velocity and spin information (set to 0)*/
        int bAid = 1;     /* default = set to particle ID */
        int bVerbose = 1; /* default = verbose */
	int bAdel = 0;	  /* default = break apart aggs */

        while ( (i = getopt(argc, argv, "vsriqa")) != EOF ) {
                switch (i) {
                case 'v':
                        bVtk = 1;
                        break;
                case 's':
                        bScale = 0;
                        break;
                case 'r':
                        bRetain = 1;
                        break;
                case 'i':
                        bAid = 0;
                        break;
                case 'q':
                        bVerbose = 0;
                        break;
		case 'a':
			bAdel = 1;
			break;
                default:
                        usage(argv[0]);
                        }
                }
        if ( optind >= argc )
                usage(argv[0]);

        /* main process */
        process_poly(argv[optind],argv[optind+1], bVtk, bScale, bRetain, bAid, bVerbose, bAdel);

	return 0;
}
