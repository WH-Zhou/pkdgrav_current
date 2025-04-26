/*
 ** walls2vtp.c -- Yun Zhang 8-11-16
 ** ========
 ** Converts walls data to .vtp format, for visualization in paraview.
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
#include "wallsio.h"

#define VTP_EXT ".vtp"
#define OUTPUT_NOSS "walls.vtp"
#define RESOLUTION 80 /* used for visualizing sphere/disk/cylinder */
#define PLANE_LENGTH 2000 /* the length of infinite plane, unit: m */

typedef struct {
  double dWallTimeOffset;
  double dTime; /* read from input data files */
  int nWalls;
  WALL_DATA *pWalls;
  } WALL_PARAMS;

static void
CoorTran(double dPhi, double dTheta, double scale, double *vIn, double *vMove, double *vOut)  /* coordinate transformation */
{
        double CosPhi,SinPhi,CosTheta,SinTheta;
        CosPhi = cos(dPhi);
        SinPhi = sin(dPhi);
        CosTheta = cos(dTheta);
        SinTheta = sin(dTheta);
        *vOut = *(vMove) + ( SinPhi * (*vIn) + CosPhi * CosTheta * (*(vIn+1)) + CosPhi * SinTheta * (*(vIn+2)) )*scale;
        *(vOut+1) = *(vMove+1) + ( -CosPhi * (*vIn) + SinPhi * CosTheta * (*(vIn+1)) + SinPhi * SinTheta * (*(vIn+2)) )*scale;
        *(vOut+2) = *(vMove+2) + ( -SinTheta * (*(vIn+1)) + CosTheta * (*(vIn+2)) )*scale;
        }

static void
Wall_rectangle_paraview(FILE *fp, const WALL_DATA *w, double *points) /* for WallRectangle type: infinite plane is presented as a rectangle with length of WallRectangle */
{
        fprintf(fp,"    <Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"1\">\n"
                   "      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(fp,"          %.16e %.16e %.16e %.16e %.16e %.16e\n",*points,*(points+1),*(points+2),*(points+3),*(points+4),*(points+5));
        fprintf(fp,"          %.16e %.16e %.16e %.16e %.16e %.16e\n",*(points+6),*(points+7),*(points+8),*(points+9),*(points+10),*(points+11));
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Polys>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
                   "          0 1 3 2\n"
                   "        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
                   "          4\n" /* the number of nodes for each element*/
                   "        </DataArray>\n"
                   "      </Polys>\n"
                   "    </Piece>\n");
        }

static void
Wall_triangle_paraview(FILE *fp, const WALL_DATA *w, double *points) /* for WallTriangle type */
{
        fprintf(fp,"    <Piece NumberOfPoints=\"3\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"1\">\n"
                   "      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(fp,"          %.16e %.16e %.16e %.16e %.16e %.16e\n",*points,*(points+1),*(points+2),*(points+3),*(points+4),*(points+5));
        fprintf(fp,"          %.16e %.16e %.16e\n",*(points+6),*(points+7),*(points+8));
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Polys>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
                   "          0 1 2\n"
                   "        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
                   "          3\n"
                   "        </DataArray>\n"
                   "      </Polys>\n"
                   "    </Piece>\n");
        }

static void
Wall_point_paraview(FILE *fp, const WALL_DATA *w,double *vOrigin)
{
        fprintf(fp,"    <Piece NumberOfPoints=\"1\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
                   "      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(fp,"          %.16e %.16e %.16e\n",*vOrigin,*(vOrigin+1),*(vOrigin+2));
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Verts>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
                   "          0\n"
                   "        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
                   "          1\n"
                   "        </DataArray>\n"
                   "      </Verts>\n"
                   "    </Piece>\n");
        }

static void
Wall_disk_paraview(FILE *fp, const WALL_DATA *w, double *vOrigin, double dPhi, double dTheta) /* for WallDisk type */
{
        VECTOR vTemp,vPoint;
        double dAngle,Arf,dTemp1,dTemp2;
        int i;
        dAngle = 2.0 * M_PI / RESOLUTION;
        vTemp[2] = 0.0;
        dTemp1 = w->dHoleRadius*L_SCALE;
        dTemp2 = w->dRadius*L_SCALE;
        fprintf(fp,"    <Piece NumberOfPoints=\"%i\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%i\">\n",RESOLUTION*2,RESOLUTION);
        fprintf(fp,"      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ ) {
                Arf = dAngle * i;
                vTemp[0] = cos(Arf); vTemp[1] = sin(Arf);
                CoorTran(dPhi, dTheta,dTemp1, &vTemp[0], &vOrigin[0], &vPoint[0]);
                fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                CoorTran(dPhi, dTheta, dTemp2, &vTemp[0], &vOrigin[0], &vPoint[0]);
                fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                }
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Polys>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION-1; i++ )
                   fprintf(fp,"          %i %i %i %i\n",i*2,i*2+1,i*2+3,i*2+2);
        fprintf(fp,"          %i %i 1 0\n",i*2,i*2+1);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for ( i=1; i<=RESOLUTION; i++ )
                   fprintf(fp,"          %i\n",4*i);
        fprintf(fp,"        </DataArray>\n"
                   "      </Polys>\n"
                   "    </Piece>\n");
        }

static void
Wall_line_paraview(FILE *fp, const WALL_DATA *w,double *vOrigin, double dPhi, double dTheta, int bInf) /* for WallCylinder type: radius = 0 case */
{
        VECTOR vTemp,vPoint1,vPoint2;
        double length;
        if (bInf) length = PLANE_LENGTH;
        else length = 0.5 * w->dLength * L_SCALE;
        vTemp[0] = 0.0; vTemp[1] = 0.0; vTemp[2] = 1.0;
        CoorTran(dPhi, dTheta,length, &vTemp[0], &vOrigin[0], &vPoint1[0]);
        vTemp[2] = -1.0;
        CoorTran(dPhi, dTheta,length, &vTemp[0], &vOrigin[0], &vPoint2[0]);
        fprintf(fp,"    <Piece NumberOfPoints=\"2\" NumberOfVerts=\"0\" NumberOfLines=\"1\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
                   "      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(fp,"          %.16e %.16e %.16e\n",vPoint1[0],vPoint1[1],vPoint1[2]);
        fprintf(fp,"          %.16e %.16e %.16e\n",vPoint2[0],vPoint2[1],vPoint2[2]);
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Lines>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
                   "          0 1\n"
                   "        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
                   "          2\n"
                   "        </DataArray>\n"
                   "      </Lines>\n"
                   "    </Piece>\n");
        }

static void
Wall_ring_paraview(FILE *fp, const WALL_DATA *w,double *vOrigin, double dPhi, double dTheta) /* for WallCylinder type: length = 0 case */
{
        VECTOR vTemp,vPoint;
        double dAngle,dTemp,Arf;
        int i;
        dAngle = 2.0 * M_PI / RESOLUTION;
        vTemp[2] = 0.0;
        dTemp = w->dRadius*L_SCALE;
        fprintf(fp,"    <Piece NumberOfPoints=\"%i\" NumberOfVerts=\"0\" NumberOfLines=\"%i\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",RESOLUTION,RESOLUTION);
        fprintf(fp,"      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ ) {
                Arf = dAngle * i;
                vTemp[0] = cos(Arf); vTemp[1] = sin(Arf);
                CoorTran(dPhi, dTheta,dTemp, &vTemp[0], &vOrigin[0], &vPoint[0]);
                fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                }
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ )
                fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ )
                fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ )
                fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Lines>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION-1; i++ )
                   fprintf(fp,"          %i %i\n",i,i+1);
        fprintf(fp,"          %i 0\n",i);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for ( i=1; i<=RESOLUTION; i++ )
                   fprintf(fp,"          %i\n",2*i);
        fprintf(fp,"        </DataArray>\n"
                   "      </Lines>\n"
                   "    </Piece>\n");
        }

static void
Wall_cylinder_paraview(FILE *fp, const WALL_DATA *w, double *vOrigin, double dPhi, double dTheta, int bInf) /* for WallCylinder type (infinite or finite) */
{
        VECTOR vTemp,vOrigin1,vOrigin2,vPoint;
        double dAngle,Arf,dRadius1,dRadius2,length;
        int i;
        dAngle = 2.0 * M_PI / RESOLUTION;
        dRadius1 = w->dRadius * L_SCALE;
        if (bInf) {
                dRadius2 = w->dRadius * L_SCALE;
                length = PLANE_LENGTH;
                }
        else {
                dRadius2 = w->dRadius * (1.0 - w->dTaper) * L_SCALE;
                length = 0.5 * w->dLength * L_SCALE;
                }
        vOrigin1[0] = vOrigin[0] - w->vOrient[0] * length;
        vOrigin1[1] = vOrigin[1] - w->vOrient[1] * length;
        vOrigin1[2] = vOrigin[2] - w->vOrient[2] * length;
        vOrigin2[0] = vOrigin[0] + w->vOrient[0] * length;
        vOrigin2[1] = vOrigin[1] + w->vOrient[1] * length;
        vOrigin2[2] = vOrigin[2] + w->vOrient[2] * length;
        vTemp[2] = 0.0;
        fprintf(fp,"    <Piece NumberOfPoints=\"%i\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%i\">\n",RESOLUTION*2,RESOLUTION);
        fprintf(fp,"      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ ) {
                Arf = dAngle * i;
                vTemp[0] = cos(Arf); vTemp[1] = sin(Arf);
                CoorTran(dPhi, dTheta, dRadius1, &vTemp[0], &vOrigin1[0], &vPoint[0]);
                fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                CoorTran(dPhi, dTheta, dRadius2, &vTemp[0], &vOrigin2[0], &vPoint[0]);
                fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                }
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*2; i++ )
                fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Polys>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION-1; i++ )
                   fprintf(fp,"          %i %i %i %i\n",i*2,i*2+1,i*2+3,i*2+2);
        fprintf(fp,"          %i %i 1 0\n",i*2,i*2+1);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for ( i=1; i<=RESOLUTION; i++ )
                   fprintf(fp,"          %i\n",4*i);
        fprintf(fp,"        </DataArray>\n"
                   "      </Polys>\n"
                   "    </Piece>\n");
        }

static void
Wall_shell_paraview(FILE *fp, const WALL_DATA *w, double *vOrigin, double dPhi, double dTheta) /* for WallShell type */
{
        VECTOR vTemp,vOriginZ,vPoint;
        double dAngle,dAngleZ,dOpenAngle,Arf,ArfZ,dRadius,dRadiusRing,dTemp;
        int i,j;
        dOpenAngle = w->dOpenAngle / 180.0 * M_PI;
        dAngle = 2.0 * M_PI / RESOLUTION;
        dAngleZ = (M_PI - dOpenAngle) / (RESOLUTION-1);
        dRadius = w->dRadius * L_SCALE;
        vTemp[2] = 0.0;
        fprintf(fp,"    <Piece NumberOfPoints=\"%i\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%i\">\n",RESOLUTION*RESOLUTION,RESOLUTION*RESOLUTION-RESOLUTION);
        fprintf(fp,"      <Points>\n"
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION; i++ ) {
                ArfZ = dOpenAngle + dAngleZ *i;
                dRadiusRing = dRadius * sin(ArfZ);
                dTemp = dRadius * cos(ArfZ);
                vOriginZ[0] = vOrigin[0] + w->vOrient[0] * dTemp;
                vOriginZ[1] = vOrigin[1] + w->vOrient[1] * dTemp;
                vOriginZ[2] = vOrigin[2] + w->vOrient[2] * dTemp;
                for ( j=0; j<RESOLUTION; j++ ) {
                        Arf = dAngle * j;
                        vTemp[0] = cos(Arf); vTemp[1] = sin(Arf);
                        CoorTran(dPhi, dTheta, dRadiusRing, &vTemp[0], &vOriginZ[0], &vPoint[0]);
                        fprintf(fp,"          %.16e %.16e %.16e\n",vPoint[0],vPoint[1],vPoint[2]);
                        }
                }
        fprintf(fp,"        </DataArray>\n"
                   "      </Points>\n"
                   "      <CellData>\n"
                   "        <DataArray type=\"Int16\" Name=\"Color\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*RESOLUTION; i++ )
                fprintf(fp,"          %i\n",w->iColor);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Float32\" Name=\"Opacity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*RESOLUTION; i++ )
                fprintf(fp,"          %.8e\n",w->dTrans);
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int16\" Name=\"Assembly\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION*RESOLUTION; i++ )
                fprintf(fp,"          %i\n",w->iAssembly);
        fprintf(fp,"        </DataArray>\n"
                   "      </CellData>\n"
                   "      <Polys>\n"
                   "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        for ( i=0; i<RESOLUTION-1; i++ ) {
                for ( j=0; j<RESOLUTION-1; j++ )
                       fprintf(fp,"          %i %i %i %i\n",i*RESOLUTION+j,i*RESOLUTION+j+1,(i+1)*RESOLUTION+j+1,(i+1)*RESOLUTION+j);
                fprintf(fp,"          %i %i %i %i\n",i*RESOLUTION+j,i*RESOLUTION,(i+1)*RESOLUTION,(i+1)*RESOLUTION+j);
                }
        fprintf(fp,"        </DataArray>\n"
                   "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for ( i=1; i<=RESOLUTION*RESOLUTION-RESOLUTION; i++ )
                   fprintf(fp,"          %i\n",4*i);
        fprintf(fp,"        </DataArray>\n"
                   "      </Polys>\n"
                   "    </Piece>\n");
        }

static void draw_wall(FILE *fp,const WALL_PARAMS *p,WALL_DATA *w,int bReact,char *WALLREACT,int nStep,int WALL_ID)
{
        VECTOR vOrigin,vTravel,vTemp;
        FILE *fpW = NULL;
        double dPhi,dTheta,dTime,points[4][3],dTemp;
        int i,iTemp;
        int c;

        if (bReact && w->dMass != 0) {
                if ((fpW = fopen(WALLREACT,"r")) == NULL) {
                        fprintf(stderr,"Unable to open \"%s\".\n",WALLREACT);
                        if (fpW != NULL)
                                fclose(fpW);
                        exit(1);
                        }
                fscanf(fpW,"%i",&iTemp);
                if (iTemp != nStep)
                        while ( (c=fgetc(fpW)) != EOF )
                                if (c == '\n'){
                                        fscanf(fpW,"%i",&iTemp);
                                        if (iTemp == nStep) break;
                                        }
                for (i=0;i<WALL_ID;i++) {
                        if (p->pWalls[i].iType == WallTriangle || p->pWalls[i].iType == WallRectangle )
                                fscanf(fpW,"%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",&iTemp,&iTemp,
                                       &vTemp[0],&vTemp[1],&vTemp[2],&vTemp[0],&vTemp[1],&vTemp[2],&vTemp[0],&vTemp[1],&vTemp[2]);
                        else
                                fscanf(fpW,"%i%i%lf%lf%lf%lf%lf%lf",&iTemp,&iTemp,
                                       &vTemp[0],&vTemp[1],&vTemp[2],&vTemp[0],&vTemp[1],&vTemp[2]);
                        }
                if (w->iType == WallTriangle || w->iType == WallRectangle )
                        fscanf(fpW,"%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",&iTemp,&iTemp,&vOrigin[0],&vOrigin[1],&vOrigin[2],
                                &w->vVertex1[0],&w->vVertex1[1],&w->vVertex1[2],&w->vVertex2[0],&w->vVertex2[1],&w->vVertex2[2]);
                else
                        fscanf(fpW,"%i%i%lf%lf%lf%lf%lf%lf",&iTemp,&iTemp,&vOrigin[0],&vOrigin[1],&vOrigin[2],
                                &w->vOrient[0],&w->vOrient[1],&w->vOrient[2]);
                /*printf("1. WALL_ID: %i, iTemp: %i, vOrigin:%.8e %.8e %.8e\n",WALL_ID,iTemp,vOrigin[0],vOrigin[1],vOrigin[2]);*/ /* DEBUG!!! Yun */
                fclose(fpW);
                }
        else {
                COPY_VEC(w->vOrigin,vOrigin); /* so we can modify it in the case of moving walls */
                /* account for any wall motion */
                dTime = (w->dStep == 0.0 ? p->dTime : ((int) (p->dTime/w->dStep))*w->dStep);
                COPY_VEC(w->vVel,vTravel);
                SCALE_VEC(vTravel,dTime);
                ADD_VEC(vOrigin,vTravel,vOrigin);
                /* and oscillation */
                COPY_VEC(w->vOscVec,vTravel);
                SCALE_VEC(vTravel,w->dOscAmp*sin(w->dOscFreq*dTime + w->dOscPhase));
                ADD_VEC(vOrigin,vTravel,vOrigin);
                }
        SCALE_VEC(vOrigin,L_SCALE); /* unit: m */

        /* rotation angles (not needed in all cases) */
        dPhi = atan2(w->vOrient[Y],w->vOrient[X]);
        dTheta = atan2(sqrt(w->vOrient[X]*w->vOrient[X] + w->vOrient[Y]*w->vOrient[Y]),w->vOrient[Z]);

        /* convert to .vtk format based on the wall's type */
        switch (w->iType) {
        case WallPlane: {
                dTemp = PLANE_LENGTH;
                vTemp[0] = -1; vTemp[1] = -1; vTemp[2] = 0;
                CoorTran(dPhi, dTheta, dTemp, &vTemp[0], &vOrigin[0], &points[0][0]);
                vTemp[0] = 1;
                CoorTran(dPhi, dTheta, dTemp, &vTemp[0], &vOrigin[0], &points[1][0]);
                vTemp[0] = -1; vTemp[1] = 1;
                CoorTran(dPhi, dTheta, dTemp, &vTemp[0], &vOrigin[0], &points[2][0]);
                vTemp[0] = 1;
                CoorTran(dPhi, dTheta, dTemp, &vTemp[0], &vOrigin[0], &points[3][0]);
                Wall_rectangle_paraview(fp, w, &points[0][0]);
                break;
                }
        case WallTriangle: {
                COPY_VEC( vOrigin, &points[0][0] );
                COPY_VEC(w->vVertex1,vTemp);
                SCALE_VEC(vTemp,L_SCALE);
                ADD_VEC( vOrigin, vTemp, &points[1][0] );
                COPY_VEC(w->vVertex2,vTemp);
                SCALE_VEC(vTemp,L_SCALE);
                ADD_VEC( vOrigin, vTemp, &points[2][0] );
                Wall_triangle_paraview(fp, w, &points[0][0]);
                break;
                }
        case WallRectangle: {
                COPY_VEC( vOrigin, &points[0][0] );
                COPY_VEC(w->vVertex1,vTemp);
                SCALE_VEC(vTemp,L_SCALE);
                ADD_VEC( vOrigin, vTemp, &points[1][0] );
                COPY_VEC(w->vVertex2,vTemp);
                SCALE_VEC(vTemp,L_SCALE);
                ADD_VEC( vOrigin, vTemp, &points[2][0] );
                ADD_VEC( &points[1][0], vTemp, &points[3][0] );
                Wall_rectangle_paraview(fp, w, &points[0][0]);
                break;
                }
        case WallDisk: {
                if (w->dRadius == 0.0)  /* degenerate case: point */
                        Wall_point_paraview(fp,w,&vOrigin[0]);
                else
                        Wall_disk_paraview(fp,w,&vOrigin[0],dPhi,dTheta);
                break;
                }
        case WallCylinderInfinite: {
                if (w->dRadius == 0.0)   /* degenerate case: line */
                        Wall_line_paraview(fp,w,&vOrigin[0],dPhi,dTheta,1);
                else
                        Wall_cylinder_paraview(fp,w,&vOrigin[0],dPhi,dTheta,1);
                break;
                }
        case WallCylinderFinite: {
                if (w->dRadius == 0.0 && w->dLength == 0.0)   /* degenerate case: line */
                        Wall_point_paraview(fp,w,&vOrigin[0]);
                else if (w->dRadius == 0.0)
                        Wall_line_paraview(fp,w,&vOrigin[0],dPhi,dTheta,0);
                else if (w->dLength == 0.0)
                        Wall_ring_paraview(fp,w,&vOrigin[0],dPhi,dTheta);
                else
                        Wall_cylinder_paraview(fp,w,&vOrigin[0],dPhi,dTheta,0);
                break;
                }
        case WallShell: {
                if (w->dRadius == 0.0 || w->dOpenAngle == 180.0)   /* degenerate case: point */
                        Wall_point_paraview(fp,w,&vOrigin[0]);
                else
                        Wall_shell_paraview(fp,w,&vOrigin[0],dPhi,dTheta);
                break;
                }
        default:
                assert(0);
                }
        }

static int
convert(int argc, char *infile[], int bReact, int bNeedss, int bVerbose)
{
        SSIO ssio;
        SSHEAD h;
        FILE *fp = NULL;
        double dTime;
        char outfile[MAXPATHLEN],CalStep[MAXPATHLEN];
        int i,j,k,nStep = 0,ibegin = 1;
        WALL_PARAMS *pW,WallP;

        /* read wall file */
        pW = &WallP;
        assert(pW != NULL);
        pW->nWalls = 0;
		/* note optind is declared globally... */
        if (bVerbose) printf("Reading wall data from \"%s\"...\n",infile[optind]);
        if ((fp = fopen(infile[optind],"r")) == NULL) {
                fprintf(stderr,"Unable to open \"%s\".\n",infile[optind]);
                if (fp != NULL)
                        fclose(fp);
                exit(1);
                }
        if (wallsParseWallsFile(fp,&pW->nWalls,&pW->pWalls,&dTime,bVerbose) != 0) {
                fprintf(stderr,"Error occurred while parsing wall data.\n");
                if (fp != NULL)
                        fclose(fp);
                exit(1);
                }
        fclose(fp);
        if (bVerbose) printf("Number of walls read = %i.\n",pW->nWalls);
        if (pW->dWallTimeOffset == 0.0)
                pW->dWallTimeOffset = dTime;
        /* main convert process */
        if (bNeedss) {
                if (bReact) ibegin = 2;
                for ( i=ibegin+optind;i<argc;++i ) {
                        /* Get the time */
                        if (ssioOpen(infile[i],&ssio,SSIO_READ)) {
                                (void) fprintf(stderr,"Unable to open \"%s\"\n",infile[i]);
                                return 1;
                                }
                        if (ssioHead(&ssio,&h) || h.n_data < 0) {
                                (void) fprintf(stderr,"Corrupt header\n");
                                (void) ssioClose(&ssio);
                                return 1;
                                }
                        if (bVerbose) printf("%s: time %g\n",infile[i],h.time);
                        pW->dTime = h.time + pW->dWallTimeOffset;
                        if ( ssioNewExt(infile[i],SS_EXT,outfile,VTP_EXT) ) {
                                fprintf(stderr,"Unable to generate output .vtp file filename.\n");
                                return 1;
                                }
                        if (bReact) {
                                k = 0;
                                for ( j=0;j<(int) strlen(infile[i]);j++) {
                                        if (infile[i][j]>='0' && infile[i][j]<='9' && k > 0) {
                                                CalStep[k] = infile[i][j];
                                                k++;
                                                }
                                        if (infile[i][j]>='1' && infile[i][j]<='9' && k == 0) {
                                                CalStep[k] = infile[i][j];
                                                k++;
                                                }
                                        }
                                nStep = atoi(CalStep);
                                }
                        if ((fp = fopen(outfile,"w")) == NULL) {
                                fprintf(stderr,"Unable to open \"%s\".\n",outfile);
                                if (fp != NULL)
                                        fclose(fp);
                                exit(1);
                                }
                        fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
                        fprintf(fp,"  <PolyData>\n");
                        for ( j=0;j<pW->nWalls;j++ )
                                draw_wall(fp,pW,&pW->pWalls[j],bReact,infile[1+optind],nStep,j);
                        fprintf(fp,"  </PolyData>\n");
                        fprintf(fp,"</VTKFile>\n");
                        fclose(fp);
                        }
                }
        else {
                pW->dTime = pW->dWallTimeOffset;
                nStep = 0;
                if ((fp = fopen(OUTPUT_NOSS,"w")) == NULL) {
                        fprintf(stderr,"Unable to open \"%s\".\n",OUTPUT_NOSS);
                        if (fp != NULL)
                                fclose(fp);
                        exit(1);
                        }
                fprintf(fp,"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
                fprintf(fp,"  <PolyData>\n");
                for ( j=0;j<pW->nWalls;j++ )
                        draw_wall(fp,pW,&pW->pWalls[j],bReact,infile[1+optind],nStep,j);
                fprintf(fp,"  </PolyData>\n");
                fprintf(fp,"</VTKFile>\n");
                fclose(fp);
                }

        return 0;
        }

static void usage(const char *achProgName)
{
        fprintf(stderr,"Usage: %s [ -r ] [ -n ] [ -q ] wall.dat [ wallreact file ] [ss.file] \n"
                "Note: the order of input file should be the same as described in the usage\n"
                "Where: -r = use walls react (must specify wallreact file)\n"
                "       -n = no ss.file is needed (time is set to 0)\n"
                "       -q = activate quiet mode\n"
                "       output wall.vtk file corresponding to the given ss.file\n",achProgName);
        exit(0);
        }

int
main(int argc,char *argv[])
{
        int i;
        int bReact = 0;    /* default = don't use walls react */
        int bNeedss = 1;   /* default = need input ss.file */
        int bVerbose = 1;  /* default = verbose */
        int flag = 0;  /* default = get correct output */

        while ( (i = getopt(argc, argv, "r::n::q::")) != EOF ) {
                switch (i) {
                case 'r':
                        bReact = 1;
                        break;
                case 'n':
                        bNeedss = 0;
                        break;
                case 'q':
                        bVerbose = 0;
                        break;
                default:
                        usage(argv[0]);
                        }
                }
        if ( optind >= argc )
                usage(argv[0]);
        flag = convert( argc, argv, bReact, bNeedss, bVerbose );
        if ( bVerbose && !flag) printf("Open paraview to check the results and have fun!\n");
        return 0;
	}

/* walls2vtp.c */
