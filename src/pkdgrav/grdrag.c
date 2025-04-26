/* suppress pedantic "empty translation unit" warning if !GR_DRAG */
typedef int make_iso_compilers_happy;

#ifdef GR_DRAG

/* Program to compute the new position and velocity for a 
test particle around a massive object, integrated correctly
according to dynamics in the Schwarzschild spacetime.  This
returns the new positions and velocities (in units of AU and
2pi AU/yr) after some time interval measured in units of yr/2pi.
If the test particle merges with the massive object, the final
argument ("merge") equals 1.  Otherwise it is 0. 

In this version I attempt to include gravitational radiation
energy and angular momentum fluxes to lowest order, from 
Kidder 1995, PRD, 52, 821, equations (3.25a) and (3.28a),
respectively. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PI 3.14159265359
#define AU 1.49598e13 /* 1 AU in cm */
#define C 2.99792458e10 /* C in cm/s */
#define YR 31558207.07 /* 1 YR in s */ /* different from google */
#define MSUN 1.98892e33 /* Solar mass in g */
#define G 6.6726e-8 /* G in CGS */

void grdrag(double px0, double py0, double pz0, double vx0, double
   vy0, double vz0, double m, double M, double interval, double *pxfin,
   double *pyfin, double *pzfin, double *vxfin, double *vyfin, double *vzfin,
	    int *merge, double *apocenter, double *pericenter, double subcode_r2, double *t_final, int returnr_p);
void rk4(double *y, double *dydt, int n, double h, double *yout);
void rkqc(double *y,double *dydt,int n,double t,double htry,double eps,
          double *yscal,double *hdid,double *hnext);
void derivs(double *y, double *dydt);
double pickpericenter(double r1, double r2, double r3);


void grdrag(double px0, double py0, double pz0, double vx0, double
   vy0, double vz0, double m, double M, double interval, double *pxfin,
   double *pyfin, double *pzfin, double *vxfin, double *vyfin, double *vzfin,
	    int *merge, double *apocenter, double *pericenter, double subcode_r2, double *t_final, int returnr_p)
{
  int i,j;
  double r0,r,phi,lx,ly,lz,u_phi,uphi,u_t,ut,ur,vr,vphi;
  double t,dt,mgrav,T;

  double p0mag,px0unit,py0unit,pz0unit,lxunit,lyunit,lzunit; /* This next group is used for coordinate transformations */
  double oxunit,oyunit,ozunit, phixunit,phiyunit,phizunit;

  double eps,hdid,hnext,*y,*dydt,*yscal;
  double *yold, *dydtold, *yscalold; /* These are to make the particle exit at *exactly* the right radius */
  double urold,ur0;
  double Enewt,f,v2,v02;
  double cubicQ, cubicR, cubicEnergy, theta, cubicA, cubicB, cubicC, root1, root2, root3; /* These are all for the pericenter calculation */
  
   y=calloc(7,sizeof(double));
   yold=calloc(7,sizeof(double));
   dydt=calloc(7,sizeof(double));
   dydtold=calloc(7,sizeof(double));
   yscal=calloc(7,sizeof(double));
   yscalold=calloc(7,sizeof(double));
   eps=1.0e-6;

/* First step: convert from pkdgrav units (AU, 2pi AU/yr, masses
in solar masses, times in yr/2pi) to geometrized units in which
distances and times are in units of M=GM/c^2 and speeds are in
units of c. */

   mgrav=M*MSUN*G/(C*C); /* Mass in cm */
   T=interval*YR/(2*PI)/(mgrav/C); /* Time in M */
   px0*=AU;
   py0*=AU;
   pz0*=AU;   /* Now distances are in cm. */
   vx0*=2.0*PI*AU/YR;
   vy0*=2.0*PI*AU/YR;
   vz0*=2.0*PI*AU/YR; /* Now speeds are in cm/s. */

   px0/=mgrav;
   py0/=mgrav;
   pz0/=mgrav;  /* Now distances are in units of M. */
   vx0/=C;
   vy0/=C;
   vz0/=C;  /* Speeds are now in units of c. */

   r0=sqrt(px0*px0+py0*py0+pz0*pz0); /* Distance from massive object in M*/
   Enewt=1.0+0.5*(vx0*vx0+vy0*vy0+vz0*vz0)-1.0/(r0-3);
   lx=py0*vz0-pz0*vy0;
   ly=pz0*vx0-px0*vz0;
   lz=px0*vy0-py0*vx0;
   u_phi=sqrt(lx*lx+ly*ly+lz*lz);   /* Specific angular momentum, units of M. */
/* Now construct two unit vectors, one along r and the other perpendicular
to r but in the orbital plane, for later use in translating the
angle traveled to Cartesian coordinates. */
   p0mag=sqrt(px0*px0+py0*py0+pz0*pz0);
   px0unit=px0/p0mag;
   py0unit=py0/p0mag;
   pz0unit=pz0/p0mag;
   lxunit=lx/u_phi;
   lyunit=ly/u_phi;
   lzunit=lz/u_phi;
   oxunit=pz0unit*lyunit-py0unit*lzunit;
   oyunit=px0unit*lzunit-pz0unit*lxunit;
   ozunit=py0unit*lxunit-px0unit*lyunit;
   
   ur=(px0*vx0+py0*vy0+pz0*vz0)/r0; /* Radial speed, units of c. */
   f=Enewt*Enewt/(1.0-2.0/r0)-1.0;
   f/=(u_phi*u_phi/(r0*r0)+ur*ur/(1.0-2.0/r0));
   f=sqrt(f);
   u_phi*=f;
   ur*=f;
   u_t=-Enewt;  /* Specific energy relative to mc^2 */
   y[0]=r0;
   y[1]=0.0;
   y[2]=ur;
   y[3]=u_phi;
   y[4]=u_t;
   y[5]=m*M/(m+M);
   y[6]=m+M;
   derivs(y, dydt);
   yscal[0]=y[0];
   yscal[1]=1.0;
   yscal[2]=1.0e-4; 
   yscal[3]=1.0;
   yscal[4]=1.0;
   yscal[5]=1.0;
   yscal[6]=1.0;
   dt=0.0001*T;
   t=0.0;
   *merge=0;
   *apocenter=0;
   *t_final=0.0;
   urold=ur;
   r=y[0];

   /* If we get the signal from the GRICP function that this is the first entry, we calculate and return the pericenter */
   if (returnr_p) {
   cubicEnergy = (u_t*u_t - 1.0)/2.0; /* Script E in Hartle's GR book */
   cubicA = 1.0/cubicEnergy; /* Values of the constants in the cubic equation */
   cubicB = -0.5*u_phi*u_phi/cubicEnergy; 
   cubicC = u_phi*u_phi/cubicEnergy;
   cubicQ = (cubicA*cubicA - 3*cubicB)/9.0;
   cubicR = (2.0*cubicA*cubicA*cubicA - 9.0*cubicA*cubicB + 27.0*cubicC)/54.0;
   /* We check to see if we have 3 real roots */
   theta = acos(cubicR/sqrt(cubicQ*cubicQ*cubicQ));
   root1 = -2.0*sqrt(cubicQ)*cos(theta/3.0) - cubicA/3.0;
   root2 = -2.0*sqrt(cubicQ)*cos((theta+2.0*PI)/3.0) - cubicA/3.0;
   root3 = -2.0*sqrt(cubicQ)*cos((theta-2.0*PI)/3.0) - cubicA/3.0;
   /* Now we choose the pericenter (the smallest positive root) */
   assert((root1 > 0) || (root2 > 0) || (root3 > 0));
   *pericenter = pickpericenter(root1, root2, root3)*mgrav/AU;
   }

   for (; t<T; )
   {
     for (i = 0 ; i < 7 ; ++i){
       yold[i] = y[i];
       dydtold[i] = dydt[i];
       yscalold[i] = yscal[i];
       }
      rkqc(y,dydt,7,t,dt,eps,yscal,&hdid,&hnext);
      r=y[0];
      ur0=y[2];
      ut=-y[4]/(1-2.0/y[0]);
      ur=(1.0-2.0/r)*(-1.0-ut*y[4]-u_phi*u_phi/(r*r));
      if (ur<0.0) 
      {
        ur=1.0e-9;
        if (r>r0) ur=-1.0e-9;
      }
      else
      {
        ur=sqrt(ur);
        if (ur0<0.0) ur*=-1.0;
      }
      y[2]=ur; 
      ur=y[2];


      if (r*mgrav/AU > sqrt(subcode_r2)) /* If we've gone back outside our entry radius we need to find out the time at which we get to exactly the entry radius */
	{
	  for (j = 0 ; j < 4 ; ++j) {
	    dt *= (sqrt(subcode_r2) - yold[0]*mgrav/AU)/(y[0]*mgrav/AU-yold[0]*mgrav/AU);
	    for (i = 0 ; i < 7 ; ++i){
	      y[i] = yold[i];
	      dydt[i] = dydtold[i];
	      yscal[i] = yscalold[i];
	    }
	    rkqc(y, dydt, 7, t, dt, eps, yscal, &hdid, &hnext);
	  }
	  y[0] = sqrt(subcode_r2)*AU/mgrav;
	  t+=dt*ut;
	  *t_final = t*mgrav*2.0*PI/(YR*C);
	  break;
	}



      /* Checks to see if we're getting apocenter inside the subcode radius - this shouldn't even be possible without the energy loss terms on,
         and it should be extremely unlikely even at that.  */
      
      if ((urold > 0) && (ur < 0))
	{
	  *apocenter=y[0];
	  *apocenter*=mgrav/AU;  /* Converts it back to AU */
 	}
	
      urold=ur;
      if ((r<3.0) || (*apocenter > 0.0))
      {
        *merge=1;
        *pxfin=0.0;
        *pyfin=0.0;
        *pzfin=0.0;
        *vxfin=0.0;
        *vyfin=0.0;
        *vzfin=0.0;
        break;
      }

      t+=dt*ut;
      dt=hnext;
      if (t+dt>T) 
	dt=T-t;
      phi=y[1];


   }
   u_phi=y[3];
   u_t=y[4];

   r=y[0];
   phi=y[1];
   ur=y[2];
   ut=-y[4]/sqrt(1-2.0/r);
   uphi=y[3]/(r*r);
   vr=ur;
   vphi=uphi;

   /* Now we convert from theta-phi space back to x-y-z coordinates */

   *pxfin=(px0unit*cos(phi)+oxunit*sin(phi));
   *pyfin=(py0unit*cos(phi)+oyunit*sin(phi));
   *pzfin=(pz0unit*cos(phi)+ozunit*sin(phi));
   phixunit=*pzfin*lyunit-*pyfin*lzunit;
   phiyunit=*pxfin*lzunit-*pzfin*lxunit;
   phizunit=*pyfin*lxunit-*pxfin*lyunit;
  
   *vxfin=vr*(*pxfin)+r*vphi*phixunit;
   *vyfin=vr*(*pyfin)+r*vphi*phiyunit;
   *vzfin=vr*(*pzfin)+r*vphi*phizunit;
   v2=2.0*(-u_t-1.0+1.0/(r-3));
   v02=*vxfin*(*vxfin)+*vyfin*(*vyfin)+*vzfin*(*vzfin);
   *vxfin*=sqrt(v2/v02);
   *vyfin*=sqrt(v2/v02);
   *vzfin*=sqrt(v2/v02);

   /* Now we convert back to code units */
   *vxfin=*vxfin*C*YR/(2*PI*AU);
   *vyfin=*vyfin*C*YR/(2*PI*AU);
   *vzfin=*vzfin*C*YR/(2*PI*AU);
   *pxfin*=r*mgrav/AU;
   *pyfin*=r*mgrav/AU;
   *pzfin*=r*mgrav/AU; 

   free(y);
   free(dydt);
   free(yscal);
   free(dydtold);
   free(yold);
}

void rk4(double *y, double *dydt, int n, double h, double *yout)
{
/* Fourth-order Runge-Kutta integrator.  *y is an array with the n variables;
*dydt is the array of their derivatives; *yout is the array of values at
t+dt. */

   int i;
   double hh,h6,*dym,*dyt,*yt;
  
   dym=malloc(n*sizeof(double));
   dyt=malloc(n*sizeof(double));
   yt=malloc(n*sizeof(double));
   
   hh=0.5*h;
   h6=h/6.0;
   
   for (i=0; i<n; i++)
      yt[i]=y[i]+hh*dydt[i];
   derivs(yt,dyt);
   for (i=0; i<n; i++)
      yt[i]=y[i]+hh*dyt[i];
   derivs(yt,dym);
   for (i=0; i<n; i++)
   {
      yt[i]=y[i]+h*dym[i];
      dym[i]+=dyt[i];
   }
   derivs(yt,dyt);
   for (i=0; i<n; i++)
      yout[i]=y[i]+h6*(dydt[i]+dyt[i]+2.0*dym[i]);
   free(dym);
   free(dyt);
   free(yt);
}

void rkqc(double *y,double *dydt,int n,double t,double htry,double eps,
          double *yscal,double *hdid,double *hnext)
{
/* Quality-controlled fifth-order Runge-Kutta. */
   
   int i;
   double tsav,hh,h,temp,errmax;
   double *dysav,*ysav,*ytemp;
   double PGROW,PSHRNK,FCOR,SAFETY,ERRCON;

   PGROW=-0.20;
   PSHRNK=-0.25;
   FCOR=1.0/15.0;
   SAFETY=0.9;
   ERRCON=6.0e-4;

   dysav=malloc(n*sizeof(double));
   ysav=malloc(n*sizeof(double));
   ytemp=malloc(n*sizeof(double));
   tsav=t;
   for (i=0; i<n; i++)
   {
      ysav[i]=y[i];
      dysav[i]=dydt[i];
   }
   h=htry;
   for (;;)
   {
      hh=0.5*h;
      rk4(ysav,dysav,n,hh,ytemp);
      t=tsav+hh;
      derivs(ytemp,dydt);
      rk4(ytemp,dydt,n,hh,y);
      t=tsav+h;
      if (t==tsav) printf("Step size too small in RKQC");
      rk4(ysav,dysav,n,h,ytemp);
      errmax=0.0;
      for (i=0; i<n; i++)
      {
         ytemp[i]=y[i]-ytemp[i];
	 temp=fabs(ytemp[i]/yscal[i]);
	 if (errmax<temp) errmax=temp;
      }
      errmax/=eps;
      if (errmax<1.0)
      {
        *hdid=h;
	*hnext=(errmax>ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
	break;
      }
      h=SAFETY*h*exp(PSHRNK*log(errmax));
   }
   for (i=0; i<n; i++) 
      y[i]+=ytemp[i]*FCOR;
      free(ytemp);
      free(dysav);
      free(ysav);
}

void derivs(double *y, double *dydt)
{
    double v;

    dydt[0]=y[2];
    dydt[1]=y[3]/(y[0]*y[0]);
    dydt[2]=-1.0/(y[0]*y[0])+(1.0-3.0/y[0])*y[3]*y[3]/(y[0]*y[0]*y[0]);
    v=sqrt(y[2]*y[2]+y[3]*y[3]/(y[0]*y[0]));
    dydt[3]=-1.6*(y[5]/y[6])*y[3]*(2.0*v*v-3.0*y[2]*y[2]+2.0/y[0]);
    dydt[3]/=(y[0]*y[0]*y[0]);
    dydt[4]=(8.0/15.0)*(y[5]/y[6])/(y[0]*y[0]*y[0]*y[0]);
    dydt[4]*=(12.0*v*v-11.0*y[2]*y[2]);
    /*    dydt[3] = 0.0;
	  dydt[4] = 0.0;  Set these to 0 to turn off the energy loss terms */

/* Note: lack of negative sign in dydt[4] is because u_t=-e. */
    dydt[5]=0.0;
    dydt[6]=0.0;
}

double pickpericenter(double r1, double r2, double r3)
{
  /* This function takes in the 3 roots from the cubic equation and figures
     out which one of them corresponds to pericenter.  We assume that 
     one of the roots is the pericenter, another the apocenter...*/

  double tempr;
  
  if (r1 < 3.0)
    tempr = r2 < r3 ? r2 : r3;
  else if (r2 < 3.0)
    tempr = r1 < r3 ? r1 : r3;
  else if (r3 < 3.0)
    tempr = r1 < r2 ? r1 : r2;

  /* If none of these are true, all three roots are > 3.0 and we return the smallest */
  else {
    if (r1 < r2)
      tempr = r1 < r3 ? r1 : r3;
    else
      tempr = r2 < r3 ? r2 : r3;
  }
   
  return tempr;
}

#endif /* GR_DRAG */
