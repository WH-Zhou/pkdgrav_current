/* suppress pedantic "empty translation unit" warning if !COLLMOD_ZML */
typedef int make_iso_compilers_happy;

#ifdef COLLMOD_ZML

/* ZML 03.20.12 - resolves collisions in a planetesimal disk using collision model from L&S ApJ 2012 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /*for getpid()*/
#include <math.h>
#include <assert.h>
#include <time.h>
#include "collmod.h"
#include "collision.h"
#include "random.h" /*DEBUG only used for some generators below---reconcile*/

#ifndef PI
#define PI M_PI
#endif

#ifndef TWO_PI
#define TWO_PI (2*PI)
#endif

#define TINY 1.0e-14 /* it's crashed with 1.0e-15 on Linux-alpha -- DCR */

#define COLLMOD_NOFRIC 0 /* No dynamical friction or gas friction */
#define COLLMOD_DYNAMICAL 1 /* Just dynamical friction */

// fixes suggested by zoe XXX
#define COLLMOD_GAS 3
#define COLLMOD_DYNANDGAS 4

#define COLLMOD_FRICTION COLLMOD_NOFRIC	
#define EXPANSION_FACTOR 6.0 /* BEWARE - if you have expanded the particles - this should not be 1 - currently only used in the drag calcs! */

#define MAX_PART 1000
#define AU_CGS 14959887070000
#define YEAR 31557600.0
#define LEN_CGS 1.5e13
#define MASS_CGS 1.99e33
#define TIME_CGS (3.14e7/TWO_PI)
#define VEL_CGS (LEN_CGS/TIME_CGS)

/* pkd routines */

void
collLocateBin(const PARTICLE *p,const DUST_BINS_PARAMS *DB,int *piBin,double *r)
{
	// Determines which dustBin a particle is in.
	
	double r_xy;

	if (DB->nDustBins == 0) { /* all dust goes to trash if no bins */
		*piBin = -1;
		return ;
		}

	r_xy = sqrt(SQ(p->r[0]) + SQ(p->r[1]));
	
	// Now calculate which DustBin particle is in
	*piBin = (r_xy - DB->dDustBinsInner)/DB->dDustBinsWidth;

	// If particle location is equal to the boundary of
	// The outer edge of the DustBin domain then set
	// DustBin value as the penultimate one
	
	if (r_xy == DB->dDustBinsOuter)
		*piBin = DB->nDustBins - 1; /* very unlikely! */

	if (r != NULL)
		*r = r_xy;
}

double pkdDustBinsInclAvg(PKD pkd)
{

    // Determines the average inclination of the DustBins
    // Based on particle eccentricities

	PARTICLE *p;
	double InclTot = 0.0,hx,hy,hz,h,incl;
	int i,nLocal = pkdLocal(pkd);
	
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		
		// Cross-Product to work out the value of
		// Specific Relative Angular Momentum
		hx = p->r[1]*p->v[2] - p->v[1]*p->r[2];
		hy = p->v[0]*p->r[2] - p->r[0]*p->v[2];
		hz = p->r[0]*p->v[1] - p->v[0]*p->r[1];
		h = sqrt(hx*hx + hy*hy + hz*hz);
		
		// Find inclination of particle by doing
		// Arc Cos of z-SPAM over SPAM total
		incl = acos(hz/h);
		InclTot += incl;
	}
	return InclTot;
}

double pkdDustBinsInclMax(PKD pkd)
{

    // Determines the maximum inclination of the DustBins
    // Based on particle eccentricities
    
	PARTICLE *p;
	double incl = 0.0, InclMax = 0.0,sin_incl2,h,hx,hy,hz;
	int i,nLocal = pkdLocal(pkd);
	
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		
		// Cross-Product to work out the value of
		// Specific Relative Angular Momentum

		hx = p->r[1]*p->v[2] - p->v[1]*p->r[2];
		hy = p->v[0]*p->r[2] - p->r[0]*p->v[2];
		hz = p->r[0]*p->v[1] - p->v[0]*p->r[1];
		h = sqrt(hx*hx + hy*hy + hz*hz);
		
		// Work out maximum inclination
		incl = acos(hz/h);
		sin_incl2 = SQ(sin(incl));
		
		// If inclination is greater than InclMax
		// Which is set in ss.par then override
		
		if (sin_incl2 > InclMax)
			InclMax = sin_incl2;
	}
	return InclMax;
}

void pkdDustBinsGetMass(PKD pkd,DUST_BINS_PARAMS *DB,DustBins pDustBins[],
						double dt,double M,double pDustBinsMassLoss[])
{

    // Function calculates the total dust lost from all DustBins
    // Due to accreting particles within them
    
	PARTICLE *p;
	double a,x;
	int i,nLocal = pkdLocal(pkd);
	
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		//printf("x %e y %e z %e\n",p->r[0],p->r[1],p->r[2]);
		// For each particle find the DustBin associated
		collLocateBin(p,DB,&p->iBin,&a); /*DEBUG for now, "a" is "r_xy"!*/

		// If particle lies outside of DustBin domain
		// Then set dDustMass to zero and return

		if (p->iBin < 0 || p->iBin >= DB->nDustBins){
			p->dDustMass = 0.0;
		}

		// Otherwise find particle parameters and dust accretion

		else {
			double omega,P,h[3],h2,e,rho; 
			assert(a > 0.0);
			
			// Angular speed
			omega = sqrt((M + p->fMass)/CUBE(a));
			
			// Orbital period
			P = TWO_PI/omega;
			assert(P > 0.0);
			
			// Eccentricity from S.P.A.M
			h[0] = p->r[1]*p->v[2] - p->r[2]*p->v[1];
			h[1] = p->r[2]*p->v[0] - p->r[0]*p->v[2];
			h[2] = p->r[0]*p->v[1] - p->r[1]*p->v[0];
			h2 = SQ(h[0]) + SQ(h[1]) + SQ(h[2]);
			x = 1.0 - h2/(a*(M + p->fMass));
			if (fabs(x) < TINY)
				e = TINY;
			else
				e = sqrt(x);
			
			// If orbit parabolic or hyperbolic then set eccentricity to circular
			
			if (isnan(e)) {
			  e=0.0;
			}
			if (e >= 1.0 || e < 0.0) e = 0.0;
			assert(e < 1.0);
			
			// Calculate the amount of dust accreted onto particle
			assert(pDustBins[p->iBin].dVolume > 0.0);
			rho = pDustBins[p->iBin].dMass/pDustBins[p->iBin].dVolume;
			p->dDustMass = e*PI*SQ(RADIUS(p))*TWO_PI*a*rho*dt/P;
			
			
			// Add up total amount accreted and set this equal
			// To the amount of dust lost from DustBins
			
			pDustBinsMassLoss[p->iBin] += p->dDustMass;
			
		}
	}
}

#undef TINY


/*DEBUG following function fails pedantic compilation checking due to
  setting unused variables, which occurs because RUBBLE_FRICTION
  should really be a variable rather than a macro...*/

void pkdDustBinsApply(PKD pkd,double M,double pMassIncrFrac[],int nBins,double dt,
					  DustBins pDustBins[],double dVdisp
#ifdef ORIGIN_HISTOGRAM
					  ,FLOAT aOriginBins[][NUM_ORIGIN_BINS]
#endif
					  )
{
    // Applies dust to the particles during accretion by
    // Knowing the total mass of dust lost in pkdDustBinsGetMass()
    // Increases particle mass, radius and adjust velocity by
    // Momentum-conserving drag force application
    
	
	static int bFirstCall = 1;
		
	PARTICLE *p;
	double dm,R,omega,vkx,vky,diff_vx,diff_vy,diff_vz;
	double r,a,R_H,v2,v_mag,incl,sin_incl2,rho,v_kep2;
	double dvx,dvy,dvz,new_vx = 0.0,new_vy = 0.0,new_vz = 0.0;
	double h,hx,hy,hz,v2_diff;
	int i;
	
#ifndef COLLMOD_FRICTION
	assert(0);
#endif	

	if (bFirstCall) {
		printf("COLLMOD FRICTION OPTION = %i\n",COLLMOD_FRICTION); 
		printf("EXPANSION FACTOR = %g\n",EXPANSION_FACTOR);
	}


    // For each particle in the system
    
	for (i=0;i<pkd->nLocal;i++) {
	
		p = &pkd->pStore[i];
		if (p->iBin >= 0 && p->iBin < nBins) {
			
			// Work out mass to be accreted onto particle
			dm = pMassIncrFrac[p->iBin]*p->dDustMass;
			assert(dm >= 0.0);
			assert(p->fMass > dm);
			
			// Compute change to orbit based on momentum conservation
			
			// Angular speed of the dust in plane
			omega = sqrt(M)/pow(SQ(p->r[0]) + SQ(p->r[1]),0.75);
			
			// Keplarian velocity components
			vkx = - omega*p->r[1];
			vky =   omega*p->r[0];
			
			// If the planetesimals instantaneous velocity is larger
			// Than the keplarian velocity at that locations
			// Then planetesimal will experience a drag force
			// If the instantaneous velocity is smaller then
			// The planetesimal will receive a kick
			// Dust should circularize planetesimal orbit
			
			// Find difference in dust and planetesimal velocity
			diff_vx = p->v[0]-vkx;
			diff_vy = p->v[1]-vky;
			diff_vz = p->v[2];
			
			// Calculate the inclination
			hx = p->r[1]*p->v[2] - p->v[1]*p->r[2];
			hy = p->v[0]*p->r[2] - p->r[0]*p->v[2];
			hz = p->r[0]*p->v[1] - p->v[0]*p->r[1];
			h = sqrt(hx*hx + hy*hy + hz*hz);
			incl = acos(hz/h);
			sin_incl2 = SQ(sin(incl));
			
			r = sqrt(SQ(p->r[0])+SQ(p->r[1])+SQ(p->r[2]));
			v2 = SQ(p->v[0])+SQ(p->v[1])+SQ(p->v[2]);
			v2_diff = SQ(diff_vx)+SQ(diff_vy)+SQ(diff_vz);
			a = 1.0/(2.0/r - v2/(M + p->fMass));
			a = r;
			
			// Calculate the 'Hill Radius' which is the
			// planetesimals maximum radius of influence
			R_H = pow(p->fMass/3,(1.0/3.0))*a;
			
			assert(pDustBins[p->iBin].dVolume > 0.0);
			rho = pDustBins[p->iBin].dMass/pDustBins[p->iBin].dVolume;
			
			// Find the new velocities
 			new_vx = vkx + p->fMass*diff_vx/(dm + p->fMass);
			new_vy = vky + p->fMass*diff_vy/(dm + p->fMass);
			new_vz = p->fMass*diff_vz/(dm + p->fMass);
			
			v_kep2 = SQ(vkx) + SQ(vky);
			v_mag = sqrt(v2_diff);
			
#if (COLLMOD_FRICTION == COLLMOD_DYNAMICAL)
			
			sigma_2 = 2*dVdisp*v_kep2;
			
			b_max = R_H + a*sqrt(dVdisp + sin_incl2);
			lamda = b_max*(v2_diff+sigma_2)/p->fMass;
			coulomb = log(1 + SQ(lamda));
			
			// Assumes Max dist for background see BT p 424
			
			if (sin_incl2 >= dVdisp) { 
				vdynx = -TWO_PI*p->fMass/v2_diff*coulomb*rho*dt*diff_vx/v_mag*EXPANSION_FACTOR*EXPANSION_FACTOR;
				vdyny = -TWO_PI*p->fMass/v2_diff*coulomb*rho*dt*diff_vy/v_mag*EXPANSION_FACTOR*EXPANSION_FACTOR;
				vdynz = -TWO_PI*p->fMass/v2_diff*coulomb*rho*dt*diff_vz/v_mag*EXPANSION_FACTOR*EXPANSION_FACTOR;
			} 
			else {
				vdynx = -2.0/3*dt*sqrt(TWO_PI)/sqrt(sigma_2*sigma_2*sigma_2)*coulomb*rho*p->fMass*diff_vx*EXPANSION_FACTOR*EXPANSION_FACTOR;
				vdyny = -2.0/3*dt*sqrt(TWO_PI)/sqrt(sigma_2*sigma_2*sigma_2)*coulomb*rho*p->fMass*diff_vy*EXPANSION_FACTOR*EXPANSION_FACTOR;
				vdynz = -2.0/3*dt*sqrt(TWO_PI)/sqrt(sigma_2*sigma_2*sigma_2)*coulomb*rho*p->fMass*diff_vz*EXPANSION_FACTOR*EXPANSION_FACTOR;
			}
			
			// Assign new velocities corrected for dynamical friction
			new_vx += vdynx;
			new_vy += vdyny;
			new_vz += vdynz;
			
#endif
			
#if (COLLMOD_FRICTION == COLLMOD_GAS || COLLMOD_FRICTION == COLLMOD_DYNANDGAS)

// Convert 1700 in gcc^-1 converted to pkdgrav units
#define SIGMA_G (1700*1.495e13*1.495e13/1.9891e33) 
// Gamma = 1.4 (for H_2), T = 400 K, m=2m_H converted to pkdgrav units
#define C_S (sqrt(1.4*1.38e-23*400/(2*1.67e-27))/30000.0) 

			{
				double omega,e,tau_damp,k;
				omega = sqrt((M + p->fMass)/(a*a*a));
				e = h*h/(a*(M + p->fMass));
				if (e > 1.0 + 1.0e-6) {
					(void) fprintf(stderr,"h^2/aGM = %g > 1.000001\n",e);
					(void) fprintf(stderr,"pos = %g %g %g\n",p->r[0],p->r[1],p->r[2]);
					(void) fprintf(stderr,"vel = %g %g %g\n",p->v[0],p->v[1],p->v[2]);
					(void) fprintf(stderr,"Kep = %g %g 0\n",vkx,vky);
					(void) fprintf(stderr,"Info: iOrder=%i, mass=%g, a=%g, h=%g, M=%g\n",p->iOrder,p->fMass,a,h,M);
					assert(0);
				}
				
				// Circularise hyperbolic orbits
				
				if (e > 1.0)
					e = 0.0;
				else
					e = sqrt(1.0 - e);
				tau_damp = (M/p->fMass)*(M/(SIGMA_G*a*a))*pow(C_S/(a*omega),4.0)*(1.0/omega)*(1.0/(EXPANSION_FACTOR*EXPANSION_FACTOR));
				
				if (bFirstCall)
					(void) printf("RUBBLE GAS ECCENTRICITY DAMPING TIMESCALE = %g yrs\n",tau_damp/TWO_PI);
					
				// The following assumes all velocity components damp equally and that the eccentricity is small (so that e ~ (v-vk)/vk)
				k = - e*(dt/tau_damp)*v_mag*sqrt(v_kep2)/(diff_vx*p->v[0] + diff_vy*p->v[1] + diff_vz*p->v[2]);
				new_vx += k*p->v[0];
				new_vy += k*p->v[1];
				new_vz += k*p->v[2];
			}
#endif
			
			// This ensures corrections are no larger than amount
			// required to drive to zero eccentricity/inclination
			// (we're assuming here we can do this one component at a
			// time!) (note, sign could go either way: planetesimal
			// may be going faster or slower than Keplerian)
			
			dvx = new_vx - p->v[0];
			dvy = new_vy - p->v[1];
			dvz = new_vz - p->v[2];
			
			if (fabs(dvx) > fabs(diff_vx))
				dvx = (dvx/fabs(dvx))*fabs(diff_vx);
			
			if (fabs(dvy) > fabs(diff_vy))
				dvy = (dvy/fabs(dvy))*fabs(diff_vy);
			
			if (fabs(dvz) > fabs(diff_vz))
				dvz = (dvz/fabs(dvz))*fabs(diff_vz);
			
			// Apply corrections
			p->v[0] += dvx;
			p->v[1] += dvy;
			p->v[2] += dvz;
			
			// Grow particle radius
			R = RADIUS(p)*pow(1.0 + dm/p->fMass,1.0/3.0);
			p->fSoft = R;
			
#ifdef ORIGIN_HISTOGRAM
			MergeHistograms(p->fMass, p->origin_bins, dm, aOriginBins[p->iBin]);
#endif 

			// Grow particle mass
			
			// Add accretion mass to particle mass
			p->fMass += dm;
			
			// Make sure the dust associated with particle returns to zero
			p->dDustMass = 0.0;
			
		} 
		else {
			;
		}
	}
	if (bFirstCall)
		bFirstCall = 0;
}

void 
pkdCollModResetColFlag(PKD pkd)
{
	// Resets the flag that calls collmod.c or
	// collision.c to collide particles together

	int i,nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;i++)
		pkd->pStore[i].bMayCollide = 0;
}


int 
pkdCollModCheckForKDKRestart(PKD pkd)
{

    // Checks to see if particle flagged for collision

	int i,nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;i++)
		if (pkd->pStore[i].bMayCollide)
			return 1;
	
	return 0;
}

void 
pkdCollModStep(PKD pkd,double dMaxStep,double dMinStep)
{

	// Particles can find themselves on two differing 'rungs'
	// Which determines the timestep applied to it.
	// Those NOT undergoing collision need to be on the top 
	// Normal rung. Those undergoing a collision event will
	// Be assigned to the lower rung of small timestep
	// ADDED: To keep rungs in sync (since collisions are determined
	// *during* the drift), planetesimals that *may* collide during
	// this interval should also be forced to the lowest rung.
	
	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	// Search through all the particles
	
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		
		// If particle is a planetesimal and is not assigned
		// A collision flag then put it in top rung
		
		if (p->bMayCollide == 1){
			pkd->pStore[i].dt = dMinStep;
			//printf("Are we taking things slowly?\n");
			}
		
		// Otherwise assign particle to lower rung
		
		else {
			pkd->pStore[i].dt = dMaxStep;
			//printf("Are we taking things fastly?\n");
	    }
	}
}

void collCreatePart(const COLLIDER *col_a,const COLLIDER *col_b,
					COLLIDER **col_c)
{
	
	// Function is called when collCollModCollide() determines
	// That a collision results in a perfect merge
	// Thus combines two colliding particles into one
	
	VECTOR r_temp_a,r_temp_b,v_temp_a,v_temp_b,rcom,vcom;
	double orig_den;
	
	// Create outcome particle named col_c
	*col_c = (COLLIDER *) malloc(sizeof(COLLIDER));
	assert(*col_c != NULL);
	
	// Assign mass of col_c as col_a + col_b
	(*col_c)->fMass = col_a->fMass + col_b->fMass;
	
	// Find the density of the colliders under the assumption
	// That particles A and B have the same value
	orig_den = col_a->fMass/(4.0/3.0*PI*pow(col_a->fRadius,3.0));
	
	// Assign col_c radius based on this density and the mass
	(*col_c)->fRadius = pow(((*col_c)->fMass/(4.0/3.0*PI*orig_den)),(1.0/3.0));
	
	// Determines the system centre of mass location
	COPY_VEC(col_a->r,r_temp_a);
	SCALE_VEC(r_temp_a,col_a->fMass);
	COPY_VEC(col_b->r,r_temp_b);
	SCALE_VEC(r_temp_b,col_b->fMass);
	ADD_VEC(r_temp_a,r_temp_b,rcom);
	COPY_VEC(rcom,(*col_c)->r);
	NORM_VEC((*col_c)->r,col_a->fMass + col_b->fMass);
	
	// Calculates the outcome velocity for col_c
	COPY_VEC(col_a->v,v_temp_a);
	SCALE_VEC(v_temp_a,col_a->fMass);
	COPY_VEC(col_b->v,v_temp_b);
	SCALE_VEC(v_temp_b,col_b->fMass);
	ADD_VEC(v_temp_a,v_temp_b,vcom);
	COPY_VEC(vcom,(*col_c)->v);
	NORM_VEC((*col_c)->v,col_a->fMass + col_b->fMass);
	
	// Collider structure contains acceleration
	// But this is not used so value is zeroed
	ZERO_VEC((*col_c)->a);
	
	// No spin information so zeroed
	// *IMPROVE THIS LATER*
	ZERO_VEC((*col_c)->w);
	
	// Set col_c to be on top rung
	
	if (col_a->iRung >= col_b->iRung)
		(*col_c)->iRung = col_a->iRung;	
	else
		(*col_c)->iRung = col_b->iRung;
		
	(*col_c)->bFrag = 0;
		
	//printf("Merge Outcome: mass %e, radius %e, position %e %e %e\n", (*col_c)->fMass, (*col_c)->fRadius, (*col_c)->r[0], (*col_c)->r[1], (*col_c)->r[2]);

#ifdef ORIGIN_HISTOGRAM
	{
		int i;
		for(i = 0; i < NUM_ORIGIN_BINS; ++i)
			(*col_c)->origin_bins[i] = col_a->origin_bins[i];
		MergeHistograms(col_a->fMass, (*col_c)->origin_bins, col_b->fMass, col_b->origin_bins);
	}
#endif
}

static void
collRanVel(double speed_max,VECTOR v)
{
    // Calculates a random velocity
    // Used for offsets in fragmentation
    
	double v_max = speed_max/sqrt(N_DIM);

	do {
	  SET_VEC(v,(randUniform()*2.-1),(randUniform()*2-1),(randUniform()*2-1));
		SCALE_VEC(v,v_max);
		//printf("v = %e %e %e\n",v[0],v[1],v[2]);
	} while (MAG_SQ(v) > SQ(speed_max));
	
}

//JD: This is here for archiving and reference purposes, delete when new version works
void collCreateColliders_OLD(int *nOut, double b, double C, double M_lr, double M_slr,
							double rhoint, const COLLIDER *targ, const COLLIDER *proj,
							double *M_tail, COLLIDER *col_c, PKD pkd, double alpha)
{

	// Called when a collision results in fragmentation
	// Produces the fragments based on power law tail
	
    VECTOR r_temp_a, r_temp_b, rcom, p_init, p_lr, p_rem, p_i,pos_prev,pos_temp;
    VECTOR v_temp_a, v_temp_b, vcom,v_unit,v_ran_offset;
    double D_prev, D_N, mag_p_unit,v_esc,part_max_speed;
    double r_com, mr, mtot, rho;
	double beta = 2.85;
	int i,j,k,s;

	rho = rhoint;
	
	// Makes largest remnant (LR) using mass calculated in collCollModCollide()
	col_c[0].fMass = (1/MASS_CGS)*M_lr;
	col_c[0].fRadius = (1/LEN_CGS)*pow((M_lr/(4./3.*PI*rho)),1./3.);

	// Calculates the centre of mass and puts the LR there
	COPY_VEC(targ->r,r_temp_a);
	SCALE_VEC(r_temp_a,targ->fMass);
	COPY_VEC(proj->r,r_temp_b);
	SCALE_VEC(r_temp_b,proj->fMass);
	ADD_VEC(r_temp_a,r_temp_b,rcom);
	COPY_VEC(rcom,col_c[0].r);
	NORM_VEC(col_c[0].r,targ->fMass + proj->fMass);

	// Calculates centre of mass velocity 
	COPY_VEC(targ->v,v_temp_a);
	SCALE_VEC(v_temp_a,targ->fMass);
	COPY_VEC(proj->v,v_temp_b);
	SCALE_VEC(v_temp_b,proj->fMass);
	ADD_VEC(v_temp_a,v_temp_b,vcom);
	
	COPY_VEC(vcom,p_init);
	NORM_VEC(vcom,targ->fMass+proj->fMass);

	printf("Impact Parameter is %e\n", b);

	// Assigns velocity of LR based on impact parameter value
	
	printf("COM Velocity %e %e %e\n",vcom[0],vcom[1],vcom[2]);
	
	if (b == 0) {
		// If head on then take COM velocity
		COPY_VEC(vcom, col_c[0].v);
	} else if (b > 0.7) {
		// If high impact angle then LR should have velocity of target
		COPY_VEC(targ->v, col_c[0].v);
	} else {
		// An intermediate impact angle value means velocity 
		// Is calculated using a quasi-linear function of b
		// At the moment assume purely linear dependance
		// Need some information about the direction of the velocity vector 
		// Just do it with components? *IMPROVE THIS LATER*
		
	  col_c[0].v[0] = 1/0.7*(targ->v[0]*b + vcom[0]*(0.7-b));
	  col_c[0].v[1] = 1/0.7*(targ->v[1]*b + vcom[1]*(0.7-b));
	  col_c[0].v[2] = 1/0.7*(targ->v[2]*b + vcom[2]*(0.7-b));
	}
	
	printf("Largest Remnant Velocity %e %e %e\n", col_c[0].v[0], col_c[0].v[1], col_c[0].v[2]);

	// Zero the particle spin and acceleration
	ZERO_VEC(col_c[0].w);
	ZERO_VEC(col_c[0].a);
	
	// A fragmentation event may result in only one resolved particle
	// With the rest going into dust due to the resolution limit
	// If there IS more than one particle then complete this loop
	
	if (*nOut > 1) {
		// ZML CHECKS
		//double pos, semi, h_2, speed_2, ecc;	
		// Set mass and radius components of the second largest remnant (SLR)
		col_c[1].fMass = (1/MASS_CGS)*M_slr;
		col_c[1].fRadius = (1/LEN_CGS)*pow((M_slr/(4./3.*PI*rho)),1./3.);
		
		// Initialise the tail mass
		*M_tail = 0.0;
		
		// If a tail is produce (more than two remnants) then complete loop
		
		if (*nOut > 2) {
			int ii;
			
			// Set D_prev the distance to previous remnant
			D_prev = col_c[1].fRadius*LEN_CGS; //CGS
			printf("Before loop M_tail = %e\n", *M_tail);
			printf("C %e, beta %e M_slr %e D_prev/R_slr %e %e, M_lr %e R_lr %e\n", C, beta, M_slr, D_prev, (LEN_CGS*col_c[1].fRadius), M_lr, (LEN_CGS*col_c[0].fRadius));  
			
			// For all remnants in tail
			
			for (ii=2;ii<*nOut;ii++) {
				// D_N in CGS
				D_N = pow((pow(D_prev,-beta) + beta/C),-(1.0/beta)); 
			  
				// D_N is a scaling factor that determines the mass and radius
				// Of all the tail remnants
				col_c[ii].fRadius = (1/LEN_CGS)*D_N;
				col_c[ii].fMass = (1/MASS_CGS)*(4./3.*PI*pow(D_N,3)*rho);
			  
				// Find the tail mass by summing over all remnants masses
				*M_tail += col_c[ii].fMass;
				D_prev = D_N; //CGS  
			}

			printf("M_tail = %e\n", (MASS_CGS*(*M_tail)));
		}
		
		// Determine momentum, velocity and location parameters for tail remnants
		// Assume velocity of remnants is mass-weighted so give all remnants same velocity
		COPY_VEC(col_c[0].v,p_lr); //SYS
		printf("col_c[0].v %e %e %e\n",col_c[0].v[0],col_c[0].v[1],col_c[0].v[2]);
		SCALE_VEC(p_lr,(1/MASS_CGS*M_lr));
		
		// Calculate the remaining linear momentum
		SUB_VEC(p_init,p_lr,p_rem); //SYS
		printf("p_init %e %e %e\n",p_init[0],p_init[1],p_init[2]);
		printf("p_lr %e %e %e\n",p_lr[0],p_lr[1],p_lr[2]);

		// Then distribute this over tail remnants + SLR
		COPY_VEC(p_rem, p_i); 
		NORM_VEC(p_i,(*nOut-1));
		
		// Label remnants for colouring in pkdPutColliderInfo
		for (j=1;j<*nOut; j++) {
			col_c[j].bFrag = 1;
		}

		// For tail remants + SLR
		
		for (j=1; j<*nOut; j++) {
		  // Assign velocity as the number weighted momentum
		  // Then normalise by the remnant mass
		  COPY_VEC(p_i,col_c[j].v);
		  NORM_VEC(col_c[j].v,col_c[j].fMass);
		  
		  //printf("Particle %i: vx %e vy %e\n", j, col_c[j].v[0],col_c[j].v[1]);
		}

		COPY_VEC(col_c[0].r,pos_prev);
		mag_p_unit = MAG(p_rem);
		COPY_VEC(p_rem,v_unit);
		printf("p_rem %e %e %e\n",p_rem[0],p_rem[1],p_rem[2]);
		NORM_VEC(v_unit,mag_p_unit);
		printf("v_unit %e %e %e\n",v_unit[0],v_unit[1],v_unit[2]);
		
		//v_slr = MAG(col_c[1].v);

		mr=0.0;
		mtot=0.0;
		
		// For all tail remnants + SLR
		
		for (k=1; k<*nOut; k++) {
			
			// Find the previous particle position
			// Then add a scaled previous particle radius to this
			// Assign next particle location as this
			COPY_VEC(v_unit,pos_temp);
			SCALE_VEC(pos_temp,4.0*col_c[k-1].fRadius);
			ADD_VEC(pos_prev,pos_temp,col_c[k].r);
			
			// For remaining tail remnants
			
			if (k>0) {
			  for (s=0; s<k; s++) {
			    mr+=col_c[s].fMass*col_c[s].r[0];
			    mtot+=col_c[s].fMass;
			  }
			  
			  // Find centre of mass based on all previous particles
			  r_com = col_c[k].r[0]-(mr/mtot);
			  //printf("mr %e\n",mr);
			  //printf("col_c[k].r[0] %e\n",col_c[k].r[0]);
			  
			  if (r_com < 0 ) {
			    r_com = r_com * -1.0;
			  }
			  
			  // Find escape velocity for particle based on previous masses
			  v_esc = sqrt(2.*col_c[0].fMass/r_com);
			  //printf("Particle %d V_esc %e\n",k,v_esc);
			  //printf("v_esc r_com %e %e\n",v_esc,r_com);
			  
			  // Run away from your father
			  col_c[k].v[0] = col_c[k].v[0]+col_c[0].v[0];
			  col_c[k].v[1] = col_c[k].v[1]+col_c[0].v[1];
			  col_c[k].v[2] = col_c[k].v[2]+col_c[0].v[2];
			}
			
			part_max_speed = 0.5*MAG(col_c[1].v);
			collRanVel(part_max_speed,v_ran_offset);
			SCALE_VEC(v_ran_offset,5.0); 
			ADD_VEC(col_c[j].v,v_ran_offset,col_c[j].v);
			
			for (s=0; s<3; s++) {
			  col_c[k].v[s]=v_esc*(rand()%(1-(-1))+(-1))+(col_c[k].v[s]);
			}
			
			// Assign particle location as the previous particle location for next step
			COPY_VEC(col_c[k].r,pos_prev);
			
			
			//printf("p %d mass %e v %.25e ev %e and r %e %e %e\n",k,col_c[k].fMass,col_c[k].v[0],v_esc,col_c[k].r[0],col_c[k].r[1],col_c[k].r[2]);
			
			// Zero spin and acceleration
			ZERO_VEC(col_c[k].a);
			ZERO_VEC(col_c[k].w);
		}
	}
	
#ifdef ORIGIN_HISTOGRAM
		for (j=0;j<*nOut;j++) {
			for(i = 0; i < NUM_ORIGIN_BINS; ++i) {
				col_c[j].origin_bins[i] = targ->origin_bins[i];
			}
			MergeHistograms(targ->fMass, col_c[j].origin_bins, proj->fMass, proj->origin_bins);
		}
#endif

}

double collModulo(double val, double period){
	//this should calculate val%period for doubles and always have sign = sign of period
	double wrap;
	wrap = fmod(val, period); // has sign of wrap
	if(wrap*period < 0){ //if wrap and period are not the same sign
		wrap += period; //add wrap to the period
	}
	return(wrap);
}

double collBoxMuller(double mu, double sigma){
	return mu + sigma * randGaussian();
}

double collPGauss(double mu, double sigma, double period){
	//Creates a periodic gaussian distribution 
	return(collModulo(collBoxMuller(mu, sigma), period));
}



double collGenTheta_NEW(double theta1, double theta2, double theta_s1, double theta_s2, double A1, double A2, int N){
	static int i=0;
	static int j=0;
	double theta;
	//printf("theta1: %G, theta2: %G\n", theta1, theta2);
	if(j >= (A1*N/(A1+A2)) ){ // if have enough theta2's
		//generate stuff for theta1
		theta = collPGauss(theta1, theta_s1, 360.0);
		i++;
	}
	if( i >= (A2*N/(A1+A2)) ){ //if have enough theta1's
		//generate stuff for theta2
		theta = collPGauss(theta2, theta_s2, 360.0);
		j++;
	}
	if( randUniform() < (A1/(A1+A2)) ){ //correct proportion of the time do this
		theta = collPGauss(theta1, theta_s1, 360.0);
		i++;
	}
	else{ //rest of time do this
		theta = collPGauss(theta2, theta_s2, 360.0);
		j++;
	}
	i++;
	return(theta);
}

double collGenPhi_NEW(double theta, double phi1, double phi2, double phi_s1, double phi_s2, double theta1, double theta2){
	double phi;
	if(theta < 180.0){
		//printf("THETA < 180\n");
		phi = collModulo((collPGauss(phi1+90.0, phi_s1, 180.0)-90.0), 360.0); 
	}
	else{
		//printf("THETA > 180\n");
		phi = (collPGauss(phi2-90.0, phi_s2, 180.0)+90.0);
	}
	return(phi);
}

double collGenTheta(double theta1, double theta2, double theta_s1, double theta_s2, double A1, double A2, int N){
	static int i=0;
	//change to not depend on order
	double theta;
	//these two should make sure that theta1 is associated with phi1, theta2 associated with phi2
	if(i < (A1*N/(A1+A2)) ){
		//generate stuff for theta1
		theta = collPGauss(theta1, theta_s1, 360.0);
	}
	else{
		//generate stuff for theta2
		theta = collPGauss(theta2, theta_s2, 360.0);
	}
	i++;
	return(theta);
}

double collGenPhi(double phi1, double phi2, double phi_s1, double phi_s2, double A1, double A2, int N){
	static int i=0;
	double phi;
	//change to not depend on order
	if(i < (A1*N/(A1+A2)) ){
		phi = collModulo((collPGauss(phi1+90.0, phi_s1, 180.0)-90.0), 360.0); //phi1 = 0.0 should only exist between 0->90, 270->369
	}
	else{
		phi = (collPGauss(phi2-90.0, phi_s2, 180.0)+90.0); //phi2 = 180.0 should only exist between 90->270
	}
	i+=1;
	return(phi);
}

void collThetaPhi2Cart(VECTOR out_v, double theta, double phi, VECTOR c, VECTOR d, VECTOR e){
	// c, d, e are the new unit co-ord 'target' system vectors in the old cartesian co-ord system
	// can get c, d, e components easily from theta and phi, can go from c, d, e to x, y, z too
	// theta is angle in c-d plane, phi is angle in d-e plane
	// so, tan(theta) = v_d/v_c, tan(phi) = v_d/v_e
	// only want direction not magnitude so let |v| = 1
	// can set mag at the end, start off assuming that |d_coeff| = 1
	double d_c, e_d, c_coeff, d_coeff, e_coeff;
	VECTOR temp;
	//printf("theta: %G, phi: %G\n", theta, phi);
	d_c = tan(theta*PI/180.0); //remember C trig functions work in radians
	e_d = tan(phi*PI/180.0);
	//if theta <180 d_coeff is +ve, else is -ve, can fix magnitude and normalise later
	if(theta<180){
		d_coeff = 1;
	}
	else{
		d_coeff = -1;
	}
	c_coeff = d_coeff/d_c; // gets c-coefficient
	e_coeff = d_coeff*e_d; // gets e-coefficient
	//printf("c: %G, %G, %G\n", c[0], c[1], c[2]); //JD: Print out vectors for debugging
	//printf("d: %G, %G, %G\n", d[0], d[1], d[2]);
	//printf("e: %G, %G, %G\n", e[0], e[1], e[2]);
	//printf("c_coeff %G, d_coeff %G, e_coeff %G\n", c_coeff, d_coeff, e_coeff);
	// v_total = v1 + v2 + v3 + v4 .... 
	// we want to find where adding the correct proportions of c,d,e get us in x,y,z space
	CSCALE_VEC(c, c_coeff, temp);
	//printf("c temp %G, %G, %G\n", temp[0], temp[1], temp[2]);
	out_v[0] = temp[0];
	out_v[1] = temp[1];
	out_v[2] = temp[2];
	CSCALE_VEC(d, d_coeff, temp);
	//printf("d temp %G, %G, %G\n", temp[0], temp[1], temp[2]);
	out_v[0] += temp[0];
	out_v[1] += temp[1];
	out_v[2] += temp[2];
	CSCALE_VEC(e, e_coeff, temp);
	//printf("e temp %G, %G, %G\n", temp[0], temp[1], temp[2]);
	out_v[0] += temp[0];
	out_v[1] += temp[1];
	out_v[2] += temp[2];
	//printf("out_v %G, %G, %G\n", out_v[0], out_v[1], out_v[2]);
	NORM_VEC(out_v, MAG(out_v)); //make vector unit length

	return;
}

int collBounce(const COLLIDER *c1,const COLLIDER *c2,
			  double dEpsN,double dEpsT,
			  COLLIDER **cOut)
{
	/* bounces colliders, preserving particle order */
	//JD: this is a copy of pkdBounce, modified slightly for use here
	COLLIDER *co1,*co2;
	FLOAT n[3],s1[3],s2[3],v[3],s[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;
	
	//cOut should already be allocated
	assert(*cOut != NULL);

	(*cOut)[0] = *c1;
	(*cOut)[1] = *c2;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	mu = m1*m2/m;
	alpha = 2.5*(1.0/m1 + 1.0/m2);
	beta = 1.0/(1.0 + alpha*mu); /* always 2/7 for uniform spheres */

	a = 0;
	for (i=0;i<3;i++) {
		n[i] = (c2->r[i] - c1->r[i]);
		a += n[i]*n[i];
		}
	a = 1.0/sqrt(a);
	for (i=0;i<3;i++)
		n[i] *= a;

	s1[0] = r1*(c1->w[1]*n[2] - c1->w[2]*n[1]);
	s1[1] = r1*(c1->w[2]*n[0] - c1->w[0]*n[2]);
	s1[2] = r1*(c1->w[0]*n[1] - c1->w[1]*n[0]);

	s2[0] = r2*(c2->w[2]*n[1] - c2->w[1]*n[2]);
	s2[1] = r2*(c2->w[0]*n[2] - c2->w[2]*n[0]);
	s2[2] = r2*(c2->w[1]*n[0] - c2->w[0]*n[1]);

	for (i=0;i<3;i++) {
		v[i] = c2->v[i] - c1->v[i];
		s[i] = s2[i] - s1[i];
		}

	for (i=0;i<3;i++)
		u[i] = v[i] + s[i];

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];

	if (a >= 0) {
		printf("WARNING: Near miss a %G\n", a);
		return(1);
	}
	for (i=0;i<3;i++) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

	a = (1.0 + dEpsN);
	b = beta*(1.0 - dEpsT);
	for (i=0;i<3;i++)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;

	co1 = &(*cOut)[0];
	co2 = &(*cOut)[1];

	for (i=0;i<3;i++) {
		co1->v[i] += a*p[i];
		co2->v[i] += b*p[i];
		co1->w[i] += c*q[i];
		co2->w[i] += d*q[i];
		}

	return(0);
	}

/*DEBUG following function fails pedantic compilation checking due to
  too many unused variables...*/

void collCreateColliders(int *nOut, double b, double C, double M_lr, double M_slr, double M_tlr,
							double rhoint, const COLLIDER *targ, const COLLIDER *proj,
							double *M_tail, COLLIDER **col_c1, PKD pkd, double alpha, double V_esc,
							int bsupercat)
{
	// Called when a collision results in fragmentation
	// Produces the fragments based on power law tail
	
    VECTOR r_temp_a, r_temp_b, rcom, p_init, p_lr, p_rem, p_i,pos_prev,pos_temp;
    const COLLIDER *tempc;
	COLLIDER *col_c;
    VECTOR v_temp_a, v_temp_b, vcom,v_unit,v_ran_offset;
    VECTOR c, d, e, targ_v_com, v_coll; //JD: Vectors for debris model
    double mass_r, b_crit, A1, A2, theta1, theta2, phi1, phi2, theta_s1, theta_s2, phi_s2, phi_s1;
    double D_prev, D_N, mag_p_unit,v_esc,v_slr,part_max_speed, E_in, E_lr, E_rem;
    double r_com, mr, mtot, rho;
	double beta = 2.85;
	double M_bin, M_avail, M_not_tail, S, A, v, dv;
	int i,j, rem_start=1, rem_i = 1;
	
	puts("HERE2");
	
	//JD: DEBUG
	col_c = malloc((*nOut)*sizeof(COLLIDER));
	assert(col_c !=NULL);
	for (j=0;j<*nOut; j++) {
			col_c[j].bFrag = 0;
		}
	if(targ->fMass < proj->fMass){ //check masss
		tempc = targ;
		targ = proj;
		proj = tempc;
	}
	//CHECKING PKD PROPERTIES
	//int idSelf;
	//int nThreads;
	printf("nThreads: %d, idSelf: %d\n", pkd->nThreads, pkd->idSelf);
	//checking impact positions
	printf("targ->r: %G, targ->v: %G\n", MAG(targ->r), MAG(targ->v));
	printf("proj->r: %G, proj->v: %G\n", MAG(proj->r), MAG(proj->v));
	
	srand(time(NULL)); // seed random number generator for rand()
	rho = rhoint;
	SUB_VEC(proj->r, targ->r, c); //JD: Get centre joining vector (target to projectile)
	//JD: Ok, so we want to define theta as being in the plane that is defined by 'c' and
	//JD: the velocity of the target in the centre of mass frame
	// Makes largest remnant (LR) using mass calculated in collCollModCollide()
	col_c[0].fMass = (1./MASS_CGS)*M_lr;
	col_c[0].fRadius = (1./LEN_CGS)*pow((M_lr/(4./3.*PI*rho)),1./3.);
	ADD_VEC(targ->v, proj->v, v_coll); // find the collision velocity
	
	// Calculates the centre of mass and puts the LR there
	COPY_VEC(targ->r,r_temp_a);
	SCALE_VEC(r_temp_a,targ->fMass);
	COPY_VEC(proj->r,r_temp_b);
	SCALE_VEC(r_temp_b,proj->fMass);
	ADD_VEC(r_temp_a,r_temp_b,rcom);
	COPY_VEC(rcom,col_c[0].r); // Note rcom is actually position*mass ...
	NORM_VEC(col_c[0].r,targ->fMass + proj->fMass);
	//ZML fix
	NORM_VEC(rcom,targ->fMass + proj->fMass);

	// Calculates centre of mass velocity 
	COPY_VEC(targ->v,v_temp_a);
	SCALE_VEC(v_temp_a,targ->fMass);
	COPY_VEC(proj->v,v_temp_b);
	SCALE_VEC(v_temp_b,proj->fMass);
	ADD_VEC(v_temp_a,v_temp_b,vcom);
	
	//JD: Find the target velocity in the COM frame
	SUB_VEC(targ->v, vcom, targ_v_com);//JD: $1 - $2 = $3
	//JD: Axis that theta rotates around, e, is the vector that defines the c-targ_v_com plane
	CROSS(c, targ_v_com, e); //JD: $1 X $2 = $3 //JD: This might be backwards, 
	//JD: remember to try and swap c, e if broken
	//JD: remember, c defines the phi-rotation axis, the +ve c-targ_v_com plane defines the zero-point
	//JD: d-vector is defined by; c dot d = 0, d dot e = 0, c cross e = d (in unit vector case)
	CROSS(c, e, d); 
	NORM_VEC(c, MAG(c)); //JD: Make sure everything is a unit vector
	NORM_VEC(d, MAG(d));
	NORM_VEC(e, MAG(e));
	printf("c: %G, %G, %G\n", c[0], c[1], c[2]); //JD: Print out vectors for debugging
	printf("d: %G, %G, %G\n", d[0], d[1], d[2]);
	printf("e: %G, %G, %G\n", e[0], e[1], e[2]);
	mass_r = proj->fMass/targ->fMass; 
	mr = mass_r;
	//JD: OK, by this point we have defined the co-ord system, target velocity in COM frame, and mass ratio
	
	COPY_VEC(vcom,p_init);
	NORM_VEC(vcom,targ->fMass+proj->fMass);

	printf("Impact Parameter is %e\n", b);

	// Assigns velocity of LR based on impact parameter value
	
	printf("COM Velocity %e %e %e\n",vcom[0],vcom[1],vcom[2]);
	E_in = 0.5*(MAG_SQ(targ->v)*targ->fMass + MAG_SQ(proj->v)*proj->fMass); //energy into system
	//JD: Leave this part the same, LR handled by fit to data
	if (b == 0 || bsupercat || *nOut==1 || (col_c[0].fMass > targ->fMass)) {
		// If head on then take COM velocity
		COPY_VEC(vcom, col_c[0].v);
	} else if (b > 0.7) {
		// If high impact angle then LR should have velocity of target
		COPY_VEC(targ->v, col_c[0].v);
	} else {
		// An intermediate impact angle value means velocity 
		// Is calculated using a quasi-linear function of b
		// At the moment assume purely linear dependance
		// Need some information about the direction of the velocity vector 
		// Just do it with components? *IMPROVE THIS LATER*
		
	  col_c[0].v[0] = 1/0.7*(targ->v[0]*b + vcom[0]*(0.7-b));
	  col_c[0].v[1] = 1/0.7*(targ->v[1]*b + vcom[1]*(0.7-b));
	  col_c[0].v[2] = 1/0.7*(targ->v[2]*b + vcom[2]*(0.7-b));
	}
	
	printf("Largest Remnant Velocity %e %e %e\n", col_c[0].v[0], col_c[0].v[1], col_c[0].v[2]);
	printf("Target velocity: %G %G %G\n", targ->v[0], targ->v[1], targ->v[2]);
	// Zero the particle spin and acceleration
	ZERO_VEC(col_c[0].w);
	ZERO_VEC(col_c[0].a);
	//JD: By this point the largest remnant 'col_c[0]' is exactly where it should be
	
	
	// A fragmentation event may result in only one resolved particle
	// With the rest going into dust due to the resolution limit
	// If there IS more than one particle then complete this loop
	
	if (*nOut > 1) { //JD: more than just the largest remnant in tail
         // ZML CHECKS
                double pos, semi, h_2, speed_2, ecc;
	
		// Set mass and radius components of the second largest remnant (SLR)
		col_c[1].fMass = (1/MASS_CGS)*M_slr;
		col_c[1].fRadius = (1/LEN_CGS)*cbrt(M_slr/(4./3.*PI*rho));
		
		// Initialise the tail mass
		*M_tail = 0.0;
		
		// If a tail is produce (more than two remnants) then complete loop
		if (*nOut > 2 && M_tlr ==0) { //JD: more than lr and slr - creates the tail
			int ii;
			
			// Set D_prev the distance to previous remnant
			D_prev = col_c[1].fRadius*LEN_CGS; //CGS
			printf("Before loop M_tail = %e\n", *M_tail);
			printf("C %e, beta %e M_slr %e D_prev/R_slr %e %e, M_lr %e R_lr %e\n", C,
						beta, M_slr, D_prev, (LEN_CGS*col_c[1].fRadius), M_lr,
						(LEN_CGS*col_c[0].fRadius));  
			
			// For all remnants in tail
			
			for (ii=2;ii<*nOut;ii++) { //just tail not slr
				// D_N in CGS
			
				D_N = pow((pow(D_prev,-beta) + beta/C),-(1.0/beta)); 
			  
				// D_N is a scaling factor that determines the mass and radius
				// Of all the tail remnants
				col_c[ii].fRadius = (1/LEN_CGS)*D_N;
				col_c[ii].fMass = (1/MASS_CGS)*(4./3.*PI*pow(D_N,3)*rho);
			  
				// Find the tail mass by summing over all remnants masses
				*M_tail += col_c[ii].fMass;
				D_prev = D_N; //CGS  
			}

			printf("M_tail = %e\n", (MASS_CGS*(*M_tail)));
		} //JD: By this point we have mass+radius properties of all remnants: lr, slr, and tail.
		else if(*nOut ==3 && M_tlr != 0.0){
			col_c[2].fRadius = cbrt(3.*M_tlr/(4.*M_PI*rho))/LEN_CGS;
			col_c[2].fMass = M_tlr/MASS_CGS;
		}
		else if(*nOut > 3 && M_tlr !=0.0){
			int ii;
			//assert(*nOut > 3 && M_tlr !=0.0);
			col_c[2].fRadius = cbrt(3.*M_tlr/(4.*M_PI*rho))/LEN_CGS;
			col_c[2].fMass = M_tlr/MASS_CGS;
			D_prev = col_c[2].fRadius*LEN_CGS; //CGS
			for (ii=3;ii<*nOut;ii++) {//just tail not tlr
				// D_N in CGS
				D_N = pow((pow(D_prev,-beta) + beta/C),-(1.0/beta)); 			  
				// D_N is a scaling factor that determines the mass and radius
				// Of all the tail remnants
				col_c[ii].fRadius = (1/LEN_CGS)*D_N;
				col_c[ii].fMass = (1/MASS_CGS)*(4./3.*PI*pow(D_N,3)*rho);  
				// Find the tail mass by summing over all remnants masses
				*M_tail += col_c[ii].fMass;
				D_prev = D_N; //CGS  
			}
			printf("M_tail = %e\n", (MASS_CGS*(*M_tail)));
		}
		
		if(M_tlr != 0.0){
			int z;
			double sum;
			rem_start = 2; //is 1 by default, only changes when we have tlr
			//M_lr = M_slr; //the 'largest remnant' is now the second largest thing around
			//treat lr (target) and slr (projectile lr) as billiard balls
			//copy pkdBounce here!! probably don't need all of it.
			//make sure to do bounce then reduce mass, put extra momentum into fragments
			//make sure momentum is conserved
			printf("CHECK: targ.fMass %G proj.fMass %G\n", targ->fMass, proj->fMass);
			//VPRNT("CHECK: targ->v", targ->v);
			//VPRNT("CHECK: proj->v", proj->v);
			//VPRNT("CHECK: col_c[0]->v", col_c[0].v);
			//VPRNT("CHECK: col_c[1]->v", col_c[1].v);
			sum = 0;
			for(z=0;z<*nOut;z++){
				sum+=col_c[z].fMass;
			}
			printf("CHECK: Mass Before collBounce %G\n", sum);
			//printf("CHECK: targ.v %G %G %G\n", targ->v[0], targ->v[1], targ->v[2]);
			//printf("CHECK: proj.v %G %G %G\n", proj->v[0], proj->v[1], proj->v[2]);
			collBounce(targ, proj,
				0.8, 1.0, //dEpsN and dEpsT should ideally be the same as ss.par file...
				&col_c);
			printf("CHECK: targ.fMass %G proj.fMass %G\n", targ->fMass, proj->fMass);
			//VPRNT("CHECK: targ->v", targ->v);
			//VPRNT("CHECK: proj->v", proj->v);
			//VPRNT("CHECK: col_c[0]->v", col_c[0].v);
			//VPRNT("CHECK: col_c[1]->v", col_c[1].v);
			printf("CHECK: Before mass scaling col_c[0].fMass %G, col_c[1].fMass %G\n", col_c[0].fMass, col_c[1].fMass);
			printf("CHECK: Before mass scaling col_c[0].fRadius %G col_c[1].fRadius %G\n", col_c[0].fRadius, col_c[1].fRadius);
			col_c[0].fMass = M_lr/MASS_CGS; //should be = mass of target - ZML - shouldn't this be done in collBounce?
			col_c[1].fMass = M_slr/MASS_CGS; //shoul be < mass of projectile
			// ZML - the mass is reset but not the radius ... the target doesn't need the radius reset but in projectile disrupt M_slr is not M_proj //
			col_c[1].fRadius = cbrt(3.*M_slr/(4.*M_PI*rho))/LEN_CGS; // ZML added 02.07.14
			if(*M_tail !=0.0){ //if we have some mass in the tail
				COPY_VEC(col_c[1].v, v_temp_a);
				SCALE_VEC(v_temp_a, proj->fMass - (M_slr/MASS_CGS)); //change in momentum
				SCALE_VEC(v_temp_a, 1.0/(*nOut-rem_start)); //divide over tail remnants
				//col_c[1].fMass -= *M_tail; //take off the mass we lost
			} // a later step should re-distribute momentum correctly
			printf("CHECK: After mass scaling col_c[0].fMass %G col_c[1].fMass %G\n", col_c[0].fMass, col_c[1].fMass);
			printf("CHECK: After mass scaling col_c[0].fRadius %G col_c[1].fRadius %G\n", col_c[0].fRadius, col_c[1].fRadius);
			//else{ //else treat as a bouncing event
			//	*nOut = 2;
			//	return;
			//}
			sum = 0;
			for(z=0;z<*nOut;z++){
				sum+=col_c[z].fMass;
			}
			printf("CHECK: Mass after collBounce %G\n", sum);
			//printf("CHECK: targ.v %G %G %G\n", targ->v[0], targ->v[1], targ->v[2]);
			//printf("CHECK: proj.v %G %G %G\n", proj->v[0], proj->v[1], proj->v[2]);
			//M_lr = M_slr; //the 'largest remnant' is now the second largest thing around
			  
			
		}
		// Determine momentum, velocity and location parameters for tail remnants
		// Assume velocity of remnants is mass-weighted so give all remnants same velocity
		//JD: I shall use this to give all remnants the same magnitude of velocity, but not direction
		COPY_VEC(col_c[0].v,p_lr); //SYS units, copies velocity or LR to p_lr
		printf("col_c[0].v %e %e %e\n",col_c[0].v[0],col_c[0].v[1],col_c[0].v[2]);
		SCALE_VEC(p_lr,(1/MASS_CGS*M_lr)); //JD: turns velocity into momentum in CGS
		
		// Calculate the remaining linear momentum
		SUB_VEC(p_init,p_lr,p_rem); //SYS
		printf("p_init %e %e %e\n",p_init[0],p_init[1],p_init[2]); //show initial momentum
		printf("p_lr %e %e %e\n",p_lr[0],p_lr[1],p_lr[2]); //show LR-momentum
		printf("p_rem %G, %G, %G\n", p_rem[0], p_rem[1], p_rem[2]); //show remaining momentum
		/*
		// Then distribute this over tail remnants + SLR
		COPY_VEC(p_rem, p_i); 
		NORM_VEC(p_i,(*nOut-1));
		*/
		
		
		// Label remnants for colouring in pkdPutColliderInfo
		printf("col_c[0].bFrag %d\n", col_c[0].bFrag);
		for (j=rem_start;j<*nOut; j++) {
			col_c[j].bFrag = 1;
		}
		printf("col_c[0].bFrag %d\n", col_c[0].bFrag);
		/*
		// For tail remants + SLR
		
		
		for (j=1; j<*nOut; j++) {
		  // Assign velocity as the number weighted momentum
		  // Then normalise by the remnant mass
		  COPY_VEC(p_i,col_c[j].v);
		  NORM_VEC(col_c[j].v,col_c[j].fMass);
		  
		  //printf("Particle %i: vx %e vy %e\n", j, col_c[j].v[0],col_c[j].v[1]);
		} //JD: By this point, all remnants have had the same jolt to conserve momentum

		COPY_VEC(col_c[0].r,pos_prev); //JD: Copy position of lr
		mag_p_unit = MAG(p_rem); //JD: find magnitude of momentum adjustment
		COPY_VEC(p_rem,v_unit); //JD: copy momentum adjustment to v_unit
		printf("p_rem %e %e %e\n",p_rem[0],p_rem[1],p_rem[2]);
		NORM_VEC(v_unit,mag_p_unit); //JD: Normalise v_unit by magnitude of momentum adjustment
		printf("v_unit %e %e %e\n",v_unit[0],v_unit[1],v_unit[2]);
		//JD: v_unit is the unit vector in the direction of the velocity 
		
		v_slr = MAG(col_c[1].v); //JD: save magnitude of slr's velocity
		//JD: I overwrite this step later, may have to put back in to ensure momentum conservation afterwards
		
		mr=0.0;
		mtot=0.0;
		*/
		//JD: Put debris model here, need to decide how to fix velocity. slr is part of tail!
		//JD: Debris model may make previous momentum balancing step obsolete, can re-write at the end
		b_crit = 1.0/(cbrt(mass_r) +1.0); // critical impact parameter
		
		if(b < b_crit){ //find theta values
			theta1 = ((70.0*b_crit - 90.0)/b_crit)*b + 90.0;
		}
		else if(b >=b_crit){
			theta1 = 70.*b;
		}
		printf("theta_adj: %G\n", theta1+asin(b)*180/M_PI);
		theta2 = theta1 + 180.0;
		theta_s1 = 10.0;
		theta_s2 = 10.0;
		
		phi1 = 0.0; //find phi values
		phi2 = 180.0;
		phi_s1 = 145.0*alpha;
		phi_s2 = 145.0*alpha;
		//JD: MAKE 'v_norm' RELATE TO NORMALISED VELOCITY SO I CAN INCLUDE DEPENDENCE!
		printf("v_coll: %G, V_esc: %G\n", MAG(v_coll), V_esc); //test that we're using correct units (want CGS)
		if ((MAG(v_coll)/V_esc) < -1.4*(1.-b_crit)+0.85){ //low velocity therefore 1:1
			A1 = 1.0; //JD: Do this for now
			A2 = 1.0;
		}
		else{
			A1 = (1-b*(1-mr)); // if mr==1 then this should always be 1:1, how do I fix this? maybe |(1-b)/(mr-b)|?
			A2 = 1.0;
		}
		v_esc = sqrt(2.0*col_c[0].fMass/col_c[0].fRadius); //calculate escape velocity from lr
		
		printf("theta1 %G, theta2 %G, theta_s1 %G, theta_s2 %G, A1 %G, A2 %G\n", theta1, theta2,
			theta_s1, theta_s2, A1, A2);
		printf("phi1 %G, phi2 %G, phi_s1 %G, phi_s2 %G\n", phi1, phi2, phi_s1, phi_s2);
		
		E_lr = 0.5*(MAG_SQ(col_c[0].v)*col_c[0].fMass); // Energy of largest remnant
		E_rem = E_in - E_lr; //Energy available to remnants
		if(M_tlr!=0.0){
			E_rem = E_in - E_lr - 0.5*col_c[1].fMass*MAG_SQ(col_c[1].v);
		}
		assert(E_rem>=0.0); //Not allowed to have -ve energy
		E_rem = E_rem/(*nOut-rem_start); //energy available per remnant
		E_rem *= 1.0; //only assign a fraction of the available energy coefficient of restitution blah blah
		
		//OK, this is going to be a massive hassle. assume that we only have 'escaped' particles
		//in this simulation (that's what the collision model was for). Using the 'cumulative mass
		//in a velocity bin' histogram in Leinhardt12, you get this relation
		// log(mass in bin) = A - Sv //mass is in mass/total velocity in in escape velocities
		//A = -0.3M_lr/M_tot + 0.3 -> mass in lowest bin
		//S = 10^A/(ln(10)(dv)(M_rem/M_tot)) -> slope that conserves mass
		//ZML - do something simple ...
		//		M_not_tail = col_c[0].fMass;
		//M_avail = proj->fMass + targ->fMass;
		//if(M_tlr!=0.0){
		//	M_not_tail = col_c[0].fMass + col_c[1].fMass;
		//	M_avail = proj->fMass;
		//}
		//A = -0.3*col_c[rem_start-1].fMass/M_avail + 0.3; //largest chunk used here
		//dv = 2E-3;
		//S = pow(10, A)/(2.3026*dv*(targ->fMass + proj->fMass - M_not_tail)/(M_avail));
		//v = 1.0; //start at v = v_esc
		//M_bin = M_avail*pow(10, A - S*v);
		//printf("\nrem_start %d col_c[rem_start-1].mass %E, M_avail %E bFrag %d", rem_start, col_c[rem_start-1].fMass, M_avail, col_c[rem_start-1].bFrag);
		//printf("\nrem_start %d col_c[rem_start-1].mass/M_avail %E, A %E, S %E, M_bin %E, v %E, dv %E\n", rem_start, col_c[rem_start-1].fMass/M_avail, A, S, M_bin, v, dv);
		
		for (i=rem_start; i<*nOut; i++){ //For all tail remnants + SLR (or TLR)
			// give particles their velocities, they are ordered in decending mass so might have to randomise
			// get theta+phi randoms, are the particles in any specific order?
		  int vesc_i;
		  double theta, phi, temp_v_mag;
		  VECTOR temp_v, temp_v2, rel_vi,rel_vslr,d_offset;
		  theta = collGenTheta_NEW(theta1, theta2, theta_s1, theta_s2, A1, A2, *nOut);
		  //phi = collGenPhi(phi1, phi2, phi_s1, phi_s2, A1, A2, *nOut);
		  phi = collGenPhi_NEW(theta, phi1, phi2, phi_s1, phi_s2, theta1, theta2);
		  //printf("theta %G, phi %G\n", theta, phi);
		  //DEBUGGING
		  //theta = 75.8322;
		  //phi = 317.881;
		  collThetaPhi2Cart(temp_v, theta, phi, c, d, e);
		  col_c[i].v[0] = temp_v[0];
		  col_c[i].v[1] = temp_v[1];
		  col_c[i].v[2] = temp_v[2]; //now has unit velocity in correct direction in COM frame
			//NEW STUFF
			//ZML - commenting out binning stuff
			//			if (M_bin < col_c[i].fMass){ //if we use up all the mass in this bin
			//	v+=dv; //now in the next bin
			//	M_bin = M_avail*pow(10, A - S*v); //work out the mass in next bin
			//	while (M_bin < col_c[i].fMass){ // if still not enough, make bin wider (the bins decrease in size)
			//		dv *= 2.; //widen the new bin
			//		S = pow(10, A)/(2.3026*dv*(M_avail-col_c[rem_start-1].fMass)/M_avail); //re-calculate mass in new bin
			//		M_bin = M_avail*pow(10, A - S*v);
			//	}
			//	printf("\nM_bin new %E dv new %E particle.mass %E\n", M_bin, dv, col_c[i].fMass);
		  //	}
		  COPY_VEC(temp_v, temp_v2); //copy vector direction
		//			SCALE_VEC(temp_v, (v*0.5*dv)*v_esc); //scale vector to bin midpoint speed
		// ZML trying to fix weird eccentricity orbit problems ...
		vesc_i = rem_i%8+1;
		SCALE_VEC(temp_v, 1.1*vesc_i*v_esc); //I think v_esc should be V_esc according to LS12 ...
		temp_v_mag = MAG(temp_v);
		assert(temp_v_mag<(10.*v_esc));
            
		if( temp_v_mag>(10.*v_esc)){//cutoff at 10 v_esc
		  SCALE_VEC(temp_v2, 10.*v_esc);
		  COPY_VEC(temp_v2, col_c[i].v);
		  //col_c[i].vel = cdata.cmvel + 10.*v_esc*temp_v;
		}
		else{
		  COPY_VEC(temp_v, col_c[i].v);
		  printf("v_esc = %f vesc_i = %i temp_v = %f %f %f V_esc = %f\n", v_esc, vesc_i, temp_v[0], temp_v[1], temp_v[2], V_esc);
		  //col_c[i].vel = cdata.cmvel + (v+0.5*dv)*v_esc*temp_v; //add mid point of velocity bin to cmvel
				//		M_bin -= col_c[i].fMass; //take off mass 'used up' by this particle
		}//end of NEW STUFF
		rem_i++;
			
		//SCALE_VEC(col_c[i].v, sqrt(2.0*E_rem/col_c[i].fMass)); //scale the velocity to desired speed, here we are using equipartition of energy
		printf("remnant %d: Mass %G, |Vel|: %G, VX: %G, VY: %G, VZ: %G, v_esc %G\n", i, col_c[i].fMass, MAG(col_c[i].v), col_c[i].v[0], col_c[i].v[1], col_c[i].v[2], v_esc);
		ADD_VEC(col_c[i].v, vcom, col_c[i].v); //add centre of mass velocity (have been working in this frame)
		printf("remnant %d: Mass %G, |Vel|: %G, E: %G, Theta: %G, Phi: %G, VX: %G, VY: %G, VZ: %G\n",
		       i, col_c[i].fMass, MAG(col_c[i].v), 0.5*col_c[i].fMass*MAG_SQ(col_c[i].v), theta, phi,col_c[i].v[0], col_c[i].v[1], col_c[i].v[2]);
		// ZML humm - I think that the distance should be calculated with the relative velocity not the total velocity ...
		SUB_VEC(col_c[i].v,col_c[0].v,rel_vi);
		SUB_VEC(col_c[rem_start].v,col_c[0].v,rel_vslr);
		//	CSCALE_VEC(col_c[i].v, 2.0*targ->fRadius/MAG(col_c[rem_start].v), temp_v); //multiply by the time taken for slr to get to 2 target radii
		CSCALE_VEC(rel_vi, 20.0*targ->fRadius/MAG(rel_vslr), d_offset); //multiply by the time taken for slr to get to 20 target radii
		//ZML
		//shouldn't it be 		CSCALE_VEC(rel_vi, 20.0*col_c[0]->fRadius/MAG(rel_vslr), d_offset); //multiply by the time taken for slr to get to 20 target radii
		printf("The relative velocity is %e %e %e\n", rel_vi[0], rel_vi[1], rel_vi[2]);
		printf("The position correct is %e %e %e\n", d_offset[0], d_offset[1], d_offset[2]);
		if(M_tlr!=0.0){ //if we are a hit and run - projectile disrupt collision
				//add back on the momentum lost from projectile
		  ADD_VEC(col_c[i].v, v_temp_a, col_c[i].v); //ZML - is this correct?
		}
		//CSCALE_VEC(col_c[i].v, i*2.0*targ->fRadius/MAG(col_c[rem_start].v), temp_v); //multiply by the time taken for slr to get to 2 target radii //if you get overlap warnings, try using this version
		//SUB_VEC(col_c[1].v, vcom, temp_v2);
		//CSCALE_VEC(col_c[1].v, 1.0E2 - targ->fRadius/MAG(col_c[1].v), temp_v2);
		//SUB_VEC(temp_v, temp_v2, temp_v);
		printf("Before com correction rcom %e %e %e, col_c[i].r %e %e %e\n", rcom[0], rcom[1], rcom[2], col_c[i].r[0], col_c[i].r[1], col_c[i].r[2]);
		ADD_VEC(d_offset, rcom, col_c[i].r); //lr is at COM, add offset and assign position of remnant ZML - but this isn't always the case ... hr proj disrupt?
		printf("After com correction rcom %e %e %e, col_c[i].r %e %e %e\n", rcom[0], rcom[1], rcom[2], col_c[i].r[0], col_c[i].r[1], col_c[i].r[2]);
		// scale velocities to some decent number (maybe escape velocity?)
		// Zero spin and acceleration
		ZERO_VEC(col_c[i].a);
		ZERO_VEC(col_c[i].w);
		//printf("col_c[i] %p\n", col_c[i]);
		//printf("col_c %p, col_c[0] %p\n", (void *)col_c, &col_c[0]);
		//printf("col_c[%d] %p, ", i, &col_c[i]);
		//printf("col_c[%d].v %p\n", i, (void *)col_c[i].v);
		//printf("col_c[%d].v %G %G %G\n\n", i, col_c[i].v[0], col_c[i].v[1], col_c[i].v[2]);
		//CSCALE_VEC(col_c[i].v, col_c[i].fMass, temp_v);
		//ADD_VEC(temp_v, p_lr, p_lr); //add each remnant momentum 
		//ZML CHECK ecc and semi
		pos = sqrt(col_c[i].r[0]*col_c[i].r[0]+col_c[i].r[1]*col_c[i].r[1]+col_c[i].r[2]*col_c[i].r[2]);
		speed_2 = (col_c[i].v[0]*col_c[i].v[0]+col_c[i].v[1]*col_c[i].v[1]+col_c[i].v[2]*col_c[i].v[2]);
		semi = 1./(2./pos - speed_2);
		h_2 = (col_c[i].r[1]*col_c[i].v[2]-col_c[i].r[2]*col_c[i].v[1])*(col_c[i].r[1]*col_c[i].v[2]-col_c[i].r[2]*col_c[i].v[1]) +
		  (col_c[i].r[2]*col_c[i].v[0]-col_c[i].r[0]*col_c[i].v[2])*(col_c[i].r[2]*col_c[i].v[0]-col_c[i].r[0]*col_c[i].v[2]) +
		  (col_c[i].r[0]*col_c[i].v[1]-col_c[i].r[1]*col_c[i].v[0])*(col_c[i].r[0]*col_c[i].v[1]-col_c[i].r[1]*col_c[i].v[0]);
		ecc = sqrt(fabs(1.-h_2/semi));
		printf("particle %i at pos %f, %e %e %e has semi %f and ecc %f, offset %e %e %e\n", i, pos, col_c[i].r[0], col_c[i].r[1], col_c[i].r[2], semi, ecc, d_offset[0], d_offset[1], d_offset[2]);
		if (ecc > 1) 
		  printf("WARNING: ecc > 1\n");
		
		//assert(ecc < 1 && semi > 0);
	}
        pos = sqrt(col_c[0].r[0]*col_c[0].r[0]+col_c[0].r[1]*col_c[0].r[1]+col_c[0].r[2]*col_c[0].r[2]);
        speed_2 = (col_c[0].v[0]*col_c[0].v[0]+col_c[0].v[1]*col_c[0].v[1]+col_c[0].v[2]*col_c[0].v[2]);
        semi = 1./(2./pos - speed_2);
        h_2 = (col_c[0].r[1]*col_c[0].v[2]-col_c[0].r[2]*col_c[0].v[1])*(col_c[0].r[1]*col_c[0].v[2]-col_c[0].r[2]*col_c[0].v[1]) +
	  (col_c[0].r[2]*col_c[0].v[0]-col_c[0].r[0]*col_c[0].v[2])*(col_c[0].r[2]*col_c[0].v[0]-col_c[0].r[0]*col_c[0].v[2]) +
	  (col_c[0].r[0]*col_c[0].v[1]-col_c[0].r[1]*col_c[0].v[0])*(col_c[0].r[0]*col_c[0].v[1]-col_c[0].r[1]*col_c[0].v[0]);
        ecc = sqrt(fabs(1.-h_2/semi));
        printf("particle 0 at pos %f, %e %e %e has semi %f and ecc %f\n", pos, col_c[0].r[0], col_c[0].r[1], col_c[0].r[2], semi, ecc);
	if (ecc > 1) 
	  printf("WARNING: ecc > 1\n");
		//        assert(ecc < 1 && semi > 0);
	//at this point energy is p6roperly balanced!
	//this is for testing momentum and energy conservation
		
	{ //MOMENTUM AND ENERGY CHECK
	  //this is for testing momentum and energy conservation
	  SCALE_VEC(p_lr, 0.0); //zero this vector
	  E_rem = 0.0;
	  for(i=0;i<*nOut;i++){
	    VECTOR temp_v;
	    CSCALE_VEC(col_c[i].v, col_c[i].fMass, temp_v);
	    E_rem+=MAG_SQ(col_c[i].v)*col_c[i].fMass/2.0;
	    ADD_VEC(p_lr, temp_v, p_lr)
	      } //p_lr is now total momentum after tail + remnants have been made
	  printf("momentum in: %G %G %G, momentum out: %G %G %G, |p_i| %G, |p_o| %G\n", p_init[0], p_init[1], p_init[2],
		 p_lr[0], p_lr[1], p_lr[2], MAG(p_init), MAG(p_lr));
	  printf("Collision energies E_in: %G, E_out: %G\n", E_in, E_rem);
	  assert(E_in>E_rem);
	}
	SUB_VEC(p_init, p_lr, p_rem);
	NORM_VEC(p_rem, (*nOut-1));
        
        double E_remII;
	
        E_lr=0;
        E_remII=MAG_SQ(col_c[0].v)*col_c[0].fMass/2.0;
        if (rem_start==2) {
	  E_remII+=MAG_SQ(col_c[1].v)*col_c[1].fMass/2.0;
        }
        //recalculate energy
	for(i=rem_start;i<(*nOut);i++){ //(i=0;i<(*nOut -1);i++){
	  VECTOR temp1;
	  //double deltaE;
	  CSCALE_VEC(p_rem, 1.0/(col_c[i].fMass), temp1); //velocity adjustment to conserve momentum
	  SUB_VEC(col_c[i].v, temp1, temp1); //temp1 is the new velocity
	  //deltaE = 0.5*col_c[i].fMass*(MAG_SQ(temp1) - MAG_SQ(col_c[i].v)); //energy change between new and old
	  E_remII+=MAG_SQ(col_c[i].v)*col_c[i].fMass/2.0;
	}
        assert(E_in>E_remII);
        printf("There is more energy going in than going out E_in = %e and E_out = %e\n", E_in, E_remII);
	/*
	  j=0;
	  while(((E_in < E_rem) || (MAG(p_lr) > 1E-30)) && j<5 ){
	  j++;
	  SUB_VEC(p_init, p_lr, p_rem); //p_rem is the remaining momentum
	  printf("p_rem: %G %G %G\n", p_rem[0], p_rem[1], p_rem[2]);
	  //Ensures Conservation of Momentum
	  SUB_VEC(p_init, p_lr, p_rem);
	  NORM_VEC(p_rem, *nOut); //assign every particle with the same proportion
	  for (i=0; i<*nOut; i++){
	  VECTOR temp_v;
	  CSCALE_VEC(p_rem, 1.0/(col_c[i].fMass), temp_v);
	  ADD_VEC(col_c[i].v, temp_v, col_c[i].v); //add the remaining momentum to each remnant
	  //E_lr += 0.5*(MAG_SQ(col_c[i].v)*col_c[i].fMass); //add each remnant energy
	  
	  }
	  { //MOMENTUM AND ENERGY CHECK
	  //this is for testing momentum and energy conservation
	  SCALE_VEC(p_lr, 0.0); //zero this vector
	  E_rem = 0.0;
	  for(i=0;i<*nOut;i++){
	  VECTOR temp_v;
	  CSCALE_VEC(col_c[i].v, col_c[i].fMass, temp_v);
	  E_rem+=MAG_SQ(temp_v)/(2.0*col_c[i].fMass);
	  ADD_VEC(p_lr, temp_v, p_lr)
	  
	  }
	  printf("momentum in: %G %G %G, momentum out: %G %G %G, |p_i| %G, |p_o| %G\n", p_init[0], p_init[1], p_init[2],
	  p_lr[0], p_lr[1], p_lr[2], MAG(p_init), MAG(p_lr));
	  printf("Collision energies E_in: %G, E_out: %G\n", E_in, E_rem);
	  }
	  
	  if(E_rem > E_in){
	  printf("WARNING: Collision energies do not balance E_in: %G, E_out: %G\n", E_in, E_rem);
	  E_rem = (E_rem - E_in)/(*nOut); //find remaining energy per particle
	  //want to take off v = sqrt(2*E_rem/m) from each particle
	  for(i=0; i< *nOut; i++){ //take off the extra energy in equal momentum chunks (so mom is still conserved)
	  VECTOR temp_v;
	  CSCALE_VEC(col_c[i].v, sqrt(2.0*E_rem/col_c[i].fMass)/MAG(col_c[i].v), temp_v);
	  SUB_VEC(col_c[i].v, temp_v, col_c[i].v);
	  }
	  }
	  { //MOMENTUM AND ENERGY CHECK
	  //this is for testing momentum and energy conservation
	  SCALE_VEC(p_lr, 0.0); //zero this vector
	  E_rem = 0.0;
	  for(i=0;i<*nOut;i++){
	  VECTOR temp_v;
	  CSCALE_VEC(col_c[i].v, col_c[i].fMass, temp_v);
	  E_rem+=MAG_SQ(col_c[i].v)*col_c[i].fMass/2.0;
	  ADD_VEC(p_lr, temp_v, p_lr)
	  
	  }
	  printf("momentum in: %G %G %G, momentum out: %G %G %G, |p_i| %G, |p_o| %G\n", p_init[0], p_init[1], p_init[2],
	  p_lr[0], p_lr[1], p_lr[2], MAG(p_init), MAG(p_lr));
	  printf("Collision energies E_in: %G, E_out: %G\n", E_in, E_rem);
	  }
	  }
	*/
	// Check for momentum conservation
	if ((E_in < E_rem) || (MAG(p_rem) > 1E-30) ){
	  VECTOR temp1, temp2;
	  CSCALE_VEC(col_c[0].v, col_c[0].fMass, temp1);
	  CSCALE_VEC(col_c[1].v, col_c[1].fMass, temp2);
	  printf("WARNING: %d Energy or Momentum is not conserved within limits E_diff %G, p_mag %G\n",
		 j, (E_in - E_rem), MAG(p_rem));
	  printf("p_init %G p_lr %G p_slr %G\n", MAG(p_init), MAG(temp1), MAG(temp2));
	  printf("p_init %e %e %e\n",p_init[0],p_init[1],p_init[2]); //show initial momentum
	  printf("p_lr %e %e %e\n",temp1[0],temp1[1],temp1[2]); //show LR-momentum
	  printf("p_slr %G, %G, %G\n", temp2[0], temp2[1], temp2[2]);
	}
}
	for(i=0;i<*nOut;i++){
		VECTOR v_diff;
		SUB_VEC(col_c[i].v, vcom, v_diff);
		printf("col_c[%d].id.iOrgIdx %d bFrag %d\n",i, col_c[i].id.iOrgIdx, col_c[i].bFrag);
		printf("col_c[%d].v - vcom %E %E %E mag %E mag/v_esc %E\n", i, v_diff[0], v_diff[1], v_diff[2], MAG(v_diff), MAG(v_diff)/v_esc);
	}
	if (targ->iRung >= proj->iRung)
		for(i=0;i<*nOut;i++)
			col_c[i].iRung = targ->iRung;	
	else
		for(i=0;i<*nOut;i++)
			col_c[i].iRung = proj->iRung;
	
#ifdef ORIGIN_HISTOGRAM
	for (j=0;j<*nOut;j++) {
			for(i = 0; i < NUM_ORIGIN_BINS; ++i) {
				col_c[j].origin_bins[i] = targ->origin_bins[i];
			}
			MergeHistograms(targ->fMass, col_c[j].origin_bins, proj->fMass, proj->origin_bins);
		}
#endif
	//JD: DEBUG
	*col_c1 = col_c;
	assert(*col_c1 != NULL);
}

/*DEBUG following function fails pedantic compilation checking due to
  too many unused variables...*/

void collCollModCollide(COLLIDER *col_a, COLLIDER *col_b,
						const COLLISION_PARAMS *cp, COLLIDER **col_c,
						double *dMassInDust,int *nOut, PKD pkd)
{
	
  // Determines collision outcome based on collision model from Leinhardt & Stewart 2012 */
  
  // Questions:
  // 1. What do we do about spin? 
  // 2. Origin histograms. 
  // 3. Check units ... right now in sys units ... 
  // 4. Deal with calculating c_star from mu_bar 
  
  COLLIDER *proj, *targ;
  VECTOR vel,loc,cross;
  BOOLEAN graze, hit_and_run = 0;
  int N_lr, N_slr, N_tail, bsupercat=0;
  int N_tlr=0;
  double M_tlr=0.0, D_tlr=0.0;
  double G = 6.67e-8;
  double loc_mag,vel_mag,l,b,b_crit,alpha,Q_R,Q_supercat, V_supercat;
  double M,R,rho,V_esc, M_tail, M_lim, C, D_slr, D_lim, rhoint;
  // C_star is a material parameter for small bodies 5 and large fluid bodies 2
  double c_star = 5, mu_bar = 0.37; 
  // Change from system mass to cgs then back again for the first pass ...
  // Careful rho_1 is 1 g/cm^3 so if not cgs should not be 1
  double R_C1,Q_PD,V_PD,rho_1=1.,Q_REr,V_Er,M_lr,M_slr; 
  double beta,eta, M_tot, mu, mu_int,gamma,Q_RDstar,V_star,Q_RDstar_prime,V_star_prime;
  double M_tot_dagger, R_C1_dagger, Q_PD_dagger, V_PD_dagger, mu_dagger, gamma_dagger;
  double V_RDstar_dagger, Q_RDstar_dagger, V_supercat_dagger, Q_REr_dagger,V_Er_dagger, Q_R_dagger, targ_radcgs, proj_radcgs, targ_masscgs,proj_masscgs;
  double mu_int_dagger, V_star_dagger, Q_RDstar_prime_dagger, V_star_prime_dagger, Q_supercat_dagger;
  double phi_dagger, M_dagger, A_int_dagger, L_int_dagger, alpha_dagger, M_int_dagger;
  assert(*col_c == NULL); /* just to be sure */


  //AB
  printf(" c1 Mass: %e c2 Mass %e \n", col_a->fMass, col_b->fMass);
  //AB

  //printf("IN %s\n", __FUNCTION__);
  printf("nThreads: %d, idSelf: %d\n", pkd->nThreads, pkd->idSelf);
  // 1. First step determine interacting mass - Eqn. 11 L&S 2012
  
  // Determine projectile and target because interaction 
  // Mass is defined in terms of projectile mass
  //printf("Collider Velocities: Proj %e %e %e Targ %e %e %e\n",col_a->v[0],col_a->v[1],col_a->v[2],col_b->v[0],col_b->v[1],col_b->v[2]);
  //printf("Collider Positions:  Proj %e %e %e Targ %e %e %e\n",col_a->r[0],col_a->r[1],col_a->r[2],col_b->r[0],col_b->r[1],col_b->r[2]);
  //printf("Collider Spins:  Proj %e %e %e Targ %e %e %e\n",col_a->w[0],col_a->w[1],col_a->w[2],col_b->w[0],col_b->w[1],col_b->w[2]);
  //printf("Collider radii: Proj %G, Targ %G\n", col_a->fRadius, col_b->fRadius);
  //printf("Collider masses: Proj %G, Targ %G\n", col_a->fMass, col_b->fMass);
  if (col_a->fMass < col_b->fMass) {
    proj = col_a;
    targ = col_b;
  } else {
    proj = col_b;
    targ = col_a;
  }
	
  printf("Collider Velocities: Proj %e %e %e Targ %e %e %e\n",proj->v[0],proj->v[1],proj->v[2],targ->v[0],targ->v[1],targ->v[2]);
  printf("Collider Positions:  Proj %e %e %e Targ %e %e %e\n",proj->r[0],proj->r[1],proj->r[2],targ->r[0],targ->r[1],targ->r[2]);
  printf("Collider Spins:  Proj %e %e %e Targ %e %e %e\n",proj->w[0],proj->w[1],proj->w[2],targ->w[0],targ->w[1],targ->w[2]);
  printf("Collider radii: Proj %G, Targ %G\n", proj->fRadius, targ->fRadius);
  printf("Collider masses: Proj %G, Targ %G\n", proj->fMass, targ->fMass);
  // Change units 
  SUB_VEC(proj->v,targ->v,vel);
  SCALE_VEC(vel,VEL_CGS);
  vel_mag = MAG(vel);
  SUB_VEC(proj->r,targ->r,loc);
  SCALE_VEC(loc,LEN_CGS);
  loc_mag = MAG(loc);
  CROSS(vel,loc,cross);
  
  // Calculate the impact parameter!
  b = MAG(cross)/(loc_mag*vel_mag);
  printf("loc_mag: %G, vel_mag: %G, b: %G\n", loc_mag, vel_mag, b);
  printf("Impact Angle: %f\n",b);
  printf("Impact Velocity: %e\n",vel_mag);
  
  targ_radcgs = targ->fRadius*LEN_CGS;
  proj_radcgs = proj->fRadius*LEN_CGS;
  l = (targ_radcgs + proj_radcgs)*(1.0-b);
  //CARFUL !!! l is not defined for proj in shadow of targ ZML
  if (targ_radcgs > (b*(proj_radcgs + targ_radcgs) + proj_radcgs)){
    alpha = 1.0;
  }else{
    alpha = (3.*proj_radcgs*SQ(l) - CUBE(l))/(4.*CUBE(proj_radcgs));
  }

  // 2. Check if in perfect merging regime
  
  targ_masscgs = targ->fMass*MASS_CGS;
  proj_masscgs = proj->fMass*MASS_CGS;
  
  M = targ_masscgs + alpha*proj_masscgs;
  // Assuming both target and projectile have same rho!
  rho = targ_masscgs/(4./3.*PI*CUBE(targ_radcgs)); 
  rhoint = rho;
  R = pow((3.*M)/(4*PI*rho), 1./3.);
  V_esc = sqrt(2*G*M/R);
  printf("Escape Velocity = %e\n", V_esc);
  printf("Projectile Mass %e Target Mass %e\n",proj_masscgs,targ_masscgs);
  
  // If velocity less that escape velocity then
  // We have a perfect merging event occuring
  
  if (vel_mag < V_esc) {
    printf("OUTCOME: PERFECT MERGING\n");
    *nOut = 1;
    
    // Make sure the units have been changed back to sys units - units of col_a, col_b, and col_c are in sys units
    // Call merging function
    
    collCreatePart(col_a,col_b,col_c);
    
    // No dust produced in perfect merger
    *dMassInDust = 0.0;
  } 

  // 3. Calculate the critical impact parameter
  
  // If velocity greater than the escape velocity
  // We have some non-perfect-merging event
  
  else { 
    b_crit = targ_radcgs/(targ_radcgs + proj_radcgs);
    M_tot = proj_masscgs + targ_masscgs;
    
    // If non-grazing collision
    if (b < b_crit)
      graze = 0;
    
    // If grazing collection
    else
      graze = 1;
    
    // 4. Calculate the catastrophic disruption criterion
    
    //    A. Calculate R_C1 from total mass and density of 1 g/cm^3
    R_C1 = pow((3.*M_tot/(4.*PI*rho_1)),(1./3.));
    
    //    B. Calculate pricipal disruption value for an equivalent equal mass
    //       At R_C1, Q^*_RD,gamma=1 (eqn. 28) and its corresponding critical
    //		 Impact velocity (eqn.30)
    Q_PD = c_star*4./5.*PI*rho_1*G*SQ(R_C1);
    V_PD = sqrt(32.*PI*c_star*rho_1*G/5)*R_C1;
    printf("Q_PD: %G, V_PD: %G\n", Q_PD, V_PD);
    //	  C. Calculate the reduced mass and the reduced mass using the interacting mass eqn.12
    mu = targ_masscgs*proj_masscgs/M_tot;
    mu_int = alpha*targ_masscgs*proj_masscgs/(targ_masscgs+alpha*proj_masscgs);
    
    
    //	  D. Calculate the disruption criterion and critical impact velocity eqn. 23 & 22
    gamma = proj_masscgs/targ_masscgs;
    Q_RDstar = Q_PD*(pow((SQ(gamma+1)/(4*gamma)),(2/(3*mu_bar)-1)));
    V_star = V_PD*pow((SQ(gamma+1)/(4*gamma)),(1/(3*mu_bar)));
    
    //	  E. disruption energy and critical velocity with impact angle
    Q_RDstar_prime = pow(mu/mu_int,(2-3*mu_bar/2))*Q_RDstar;
    V_star_prime = sqrt(2.0*Q_RDstar_prime*M_tot/mu);
    
    // 5. Calculate Q_R for erosion
    
    Q_REr = Q_RDstar_prime*(-2.0*(targ_masscgs/M_tot)+2.0);
    V_Er = sqrt(2.0*Q_REr*M_tot/mu);
    
    // 6. Hit-and-run regime
    //DEBUG
    Q_supercat = 1.8*Q_RDstar_prime;
    V_supercat = sqrt(2.0*Q_supercat*(M_tot)/mu);
    printf("V_supercat: %G\n", V_supercat);
    // If collision is grazing but also the velocity
    // Is less than that required to erode the target
    // Then initialise the hit and run regime
    printf("V_Er: %G\n", V_Er);
    if (graze && (vel_mag < V_Er)) {
      // This has to be symmetric i.e. if mr =1 should give same result as normal collison
      // Mass of largest remnant = mass of target - proj: may or may not suffer erosion ...
      M_lr = targ_masscgs; 
      // What happens to the projectile? Need to calculate reverse impact.
      if(gamma<=0.5){ //do backwards - carve tunnel - calculation
	phi_dagger = 2.0*acos((l-proj_radcgs)/proj_radcgs); //just above eq. 46
	A_int_dagger = proj_radcgs*proj_radcgs*(M_PI-(phi_dagger - sin(phi_dagger))/2.0);
	L_int_dagger = 2.0*sqrt(targ_radcgs*targ_radcgs - (targ_radcgs - l/2.0)*(targ_radcgs - l/2.0));
	M_int_dagger = A_int_dagger*L_int_dagger*rho;
	alpha_dagger = M_int_dagger/targ_masscgs;
	if (alpha_dagger < 0.0) {
	  printf("WARNING: alpha_dagger < 0 %G\n", alpha_dagger);
	  alpha_dagger = 1.0;
	}
	// Assert: grazing collision and projectile is never fully in target shadow ZML
	assert(targ_radcgs < (b*(proj_radcgs + targ_radcgs) + proj_radcgs));
	R_C1 = pow((3.*(M_int_dagger+proj_masscgs)/(4.*PI*rho_1)),(1./3.));
	Q_PD = c_star*4./5.*PI*rho_1*G*SQ(R_C1);
	V_PD = sqrt(32.*PI*c_star*rho_1*G/5)*R_C1;
	M_tot_dagger = M_int_dagger + proj_masscgs;
				//alpha = (3.*targ_radcgs*SQ(l) - CUBE(l))/(4.*CUBE(targ_radcgs));
				//if (alpha < 0.0) {
				//	alpha = 1.0;
				//}
	
				//M = proj_masscgs + alpha*targ_masscgs;
				//M_tot_dagger = M;
				/*
				M_tot_dagger = M_tot
				R_C1_dagger = pow(3.*M_tot_dagger/(4.*PI*rho_1),(1.0/3.0));
				Q_PD_dagger = c_star*4./5.*PI*rho_1*G*SQ(R_C1_dagger);
				V_PD_dagger = sqrt(32.*PI*c_star*rho_1*G/5)*R_C1_dagger;
				mu_dagger = M*proj_masscgs/M_tot_dagger;
				gamma_dagger = M/proj_masscgs;
				V_RDstar_dagger = pow((0.25*(gamma_dagger+1)*(gamma_dagger+1)/gamma_dagger),1/(3*mu_bar))*V_PD_dagger;
				Q_RDstar_dagger = Q_PD_dagger*pow((0.25*(gamma_dagger + 1)*(gamma_dagger + 1)/gamma_dagger),(2/(3*mu_bar)-1));
				*/
	
	//JD: DEBUGGING REVERSE IMPACT
	mu_int_dagger = alpha_dagger*targ_masscgs*proj_masscgs/(alpha_dagger*targ_masscgs + proj_masscgs);
	printf("alpha_dagger: %G, alpha: %G\n", alpha_dagger, alpha);
	printf("mu_int_dagger: %G, mu_int: %G\n", mu_int_dagger, mu_int);
	//	  D. Calculate the disruption criterion and critical impact velocity eqn. 50
	gamma_dagger = M_int_dagger/proj_masscgs;
	printf("gamma_dagger: %G, gamma: %G\n", gamma_dagger, gamma);
	Q_RDstar_dagger = Q_PD*(pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(2/(3*mu_bar)-1)));
	V_star_dagger = V_PD*pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(1/(3*mu_bar)));
	printf("Q_RDstar_dagger: %G, Q_RDstar: %G\n", Q_RDstar_dagger, Q_RDstar);
				printf("V_star_dagger: %G, V_star: %G\n", V_star_dagger, V_star);
			
				//	  E. disruption energy and critical velocity with impact angle
				Q_RDstar_prime_dagger = Q_RDstar_dagger; //only in backwards case see eq. 17
				printf("Q_RDstar_prime_dagger: %G, Q_RDstar_prime: %G\n", Q_RDstar_prime_dagger, Q_RDstar_prime);
				V_star_prime_dagger = sqrt(2.0*Q_RDstar_prime_dagger*M_tot_dagger/mu_int_dagger);
				printf("V_star_prime_dagger: %G, V_star_prime: %G\n", V_star_prime_dagger, V_star_prime);
				// 5. Calculate Q_R for erosion

				Q_REr_dagger = Q_RDstar_prime_dagger*(-2.0*(proj_masscgs/M_tot_dagger)+2.0);
				printf("Q_REr_dagger: %G, Q_REr: %G\n", Q_REr_dagger, Q_REr);
				V_Er_dagger = sqrt(2.0*Q_REr_dagger*M_tot_dagger/mu_int_dagger);
				//JD: DEBUG ENDS HERE
			
				//Just test for disruption of proj
				//V_supercat_dagger = 1.8*V_star_dagger;
				Q_supercat_dagger = 1.8*Q_RDstar_prime_dagger;
				V_supercat_dagger = sqrt(2.0*Q_supercat_dagger*(M_tot_dagger)/mu_int_dagger);
				printf("V_supercat_dagger: %G, vel_mag: %G\n", V_supercat_dagger, vel_mag);
				printf("V_ER_dagger: %G\n", V_Er_dagger);
				mu_dagger = mu_int_dagger;
			}
			else{ //do forwards calculation
				mu_dagger = mu;
				M_tot_dagger = M_tot;
				alpha_dagger = (3.*targ_radcgs*SQ(l) - CUBE(l))/(4.*CUBE(targ_radcgs));
				//M_dagger = alpha_dagger*targ_masscgs + proj_masscgs;
				M_int_dagger = alpha_dagger*targ_masscgs;
				if (alpha_dagger < 0.0) {
					printf("WARNING: alpha_dagger < 0 %G\n", alpha_dagger);
					alpha_dagger = 1.0;
				}
				// Assert: grazing collision and projectile is never fully in target shadow ZML
				assert(targ_radcgs < (b*(proj_radcgs + targ_radcgs) + proj_radcgs));
				mu_int_dagger = alpha_dagger*targ_masscgs*proj_masscgs/(alpha_dagger*targ_masscgs+proj_masscgs);

		
				//	  D. Calculate the disruption criterion and critical impact velocity eqn. 23 & 22
				gamma_dagger = targ_masscgs/proj_masscgs;
				Q_RDstar_dagger = Q_PD*(pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(2/(3*mu_bar)-1)));
				V_star_dagger = V_PD*pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(1/(3*mu_bar)));
		
				//	  E. disruption energy and critical velocity with impact angle
				Q_RDstar_prime_dagger = pow(mu/mu_int_dagger,(2-3*mu_bar/2))*Q_RDstar_dagger;
				V_star_prime_dagger = sqrt(2.0*Q_RDstar_prime_dagger*M_tot/mu);

				// 5. Calculate Q_R for erosion

				Q_REr_dagger = Q_RDstar_prime_dagger*(-2.0*(proj_masscgs/M_tot)+2.0);
				V_Er_dagger = sqrt(2.0*Q_REr_dagger*M_tot/mu);

				// 6. Hit-and-run regime
				//DEBUG
				Q_supercat_dagger = 1.8*Q_RDstar_prime_dagger;
				V_supercat_dagger = sqrt(2.0*Q_supercat_dagger*(M_tot)/mu);
				printf("V_supercat_dagger: %G, vel_mag: %G\n", V_supercat_dagger, vel_mag);
				printf("V_ER_dagger: %G\n", V_Er_dagger);
			}
			/*
			R_C1 = pow((3.*(M_int_dagger+proj_masscgs)/(4.*PI*rho_1)),(1./3.));
			Q_PD = c_star*4./5.*PI*rho_1*G*SQ(R_C1);
			V_PD = sqrt(32.*PI*c_star*rho_1*G/5)*R_C1;
			M_tot_dagger = M_int_dagger + proj_masscgs;
			//alpha = (3.*targ_radcgs*SQ(l) - CUBE(l))/(4.*CUBE(targ_radcgs));
			//if (alpha < 0.0) {
			//	alpha = 1.0;
			//}
	
			//M = proj_masscgs + alpha*targ_masscgs;
			//M_tot_dagger = M;
			//JD: DEBUGGING REVERSE IMPACT
			mu_int_dagger = alpha_dagger*targ_masscgs*proj_masscgs/(alpha_dagger*targ_masscgs + proj_masscgs);
			printf("alpha_dagger: %G, alpha: %G\n", alpha_dagger, alpha);
			printf("mu_int_dagger: %G, mu_int: %G\n", mu_int_dagger, mu_int);
			//	  D. Calculate the disruption criterion and critical impact velocity eqn. 50
			gamma_dagger = M_int_dagger/proj_masscgs;
			printf("gamma_dagger: %G, gamma: %G\n", gamma_dagger, gamma);
			Q_RDstar_dagger = Q_PD*(pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(2/(3*mu_bar)-1)));
			V_star_dagger = V_PD*pow((SQ(gamma_dagger+1)/(4*gamma_dagger)),(1/(3*mu_bar)));
			printf("Q_RDstar_dagger: %G, Q_RDstar: %G\n", Q_RDstar_dagger, Q_RDstar);
			printf("V_star_dagger: %G, V_star: %G\n", V_star_dagger, V_star);
			
			//	  E. disruption energy and critical velocity with impact angle
			Q_RDstar_prime_dagger = pow(mu/mu_int_dagger, 2.0-(3.0*mu_bar/2.0))*Q_RDstar_dagger;
			printf("Q_RDstar_prime_dagger: %G, Q_RDstar_prime: %G\n", Q_RDstar_prime_dagger, Q_RDstar_prime);
			V_star_prime_dagger = sqrt(2.0*Q_RDstar_prime_dagger*M_tot_dagger/mu_int_dagger);
			printf("V_star_prime_dagger: %G, V_star_prime: %G\n", V_star_prime_dagger, V_star_prime);
			// 5. Calculate Q_R for erosion

			Q_REr_dagger = Q_RDstar_prime_dagger*(-2.0*(proj_masscgs/M_tot_dagger)+2.0);
			printf("Q_REr_dagger: %G, Q_REr: %G\n", Q_REr_dagger, Q_REr);
			V_Er_dagger = sqrt(2.0*Q_REr_dagger*M_tot_dagger/mu_int_dagger);
			//JD: DEBUG ENDS HERE
			
			//Just test for disruption of proj
			//V_supercat_dagger = 1.8*V_star_dagger;
			Q_supercat_dagger = 1.8*Q_RDstar_prime_dagger;
			V_supercat_dagger = sqrt(2.0*Q_supercat_dagger*(M_tot_dagger)/mu_int_dagger);
			printf("V_supercat_dagger: %G, vel_mag: %G\n", V_supercat_dagger, vel_mag);
			
			printf("V_ER_dagger: %G\n", V_Er_dagger);
			// Determine if projectile is disrupted or not
			*/
			if (vel_mag > V_Er_dagger) {
				Q_R_dagger = 0.5*mu_dagger*vel_mag*vel_mag/M_tot_dagger;
				
				// Determine is catastrophic or super catastrophic
				
				if (vel_mag > V_supercat_dagger) {
					printf("Outcome: HIT & RUN - PROJECTILE SUPERCATASTROPHICALLY DISRUPTED\n");
					eta = -1.5;
					M_slr = M_tot_dagger*0.1/pow(1.8,eta)*pow((Q_R_dagger/Q_RDstar_dagger),eta);
					printf("M_slr %G, M_tot_dagger: %G, eta: %G, Q_R_dagger: %G, Q_RDstar_dagger: %G, \n",
						M_slr, M_tot_dagger, eta, Q_R_dagger, Q_RDstar_dagger);
				} else {
					printf("Outcome: HIT & RUN - PROJECTILE DISRUPTED\n");
					M_slr = M_tot_dagger*(-0.5*(Q_R_dagger/Q_RDstar_dagger - 1) + 0.5);
					printf("M_slr %G, M_tot_dagger: %G, eta: %G, Q_R_dagger: %G, Q_RDstar_dagger: %G, \n",
						M_slr, M_tot_dagger, eta, Q_R_dagger, Q_RDstar_dagger);
				}
			} else {
				printf("Outcome: HIT & RUN - PROJECTILE INTACT\n");
				M_slr = proj_masscgs;
				hit_and_run = 1;
			}

			
		} else { 

	        // 7. Determine super-catastrophic disruption regime
			
			//printf("Determine if in super-cat regime\n");
			Q_supercat = 1.8*Q_RDstar_prime;
			V_supercat = sqrt(2.0*Q_supercat*(M_tot)/mu);

			if (vel_mag > V_Er) {
				// But what if this is not the case either. If vel_mag > V_esc but < V_Er - partial accretion...
				//printf("In erosion regime need to determine if supercatastrophic or just disruptive, vel_mag = %e, V_Er = %e\n", vel_mag, V_Er);

	     	    // 8. In Erosion regime
				
				// In erosion regime regardless of impact angle!
				
				Q_R = 0.5*mu*vel_mag*vel_mag/M_tot;
				if (vel_mag > V_supercat) {

	        		// 9. In Super-catastrophic regime

					// In supercatastrophic regime regardless of impact angle!
					
					printf("Outcome: SUPERCATASTROPHIC DISRUPTION\n");
					eta = -1.5;
					M_lr = M_tot*0.1/pow(1.8,eta)*pow((Q_R/Q_RDstar_prime),eta);
					bsupercat=1;
					
				} else {

	        		// 10. In non-grazing disruption regime

        			// Not seperately identifying partial accretion ...

	        		// 11. In grazing disruption regime

				    // Largest rem is determined by universal law
				    
					printf("Outcome: EROSIVE DISRUPTION\n");	
					M_lr = M_tot*(-0.5*(Q_R/Q_RDstar_prime - 1) + 0.5);
				}
				
				N_lr = 1;
				N_slr = 2;
				beta = 2.85;
				M_slr = M_tot*(3-beta)*(1-N_lr*M_lr/M_tot)/(N_slr*beta);
			
			// If velocity is not large enough to erode the structure then
			// It will cause some of the fragmented particle to accrete onto the second
			} else {
				Q_R = 0.5*mu*vel_mag*vel_mag/M_tot;
				printf("Outcome: PARTIAL ACCRETION\n");
				M_lr = M_tot*(-0.5*(Q_R/Q_RDstar_prime - 1) + 0.5);
				N_lr = 1;
				N_slr = 2;
				beta = 2.85;
				M_slr = M_tot*(3-beta)*(1-N_lr*M_lr/M_tot)/(N_slr*beta);
			}
		}
		printf("M_lr: %E\n", M_lr);
		printf("M_slr: %E\n", M_slr);
		M_tail = 0.0;

		// Supercatastrophic Fail Hackery
		if (M_lr < M_slr) {
			M_slr = M_lr;
		}

		M_lim = cp->DB.dCollMinMass*MASS_CGS; //CGS
		printf("nOut: %d, M_lim: %E\n", *nOut, M_lim);
		

		// If even the mass of the largest remnant lies below the mass resolution limit
		// Then all collisional mass goes into dust

		if (M_lr < M_lim) {

			// No physical particles leave CollMod
			*nOut = 0;
			*dMassInDust = M_tot/MASS_CGS;
			*col_c = NULL;
			printf("nOut: %d, M_lim: %E\n", *nOut, M_lim);
		} 

		// If the second largest mass lies below resolution then put all mass
		// Except the LR into dust

		else if (M_slr < M_lim) {

			// Only LR leaves CollMod
			*nOut = 1;
			// If ex-tail remnant then merge! /*DEBUG this shouldn't happen anymore!*/
			if (col_a->iColor != 3 && col_b->iColor != 3) {
				printf("MERGING TAIL REMNANTS\n");
			    *dMassInDust = 0.0;
			    *nOut = 1;
			    collCreatePart(col_a,col_b,col_c);
			    goto endthis;
			}
			else {
				puts("DEBUG 1");
				*dMassInDust = (M_tot - M_lr)/MASS_CGS;
				//JD: DEBUG
				//*col_c = (COLLIDER *) malloc(*nOut*sizeof(COLLIDER));
				//assert(*col_c != NULL);
				// C is not needed if the only remnant resolved is the largest one
				C = 0.0; 
				//collCreateColliders(nOut,b,C,M_lr,M_slr,rhoint,col_a,col_b,&M_tail,*col_c,pkd, alpha, V_esc, bsupercat);
				collCreateColliders(nOut,b,C,M_lr,M_slr,M_tlr,rhoint,col_a,col_b,&M_tail,col_c,pkd, alpha, V_esc, bsupercat);
				puts("HERE1");
			}

		}
		// Otherwise deal with 2+ particle system
		else {
		  // If hit and run event then just set two particles out
		  if (hit_and_run == 1) {
		    *nOut = 2;
		    *dMassInDust = 0;
			/*DEBUG malloc handled by call to pkdBounce() in this case...
		    *col_c = (COLLIDER *) malloc(*nOut*sizeof(COLLIDER));
		    assert(*col_c != NULL);
		    (*col_c)[1]->bFrag = 0;
			*/
		    C = 0.0;
		    
		    // Calculate tail distribution
		  } 
		  else {
		    double D_truelim;
			
			// If ex-tail remnant then merge!
			if (col_a->iColor != 3 && col_b->iColor != 3) {
				printf("MERGING EX_TAIL REMNANT\n");
			    *dMassInDust = 0.0;
			    *nOut = 1;
			    collCreatePart(col_a,col_b,col_c);
			    goto endthis;
			}
			else {
				
				// Calculate number of particles larger than resolution limit
                // commented out on suggestion of zoe XXX
				//if (rho < 1) rho = 2.0;
				if(graze && M_lr==targ_masscgs && M_slr!=proj_masscgs){
					//if grazing collision, M_slr is calculated as the largest remnant
					//need to calculate third largest remnant + tail
					
					N_slr = 1;
					N_tlr = 2;
					beta = 2.85;
					M_tlr = M_tot_dagger*(3-beta)*(1-N_slr*M_slr/proj_masscgs)/(N_tlr*beta);
					if (M_tlr < M_lim) {
						// Only LR leaves CollMod
						*nOut = 2;
						puts("DEBUG M_tlr < M_lim");
						*dMassInDust = (M_tot - M_lr - M_slr)/MASS_CGS;
						// C is not needed if the only M_lr and M_slr resolved
						C = 0.0; 
						collCreateColliders(nOut,b,C,M_lr,M_slr,M_tlr,rhoint,col_a,col_b,&M_tail,col_c,pkd, alpha, V_esc, bsupercat);
						//Check M_tail is zeroed in function
						goto endthis;
					}
					C = N_tlr*beta*(pow(((3-beta)*(proj_masscgs-N_slr*M_slr)/((4/3.)*PI*rho*N_tlr*beta)),beta/3));
					D_tlr = pow(3.*M_tlr/(4.*PI*rho),(1./3.));
					//D_slr = D_tlr; // In the case of hit and run proj disruption the tail is from the projectile. ZML
					D_truelim = pow((M_lim/(4./3.*PI*rho)),(1/3.));
					D_lim = pow((pow(D_tlr,(3-beta))-(proj_masscgs-(M_slr))*(3-beta)/(4./3.*PI*C*rho)),1/(3-beta));	//PJC: M_tlr should not be included in this
				}
				else{
					C = N_slr*beta*(pow(((3-beta)*(M_tot-N_lr*M_lr)/((4/3.)*PI*rho*N_slr*beta)),beta/3));	//PJC 24/07/2014: Was missing brackets around the denominator, was only dividing by the 4/3.
					D_slr = cbrt(3.*M_slr/(4.*PI*rho));
					printf("D_slr(from M_slr): %e, D_slr: %e\n",D_slr,pow((C/(N_slr*beta)),(1./beta)));
					//D_slr=pow((N_slr*beta/C),(-1./beta));
					D_truelim = cbrt((M_lim/(4./3.*PI*rho)));
					D_lim = pow((pow(D_slr,(3-beta))-(M_tot-(N_lr*M_lr))*(3-beta)/(4./3.*PI*C*rho)),1/(3-beta));	//PJC 25/07/2014: M_slr should not be included in this calculation. D_lim should be zero.
					printf("D_lim: %e, %e\n",D_lim,pow((pow(D_slr,(3-beta))-(M_tot-(M_lr))*(3-beta)/(4./3.*PI*C*rho)),1/(3-beta)));
					printf("%.20e %.20e %e\n",pow(D_slr,(3-beta)),(M_tot-(M_lr))*(3-beta)/(4./3.*PI*C*rho),1/(3-beta));
				}
				if (D_lim>D_truelim) {
					;
				} 
				else {
					D_lim = D_truelim;
				}
				if(M_tlr==0.0)
					N_tail = (int)(C/beta*(pow(D_lim,-beta) - pow(D_slr,-beta)));
				else
					N_tail = (int)(C/beta*(pow(D_lim,-beta) - pow(D_tlr,-beta)));
				printf("N_tail: %d\n",N_tail);
				assert(N_tail < (2*MAX_PART));
				assert(N_tail >= 0);
				if (N_tail == 0 && M_tlr==0.0) { //JD: what about M_tlr >0? //PJC: This bit is redundant now, but keep it for now to be safe
					*nOut = 2; //PJC 1
					*dMassInDust = (M_tot - M_lr - M_slr)/MASS_CGS;	//PJC (M_tot - M_lr)/MASS_CGS
					//JD: DEBUG
					puts("DEBUG 2");
					//*col_c = (COLLIDER *) malloc(*nOut*sizeof(COLLIDER));
					//assert(*col_c != NULL);
					C = 0.0; 
					//collCreateColliders(nOut,b,C,M_lr,M_slr,rhoint,col_a,col_b,&M_tail,*col_c,pkd, alpha, V_esc, bsupercat);
					collCreateColliders(nOut,b,C,M_lr,M_slr,M_tlr,rhoint,col_a,col_b,&M_tail,col_c,pkd, alpha, V_esc, bsupercat);
					goto endthis;
				}
				else if(N_tail==0 && M_tlr !=0.0){ //PJC: This bit is redundant now, but keep it for now to be safe
					//we now have 2 remnants to output
					*nOut = 3; //PJC 2
					*dMassInDust = (M_tot - M_lr - M_slr - M_tlr)/MASS_CGS; //PJC (M_tot - M_lr - M_slr)/MASS_CGS
					C = 0.0; //not resolving tail remnants
					collCreateColliders(nOut,b,C,M_lr,M_slr,M_tlr,rhoint,col_a,col_b,&M_tail,col_c,pkd, alpha, V_esc, bsupercat);
					goto endthis;
				}
				else{
				
				if(M_tlr==0.0)
					*nOut = N_tail + 2; //PJC 1
				else
					*nOut = N_tail + 3; //PJC 2

				
				// Need to know mass in resolved tail in order to set MassInDust
				//JD: DEBUG
				puts("DEBUG 3");
				//*col_c = (COLLIDER *) malloc(*nOut*sizeof(COLLIDER));
				//assert(*col_c != NULL);
				//collCreateColliders(nOut,b,C,M_lr,M_slr,rhoint,col_a,col_b,&M_tail,*col_c,pkd, alpha, V_esc, bsupercat);
				collCreateColliders(nOut,b,C,M_lr,M_slr,M_tlr,rhoint,col_a,col_b,&M_tail,col_c,pkd, alpha, V_esc, bsupercat);
				//M_temp = 4./3.*PI*C*rho/(3-beta)*(pow(D_slr,(3-beta))-pow(D_lim,(3-beta)));
				*dMassInDust = 1.0/MASS_CGS*(M_tot - (M_lr + M_slr + M_tlr + MASS_CGS*M_tail)); //SYS	
				}
				assert(*dMassInDust>0);
			}
		  }
		}
	}
	endthis: ;
	printf("After endthis nOut: %d, M_lim: %E\n", (*nOut), M_lim);
	//	{
	//int z;
	//double density_postcoll;

	//for (z=0;z<(*nOut);z++) {
	//  density_postcoll = ((*col_c)[z].fMass*2.0e33)/((4.0/3.0)*M_PI*CUBE((*col_c)[z].fRadius*1.5e13));
	//  printf("COLLMOD: The density of particle %i is %e\n", z, density_postcoll);
	  //if (density_postcoll < 0.009)
	  //  printf("COLLMOD DENSITY ERROR!\n");
	  //}
	//}
}

#endif /* COLLMOD_ZML*/
