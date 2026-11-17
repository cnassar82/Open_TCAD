/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *     Copyright c 1991 The board of regents of the University of 	 *
 *                      Florida.  All rights reserved.			 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of the University of 	 *
 *     Florida.								 *
 *									 *
 *									 *
 *************************************************************************/
/*   solve.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:50:15 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <sys/param.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "material.h"
#include "matrix.h"
#include "defect.h"

/* HZ is 100 per Mike Eldredge on the Convex,
   which should be included in the file sys/param, above*/
#ifdef CONVEX
#define HZ 100
#endif
#ifndef HZ
#define HZ 60
#endif

#define MSIP 1.47334
#define BSIP 22.6746
#define MMGS -3.08542
#define MBGS 15.9161
#define BMGS -6.15469
#define BBGS 31.0333

/*the amount of over solve needed for this inner loop*/
#define NEWT 1.0

extern do_setup();

/* Backward Euler for generic conservative chemical time terms and
 * algebraic pseudo species. Unsupported custom transient terms are refused
 * by the controller before any trial mutation. */
int soldif_be(double dt, float temp, double **before, double **after,
              double **before_active, double **after_active,
              double *before_area, double *after_area)
{
    bval_compute(temp,after,dt);
    cs.old_del=0; cs.new_del=dt; cs.temp=temp;
    cs.old=before; cs.new=after;
    cs.olda=before_active; cs.newa=after_active;
    cs.oco=before_area; cs.nco=after_area; cs.type=BE;
    cs.sol=tm_sol; cs.nsol=tm_nsol; cs.elim=tm_elim;
    cs.il=tm_il; cs.l=tm_l; cs.loff=tm_loff;
    return soldif(TRUE,"BE",do_setup);
}


/************************************************************************
 *									*
 *	soldif_tr( del_t, temp, old, new, vel, iter ) -  This 		*
 *  routine gets the answer at a single time step.  It does this for 	*
 *  all the impurities and at the temp and velocity of the surface 	*
 *  indicated.	If the iter flag is true, it solves the linear system	*
 *  with an iterative method.						*
 *									*
 *  Original:	MEL	12/84						*
 *									*
 ************************************************************************/
soldif_tr( del_t, temp, old, new, olda, newa, oldarea, newarea, init)
double del_t;		/*the amount of time step*/
float temp;		/*the temperature of the solve*/
double **old;		/*the answer of the previous timestep*/
double **new;		/*the answer of this solution*/
double **olda;		/*the answer of the previous timestep*/
double **newa;		/*the answer of this solution*/
double *oldarea;	/*the area of the previous timestep*/
double *newarea;	/*the area of this solution*/
int init;		/*whether or not to reuse the last factorization*/
{

    /*compute the appropriate boundary condition elements*/
    bval_compute( temp, new, del_t );

    /*load the setup structure*/
    cs.old_del = 0.0;
    cs.new_del = del_t;
    cs.temp = temp;
    cs.old = old;
    cs.new = new;
    cs.olda = olda;
    cs.newa = newa;
    cs.nco = newarea;
    cs.oco = oldarea;
    cs.type = TR;
    cs.sol = tm_sol;
    cs.nsol = tm_nsol;
    cs.elim = tm_elim;
    cs.il = tm_il;
    cs.l = tm_l;
    cs.loff = tm_loff;

    /*now call the real solve routine*/
    return( soldif( init, "TR", do_setup ) );
}


/************************************************************************
 *									*
 *	soldif_bdf( del_t, temp, old, new, vel, iter ) -  This routine	*
 *  solves a single step of the diffusion equation using BDF2.  If a	*
 *  negative concentrations are encountered, the routine returns true.	*
 *									*
 *  Original:	MEL	1/86						*
 *									*
 ************************************************************************/
soldif_bdf( del_t, old_t, temp, new, mid, old, newa, mida, olda,
	    newarea, midarea, oldarea, init)
double del_t;		/*the amount of time step*/
double old_t;		/*last time step*/
float temp;		/*the temperature of the solve*/
double **new;		/*the answer of this solution*/
double **mid;		/*the answer of the previous timestep*/
double **old;		/*the answer of the really old timestep*/
double **newa;		/*the answer of this solution*/
double **mida;		/*the answer of the previous timestep*/
double **olda;		/*the answer of the really old timestep*/
double *newarea;	/*the area of this solution*/
double *midarea;	/*the area of the previous timestep*/
double *oldarea;	/*the area of the really old timestep*/
int init;		/*whether or not a setup has been done*/
{

    /*compute the appropriate boundary condition elements*/
    bval_compute( temp, new, del_t );

    /*load the setup structure*/
    cs.new_del = del_t;
    cs.old_del = old_t;
    cs.temp = temp;
    cs.old = old;
    cs.mid = mid;
    cs.new = new;
    cs.olda = olda;
    cs.mida = mida;
    cs.newa = newa;
    cs.nco = newarea;
    cs.mco = midarea;
    cs.oco = oldarea;
    cs.type = BDF;
    cs.sol = tm_sol;
    cs.nsol = tm_nsol;
    cs.elim = tm_elim;
    cs.il = tm_il;
    cs.l = tm_l;
    cs.loff = tm_loff;

    /*now call the real solve routine*/
    return( soldif( init , "BDF", do_setup) );
}




/************************************************************************
 *									*
 *	steady_state( ) - This routine calculates the steady state 	*
 *  variables.								*
 *									*
 *  Original:	MEL	11/85						*
 *									*
 ************************************************************************/
steady_state( temp, nsol, sol, new, newarea )
float temp;
int nsol;
int *sol;
double **new;
double *newarea;
{
    register int i, k, sk;
    int elim[MAXIMP];

    if ( nsol == 0 ) return;

    for(i = 0; i < nsol; i++) {
	if ( soltoimp[sol[i]] != Psi )
	    elim[i] = FALSE;
	else
	    elim[i] = TRUE;
    }

    /*load up the solution variables*/
    for(k = 0; k < nsol; k++) {
	sk = sol[k];
	for(i = 0; i < nn; i++) new[sk][i] = nd[i]->sol[sk];
    }
    /*load up concentration variables for potential solution*/
    for(k = 0; k < n_imp; k++) {
	sk = soltoimp[k];
	switch( sk ) {
	case As :
	case B  :
	case Sb :
	case P  :
	case I  :
	case V  :
	case Psi:
	case iBe:
	case iMg:
	case iSe:
	case iSi:
	case iSn:
	case iGe:
	case iZn:
	case iC :
	case iG :
	    for(i = 0; i < nn; i++) new[k][i] = nd[i]->sol[k];
	    break;
	}
    }

    /*compute the appropriate boundary condition elements*/
    bval_compute( temp, new, 0.0 );

    /*load the setup structure*/
    cs.old_del = 0.0;
    cs.new_del = 0.0;
    cs.temp = temp;
    cs.new = new;
    cs.newa = newa;
    cs.nco = newarea;
    cs.type = SS;
    cs.sol = sol;
    cs.nsol = nsol;
    cs.elim = elim;
    cs.il = ss_il;
    cs.l = ss_l;
    cs.loff = ss_loff;

    /*now call the real solve routine*/
    k = soldif( TRUE , "SS", do_setup );

    /*save away the calculated values*/
    for(k = 0; k < nsol; k++) {
	sk = sol[k];
	for(i = 0; i < nn; i++) nd[i]->sol[sk] = new[sk][i];
    }

}



/************************************************************************
 *									*
 *	soldif() - This routine solves the equations as setup for which *
 *  ever method is appropriate.  It uses a newton loop and refactors 	*
 *  the jacobian when appropriate.					*
 *									*
 *  Original:	MEL	5/88						*
 *									*
 ************************************************************************/
soldif( init, label, do_setup )
int init;
char *label;
PTR_FNC do_setup;
{
    register int i, si, j;			/*standard flame about indices*/
    int converge = FALSE;		/*are we converged yet??*/
    int count, blk;				/*loop count*/
    struct tms before, after;		/*for the system time call*/
    double nm[MAXIMP], maxnorm; 	/*norms for the newton iteration*/
    float tset, tsol;			/*time for setup, time for solution*/
    double rhsnm[MAXIMP], rhs2, maxrhs;	/*two norm of the right hand side*/
    double scrhs, lstnm = 1.0e37;
    int negat = FALSE;
    int factor = TRUE;
    double t1, norm2();
    double absrhserr = (mode==ONED)?(1.0e5):(10);
    char *ans = "  %c%-4d    %-11.5g  %-6.4g  %-6.4g  %-7.5g  %4d\n";
    char *s1  = "   iter    newton       setup   solve   ln rhs   Blk   (%s)\n";
    char *s2  = "           error        time    time    norm\n";

    if (cs.nsol == 0 ) return(FALSE);

    /*print the header if need be*/
    if ( verbose >= V_NORMAL) { printf(s1, label); printf(s2); }

    /*compute the initial matrix terms*/
    times(&before);
    do_setup( rhsnm );
    for(rhs2=0.0, i = 0; i < cs.nsol; i++) rhs2 += rhsnm[cs.sol[i]];
    rhs2 = sqrt( rhs2 );
    if (!isfinite(rhs2) && getenv("SUPREM4_SOLVE_DIAG")) {
        static int reports=0;
        for(i=0;i<cs.nsol;i++)for(j=0;j<nn;j++) {
            si=cs.sol[i];
            if(!isfinite(rhs[si][j]) && reports++<24)
                fprintf(stderr,"[solve-diag] initial stage=%s imp=%d node=%d mat=%d rhs=%.17g value=%.17g area=%.17g dt=%.17g\n",
                        label,soltoimp[si],j,mat_nd(j),rhs[si][j],cs.new[si][j],cs.nco[j],cs.new_del);
        }
    }
    times(&after);
    tset = (after.tms_utime - before.tms_utime) / (HZ * 1.0);

    count = 0;
    while ( (! converge) && (!negat) ) {
	count++;

	times(&before);
	/*figure if we need to factor the matrix*/
        /* alway do LU factor, CPU is faster now. fix by gongding 2013.03.18 */

	/*solve it with the appropriate method*/
	solve_blocks(nn, cs.sol, cs.nsol, cs.elim, cs.il, cs.l, cs.loff, newb, factor );
	times(&after);
        tsol = (after.tms_utime - before.tms_utime) / (HZ * 1.0);

	/*compute the update dot products*/

	/*update the set of variables*/
	negat = update_vars( nm, cs.nsol, cs.sol, nn, cs.new, newb );

	if ( negat ) continue;

	/*compute the initial matrix terms*/
	times(&before);
	do_setup( rhsnm );
	times(&after);
        tset = (after.tms_utime - before.tms_utime) / (HZ * 1.0);

	lstnm = rhs2;
	for(scrhs = rhs2 = 0.0, maxrhs = maxnorm = -1.0, i = 0; i < cs.nsol; i++)  {
	    si = cs.sol[i];
	    if ( nm[si] > maxnorm ) maxnorm = nm[si];
	    if ( rhsnm[si] > maxrhs ) {blk = soltoimp[si]; maxrhs = rhsnm[si]; }
	    rhs2 += rhsnm[si];
	    for(j = 0; j < nn; j++) {
		t1 = rhs[si][j] / (cs.new[si][j] + 1.0);
		scrhs += t1 * t1;
	    }
	}
	rhs2 = sqrt( rhs2 );
	scrhs = sqrt( scrhs );

	if ( cs.type == SS )
	    converge = maxnorm < 1.0e-3;
	else
	    converge = (rhs2 < absrhserr) || (maxnorm*NEWT < 1.0e-5) ;

	/*check to make sure we are making a least a little progress*/
        if ( (count > 10) && (rhs2 > 0.999*lstnm) ) negat = TRUE;
        if ( ( rhs2 > 1.0e10*lstnm ) ) negat = TRUE;
	if ( isnan( rhs2 ) ) negat = TRUE;

	/*print out the loop data*/
	if ( verbose >= V_NORMAL) {
	    if ( rhs2 == 0.0 )
		printf(ans, ' ', count, maxnorm*NEWT, tset, tsol, -10.0, blk);
	    else
		printf(ans, ' ', count, maxnorm*NEWT, tset, tsol, log10(rhs2), blk);
	    fflush(stdout);
	}

    }
    if(negat && getenv("SUPREM4_SOLVE_DIAG"))
        fprintf(stderr,"[solve-diag] rejected stage=%s iterations=%d rhs=%.17g previous_rhs=%.17g dt=%.17g\n",
                label,count,rhs2,lstnm,cs.new_del);
    return( negat );
}





/************************************************************************
 *									*
 *	update_vars( ) - This routine takes the update vector and 	*
 *  corrects the current guess.  It also calculates the norm of the 	*
 *  the update.  If any updates take the concentration negative, the	*
 *  routine returns true.						*
 *									*
 *  Original:	MEL	11/85  (Older than this, subroutine date)	*
 *									*
 ************************************************************************/
update_vars( nm, nsol, sol, nv, vars, update )
double *nm;
int nsol;
int *sol;
int nv;
double **vars;
double **update;
{
    register int bi, si, j;		/*loop counters*/
    double tmp;
    int imp;
    int negat = FALSE;
    int reports = 0;

    for(bi = 0; bi < nsol; bi++) {
	si = sol[bi];
	nm[si] = 0.0;
	imp = soltoimp[si];

	if (imp != Psi ) {
	    for(j = 0; j < nv; j++) {

		/*calculate the value of the update*/
		tmp = vars[si][j] + update[si][j];
                if ( fabs(tmp) < 1.0 )
                {
                  tmp = 1.0;
                }

		/*check to see if the new conc will be negative*/
		if ( tmp < 0 ) {
		    if (getenv("SUPREM4_SOLVE_DIAG")) {
                        if(reports++<1) {
                        int ti, vi, t, p = pt_nd(j);
                        fprintf(stderr,
                            "[solve-diag] negative imp=%d node=%d mat=%d value=%.17g update=%.17g dt=%.17g area=%.17g\n",
                            imp,j,mat_nd(j),vars[si][j],update[si][j],cs.new_del,cs.nco[j]);
                        fprintf(stderr,"[solve-diag] point=%d x=%.17g y=%.17g\n",
                                p,cordinate(p,0),cordinate(p,1));
                        if (reports == 1) for(ti=0;ti<num_tri_nd(j);ti++) {
                            t=tri_nd(j,ti);
                            for(vi=0;vi<num_vert(t);vi++) {
                                p=pt_nd(vert_tri(t,vi));
                                fprintf(stderr,"[solve-diag] triangle=%d vertex=%d point=%d x=%.17g y=%.17g area=%.17g\n",
                                        t,vi,p,cordinate(p,0),cordinate(p,1),tri[t]->earea[vi]);
                            }
                        }
                        }
                    }
		    /*only sweat non pseudo vars*/
		    if ( ! IS_PSEUDO( imp ) ) negat++;

		    /*if the amount we are going negative is less than ABE*/
		    if ( vars[si][j] > ABE[imp] )
			vars[si][j] = LTE[imp] * ABE[imp];
		    else
			vars[si][j] = LTE[imp] * LTE[imp] * vars[si][j];

		    tmp = update[si][j] / (LTE[imp] * vars[si][j] + ABE[imp]);
		}
		else {
		    /*normal updating*/
		    vars[si][j] = tmp;
		    tmp = update[si][j] / (LTE[imp] * vars[si][j] + ABE[imp]);
		}
		/*norm calculation*/
		if ( tmp > nm[si] )
		    nm[si] = tmp;
		else if ( - tmp > nm[si])
		    nm[si] = - tmp;
	    }
	}
	else {
	    /*damp the godammn poisson equations*/
	    for(j = 0; j < nv; j++) {
		tmp = vars[si][j] + update[si][j];

		if ( tmp < -1.0 ) tmp = -1.0;
		if ( tmp > 1.0 ) tmp = 1.0;

		vars[si][j] = tmp;
		tmp = update[si][j] / (LTE[imp] * vars[si][j] + ABE[imp]);

		/*norm calculation*/
		if ( tmp > nm[si] )
		    nm[si] = tmp;
		else if ( - tmp > nm[si])
		    nm[si] = - tmp;
	    }
	}
    }
    return( negat );
}



