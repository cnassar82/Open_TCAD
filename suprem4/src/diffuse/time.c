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
/*   time.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:50:23 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "matrix.h"
#include "mesh_trial.h"

extern double timed_trbdf();
static double timed_backward_euler(double delta, float temp);
extern int soldif_be(double,float,double **,double **,double **,double **,double *,double *);
extern int time_val();


/************************************************************************
 *									*
 *	time_diff( delta, vel, temp, iter ) - This routine solves a 	*
 *  single timestep for the diffusion equation.  It solves at 		*
 *  temperature temp and velocity of the interface vector vel.  This 	*
 *  routine uses the Yeager timestep method.  It calculates the new 	*
 *  timestep and puts it into delta for the calling program to use in 	*
 *  synchronizing with oxide growth.  The main trick in this routine is *
 *  keeping the right vectors in the right locations.			*
 *									*
 *  Original:	MEL	1/86						*
 *									*
 ************************************************************************/
double timed_trbdf( delta, temp, compute )
double delta;
float temp;
int compute;	/*flag for computing the comparison*/
{
    double del_t, old_t;	/*temp holder of the timestep*/
    double Kt;		/*the constant in the error term*/
    int xv= imptosol[XVEL], yv = imptosol[YVEL];
    int dela = imptosol[DELA];
    register int i, j, k, sk;	/*indices, what else?*/
    double *oarea;
    double *marea;
    double *narea;
    double nm;
    double tsave;
    double gam, igam;
    int fact;
    struct mesh_trial *trial = NULL;
    const char *stage = "initial";

    {
        const char *scheme=getenv("SUPREM4_TIME_INTEGRATOR");
        if (scheme && !strcmp(scheme,"be"))
            return timed_backward_euler(delta,temp);
    }

    gam = 2.0 - sqrt( 2.0 );
    igam = 1.0 - gam;

    if ( tm_nsol == 0.0 ) return( 2.0 );

    if ( verbose >= V_NORMAL) printf("solving at update = %g\n", delta);

    oarea = salloc(double, (ned>nn)?(ned):(nn));
    marea = salloc(double, (ned>nn)?(ned):(nn));
    narea = salloc(double, (ned>nn)?(ned):(nn));

    del_t = delta * gam;
    tsave = total;

    /*compute the initial symbolic*/
    fact = upd_diff( 0.0, oarea, FALSE, TRUE );

    /* Reordering is complete. Save the mesh before any tentative motion
     * or edge flips, so a rejected step can restore the same geometry. */
    if (xv != -1) trial = mesh_trial_save();

    /*load the values*/
    for( j = 0; j < n_imp; j++ ) {
	if ( (!IS_DIFFUSE(soltoimp[j])) && (!IS_PSEUDO(soltoimp[j])) ) continue;
	if ( (IS_STEADY( soltoimp[j] )) && (!IS_LOCKED( soltoimp[j]) )) continue;
	for(i=0; i<nn; i++) new[j][i] = mid[j][i] = old[j][i] = nd[i]->sol[j];
    }

    /*solve the diffusion problem if neccesary*/

    /*compute the initial right hand side*/
    do_vel(delta,gam);
    init_rhs( temp, old, olda, oarea );

    /*set up for the first TR step*/
    for(k = 0; k < tm_nsol; k++)  {
	sk = tm_sol[k];
	for(i = 0; i < nn; i++) oldb[sk][i] = newb[sk][i];
    }

    /*solve TR step to del_t*/
    stage = "TR motion";
    if ( xv != -1 ) {
	/*update the grid locations*/
	fact = upd_diff( del_t, marea, FALSE, FALSE );
        if (fact < 0) goto failed;
	do_vel(delta,gam);
	stage = "TR solve";
	if (soldif_tr(del_t,temp,old,mid,olda,mida,oarea,marea,TRUE)) {
	    goto failed;
	}
    }
    else {
	if (soldif_tr(del_t,temp,old,mid,olda,mida,oarea,oarea,TRUE)) {
	    goto failed;
	}
    }

    /*set up for the BDF step*/
    total += del_t;
    old_t = del_t;
    del_t = igam * delta;
    for(k = 0; k < tm_nsol; k++) {
	sk = tm_sol[k];
	for(i = 0; i < nn; i++) new[sk][i] = mid[sk][i];
    }

    /*now solve BDF2 for (1.0 - gam) * delta*/
    stage = "BDF motion";
    if ( xv != -1 ) {
	fact = upd_diff( del_t, narea, FALSE, FALSE );
        if (fact < 0) goto failed;
	do_vel(delta,gam);
	stage = "BDF solve";
	if (soldif_bdf(del_t,old_t,temp,new,mid,old,newa,mida,olda,
			narea, marea, oarea, fact) ) {
	    goto failed;
	}
    }
    else {
	if (soldif_bdf(del_t,old_t,temp,new,mid,old,newa,mida,olda,
			oarea, oarea, oarea, FALSE) ) {
	    goto failed;
	}
    }

    if ( compute ) {
	/*set up for a big step*/
	for(k = 0; k < tm_nsol; k++) {
	    sk = tm_sol[k];
	    for(i = 0; i < nn; i++) mid[sk][i] = new[sk][i];
	}
	if ( xv != -1 ) {
	    for(i = 0; i < nn; i++) {
		mid[xv][i] = new[xv][i];
		mid[yv][i] = new[yv][i];
		mid[dela][i] = new[dela][i];
	    }
	}

	/*Solve the big step*/
	stage = "comparison solve";
	if ( xv != -1 ) {
	    /*update the grid locations*/
	    i = soldif_tr(delta,temp,old,mid,olda,mida,oarea,narea,TRUE);
	}
	else {
	    i = soldif_tr(delta,temp,old,mid,olda,mida,oarea,oarea,TRUE);
	}
        /* A failed comparison solve is not an error estimate. Reject the
         * entire trial before publishing fields or recommending a timestep. */
        if (i) goto failed;

	/*now compute the norm of error for the new timestep*/
	for(nm = 0.0, k = 0; k < tm_nsol; k++) {
	    sk = tm_sol[k];
	    for(i = 0; i < nn; i++) {
		if ( oarea[i] != 0.0 ) {
		    Kt = (new[sk][i] - mid[sk][i]) /
			 (LTE[soltoimp[sk]] * new[sk][i] + ABE[soltoimp[sk]]);
                    if (!isfinite(Kt)) goto failed;
		    nm += Kt * Kt;
		}
	    }
	}

	total = tsave;
        if (!isfinite(nm)) goto failed;
	nm = sqrt(nm / ((double)tm_nsol * nn));

	Kt = (nm > 0) ? exp( log( 60.0 / (9.0 * nm) ) / 3.0 ) : MAXFLOAT;
    }
    else {
	Kt = 1.0;
	for(nm = 0.0, k = 0; k < tm_nsol; k++) {
	    sk = tm_sol[k];
	    if ( !IS_LOCKED(soltoimp[sk]) )
		for(i = 0; i < nn; i++) nd[i]->sol[sk] = new[sk][i];
	}
    }
    /*save the values, and continue*/
    for( j = 0; j < n_imp; j++ ) {
	if ( (!IS_DIFFUSE(soltoimp[j])) && (!IS_PSEUDO(soltoimp[j])) ) continue;
	if ( (IS_STEADY( soltoimp[j] )) && (!IS_LOCKED( soltoimp[j]) )) continue;
	for(i=0; i<nn; i++) nd[i]->sol[j] = new[j][i];
    }

    /*total the areas up*/
    if ((verbose >= V_CHAT) && (xv != -1)) tot_areas( oarea, narea );

    /*save away the active concentrations*/
    for ( j = 0; j < n_imp; j++ ) {
	switch ( soltoimp[j] ) {
	case Asa   : i = As;   break;
	case Ba    : i = B;    break;
	case Sba   : i = Sb;   break;
	case Pa    : i = P;    break;
	case iBea  : i = iBe;  break;
	case iMga  : i = iMg;  break;
	case iSea  : i = iSe;  break;
	case iSia  : i = iSi;  break;
	case iSna  : i = iSn;  break;
	case iGea  : i = iGe;  break;
	case iZna  : i = iZn;  break;
	case iCa   : i = iC;   break;
	case iGa   : i = iG;   break;
	default    : i = 0;    break;
	}
	if (i) {
	    sk = imptosol[i];
	    for(i = 0; i < nn; i++) set_sol_nd(i, j, newa[sk][i]);
	}
    }

    mesh_trial_discard(trial);
    free(oarea); free(marea); free(narea);
    return( Kt );

failed:
    if (getenv("SUPREM4_SOLVE_DIAG"))
        fprintf(stderr,"[time-diag] rejected stage=%s time=%.17g dt=%.17g\n",stage,tsave,delta);
    mesh_trial_restore(trial);
    mesh_trial_discard(trial);
    total = tsave;
    free(oarea); free(marea); free(narea);
    return(-1.0);
}


/*
 * update the velocities
 */
do_vel(delta,gam)
double delta;
double gam;
{
    register int i;
    int xvst= MAXIMP - 5, yvst = MAXIMP - 6;
    int xv= imptosol[XVEL], yv = imptosol[YVEL];
    double igam = 1.0 - gam;

    if (xv == -1) return;

    /*set up the old velocities*/
    for(i = 0; i < nn; i++) {
	old[xv][i] = nd[i]->sol[xvst];
	old[yv][i] = nd[i]->sol[yvst];
	mid[xv][i] = igam*nd[i]->sol[xvst]+gam*nd[i]->sol[xv];
	mid[yv][i] = igam*nd[i]->sol[yvst]+gam*nd[i]->sol[yv];
	new[xv][i] = nd[i]->sol[xv];
	new[yv][i] = nd[i]->sol[yv];
    }

    /*use the velocity and location data to get the delta areas*/
    moving_bound( delta, gam );
}


/* An independently estimated, first-order implicit integrator. The accepted
 * value is two half BE steps; one full BE step is used only for LTE. Rejected
 * trials restore topology and geometry and never publish concentrations. */
static double timed_backward_euler(double delta, float temp)
{
    int i,j,k,sk,imp,error_node=-1,error_sol=-1;
    int xv=imptosol[XVEL];
    double *oarea,*marea,*narea,*coarse;
    double error=0,den,diff,tsave=total,half=.5*delta;
    struct mesh_trial *trial=NULL;
    const char *stage="initial";
    if (!tm_nsol) return 2.0;
    for (k=0;k<tm_nsol;k++) {
        imp=soltoimp[tm_sol[k]];
        if (!IS_PSEUDO(imp) && impur[imp].time_val!=time_val)
            panic("BE requires generic chemical time terms or algebraic pseudo species");
    }
    oarea=salloc(double,(ned>nn)?ned:nn);
    marea=salloc(double,(ned>nn)?ned:nn);
    narea=salloc(double,(ned>nn)?ned:nn);
    (void)upd_diff(0.0,oarea,FALSE,TRUE);
    coarse=salloc(double,(size_t)tm_nsol*nn);
    if(xv!=-1)trial=mesh_trial_save();
    for(j=0;j<n_imp;j++) {
        if(!IS_DIFFUSE(soltoimp[j]) && !IS_PSEUDO(soltoimp[j]))continue;
        if(IS_STEADY(soltoimp[j]) && !IS_LOCKED(soltoimp[j]))continue;
        for(i=0;i<nn;i++)old[j][i]=mid[j][i]=new[j][i]=nd[i]->sol[j];
    }
    do_vel(delta,.5);
    init_rhs(temp,old,olda,oarea);
    stage="coarse motion";
    if(xv!=-1 && upd_diff(delta,narea,FALSE,FALSE)<0)goto reject;
    if(xv==-1)memcpy(narea,oarea,nn*sizeof(double));
    do_vel(delta,.5);
    stage="coarse BE";
    if(soldif_be(delta,temp,old,new,olda,newa,oarea,narea))goto reject;
    for(k=0;k<tm_nsol;k++)for(i=0;i<nn;i++)
        coarse[(size_t)k*nn+i]=new[tm_sol[k]][i];
    mesh_trial_restore(trial);
    if(xv!=-1 && upd_diff(0.0,oarea,FALSE,FALSE)<0)goto reject;
    for(k=0;k<tm_nsol;k++)for(i=0;i<nn;i++)
        mid[tm_sol[k]][i]=old[tm_sol[k]][i];
    do_vel(delta,.5);
    stage="first half motion";
    if(xv!=-1 && upd_diff(half,marea,FALSE,FALSE)<0)goto reject;
    if(xv==-1)memcpy(marea,oarea,nn*sizeof(double));
    do_vel(delta,.5);
    stage="first half BE";
    if(soldif_be(half,temp,old,mid,olda,mida,oarea,marea))goto reject;
    for(k=0;k<tm_nsol;k++)for(i=0;i<nn;i++)
        new[tm_sol[k]][i]=mid[tm_sol[k]][i];
    total=tsave+half;
    stage="second half motion";
    if(xv!=-1 && upd_diff(half,narea,FALSE,FALSE)<0)goto reject;
    if(xv==-1)memcpy(narea,oarea,nn*sizeof(double));
    do_vel(delta,.5);
    stage="second half BE";
    if(soldif_be(half,temp,mid,new,mida,newa,marea,narea))goto reject;
    for(k=0;k<tm_nsol;k++) {
        sk=tm_sol[k];imp=soltoimp[sk];
        for(i=0;i<nn;i++)if(oarea[i]>0 || narea[i]>0) {
            den=ABE[imp]+LTE[imp]*fmax(fabs(new[sk][i]),fabs(old[sk][i]));
            diff=fabs(new[sk][i]-coarse[(size_t)k*nn+i])/den;
            if(!isfinite(diff))goto reject;
            if(diff>error){error=diff;error_node=i;error_sol=k;}
        }
    }
    stage="BE error estimate";
    if(error>1.0)goto reject;
    total=tsave;
    /*save the values, and continue*/
    for( j = 0; j < n_imp; j++ ) {
	if ( (!IS_DIFFUSE(soltoimp[j])) && (!IS_PSEUDO(soltoimp[j])) ) continue;
	if ( (IS_STEADY( soltoimp[j] )) && (!IS_LOCKED( soltoimp[j]) )) continue;
	for(i=0; i<nn; i++) nd[i]->sol[j] = new[j][i];
    }

    /*total the areas up*/


    /*save away the active concentrations*/
    for ( j = 0; j < n_imp; j++ ) {
	switch ( soltoimp[j] ) {
	case Asa   : i = As;   break;
	case Ba    : i = B;    break;
	case Sba   : i = Sb;   break;
	case Pa    : i = P;    break;
	case iBea  : i = iBe;  break;
	case iMga  : i = iMg;  break;
	case iSea  : i = iSe;  break;
	case iSia  : i = iSi;  break;
	case iSna  : i = iSn;  break;
	case iGea  : i = iGe;  break;
	case iZna  : i = iZn;  break;
	case iCa   : i = iC;   break;
	case iGa   : i = iG;   break;
	default    : i = 0;    break;
	}
	if (i) {
	    sk = imptosol[i];
	    for(i = 0; i < nn; i++) set_sol_nd(i, j, newa[sk][i]);
	}
    }


    mesh_trial_discard(trial);
    free(oarea);free(marea);free(narea);free(coarse);
    return error>0 ? fmin(2.0,.9/sqrt(error)):2.0;
reject:
    if(getenv("SUPREM4_SOLVE_DIAG"))
        fprintf(stderr,"[be-diag] reject stage=%s time=%.17g dt=%.17g error=%.17g\n",stage,tsave,delta,error);
    if(getenv("SUPREM4_SOLVE_DIAG") && error_node>=0) {
        sk=tm_sol[error_sol];i=error_node;
        fprintf(stderr,"[be-diag] imp=%d point=%d old=%.17g fine=%.17g coarse=%.17g oldarea=%.17g newarea=%.17g\n",
                soltoimp[sk],pt_nd(i),old[sk][i],new[sk][i],coarse[(size_t)error_sol*nn+i],oarea[i],narea[i]);
    }
    mesh_trial_restore(trial);mesh_trial_discard(trial);total=tsave;
    free(oarea);free(marea);free(narea);free(coarse);
    return -1.0;
}
