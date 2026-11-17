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
/*   solve_time.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:50:18 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "material.h"	/* for the substrate orientation */
#include "matrix.h"
#include "sysdep.h"
#include "device.h"
#include "poly.h"
#include "defect.h"
#include "time_control.h"

#define DAMP 0.9

extern double timed_milne(), timed_trbdf(), timed_hybrid();
extern double start_step(), solstep(), do_oxide(), grid_dt();
extern double do_poly();
extern int time_val();
static double state_loss_step(float temp,double *area);
static double rezone_loss_step(float temp);


/************************************************************************
 *									*
 *	solve_time( time, temp, etc. ) - This routine solves the 	*
 *  diffusion equations.  It steps the time algorthim and handles the 	*
 *  coordination of oxide growth and diffusion.	  The variable time 	*
 *  is the length of this temperature step.				*
 *									*
 *  Original:	MEL	12/84						*
 *									*
 ************************************************************************/
solve_time( time, temp, itime, timemeth, dump, movie, cont )
float time;		/*the length of time to take this diffusion/oxidation*/
float temp;		/*the temperature at which this is performed*/
double itime;
int timemeth; 		/*time step method*/
int dump;	/*frequency of file writes*/
char *movie;		/*the movie command string*/
int cont;
{
    double t;			/*the current time, total*/
    static double del_df;	/*the next delta for the diffusion*/
    double last_ddf;		/* the last delta for the diffusion */
    double Ktd, Kto, Ktp;
    float time_form();
    struct tms before, after;
    int wrcnt = 1;
    register int i;
    double err = 1.0;
    double *oarea;

    /*set up the loop*/
    t = 0;
    UpdateSymbolic = 0;


    /*initialize the equilibrium defect concentrations*/
    err = start_step( temp, cont );
    printf("estimated first time step %e\n", err);

    /*the time dependent writes*/
    dump_data( total, dump, movie, &wrcnt );

    if ( itime != 0.0 ) del_df = itime;
    /* A structure restart can begin with `continue` before this process has
     * taken any adaptive step. Its static history is then still zero. */
    if(!isfinite(del_df) || del_df<=0)del_df=methdata.init_time;
    /* A restart after implantation/deposition can contain a fast interface
     * relaxation. The requested first step is an upper bound, not a reason
     * to discard the inventory-based startup estimate. */
    if(isfinite(err) && err>0 && err<del_df)del_df=err;
    last_ddf = del_df;

    /*loop until t = the time at the end of the ride*/
    while ( t < time ) {

	if ( (t + del_df) > time) del_df = time - t;

	/*Print some sort of header for the timestep*/
	if (verbose >= V_NORMAL)
	    printf("--------------------\nelapsed time = %g\n", total);
	else
	    printf("Solving %12g + %12g = %12g, %8g%%, np %d\n",
		   total,del_df, total+del_df, 100*del_df/last_ddf, np);

	times(&before);

	/*The oxidation step.  It may cut back the time step */
	if (oxidizing) {
	    GridSave();
	    Kto = do_oxide(temp,sub_ornt,(gas_type==GAS_WETO2),&del_df);
	}
	else
	    Kto = time - t + 1;

	if (grain_growing) {
            GridSave();
            Ktp = do_poly(temp,&del_df);
        }
        else
            Ktp = time - t + 1;

	/*solve the time step and get the increase*/
	/* Choose the more conservative of the oxide and diffusion increases */
        {
            const char *control=getenv("SUPREM4_POSITIVITY_STEP");
            if(control && !strcmp(control,"1") && UpdateSymbolic) {
                double bound=rezone_loss_step(temp);
                if(isfinite(bound) && bound>0 && bound<del_df) {
                    if(getenv("SUPREM4_SOLVE_DIAG"))
                        fprintf(stderr,"[time-bound] time=%.17g requested=%.17g depletion=%.17g\n",total,del_df,bound);
                    del_df=bound;
                }
            }
        }
	Ktd = solstep( del_df, temp, timemeth, TRUE );

	total += del_df;
	t += del_df;

	last_ddf = del_df;
        del_df = min(Ktd, min(Kto, Ktp));

	/* But if given a formula, take that instead */
	if( timemeth == FORM)
	    del_df = time_form( methdata.dtform, total, del_df);

	if( verbose >= V_NORMAL)
	    printf("    new timestep is %g increase %g [Ktd %.2e Kto %.2e]\n", del_df, del_df/last_ddf, Ktd/last_ddf, Kto/last_ddf);

	/* if there is polylayer, calculate grain size */
	if (n_poly>0)
	    for (i=0;i<n_poly;i++){
		poly_gs(i,temp,last_ddf);
		grain_gs(i,temp,last_ddf);
	    }


	/* now remove the nodes we wanted to kill but kept around  */
	/* so their dose would redistribute */
	if (oxidizing || grain_growing) {
	    GridSave();
	    grid_sub( last_ddf );
	}

	dump_data( total, dump, movie, &wrcnt );

	times(&after);
	print_time("Time for step", &before, &after);

	/* Flush those buffers, you sucker! */
	flushbuf();
    }

}

flushbuf() {fflush( stdout); fflush( stderr); }



/************************************************************************
 *									*
 *	do_oxide() - This routine manages the velocities and the moving	*
 *  boundary set up.							*
 *								 	*
 ************************************************************************/
double do_oxide( temp, ornt, oxhow, dt )
float temp;
int ornt;
int oxhow;
double *dt;
{
    int xv=imptosol[XVEL];
    int yv=imptosol[YVEL];
    int xvst= MAXIMP - 5, yvst = MAXIMP - 6;
    register int i;

    /*just a quick check...*/
    if ( yvst < n_imp ) panic("stack overflowed heap");

    /*save away the old velocities*/
    for(i = 0; i < nn; i++) {
	nd[i]->sol[xvst] = nd[i]->sol[xv];
	nd[i]->sol[yvst] = nd[i]->sol[yv];
    }

    /*solve the growth*/
    oxgrow(temp, ornt, oxhow, dt );

    return( grid_dt ());
}



/************************************************************************
 *									*
 *	solstep() - This was made a separate routine so it could be	*
 *  recursive.  In the case of negative concentrations, the routine	*
 *  splits the timestep in half and takes two steps.			*
 *									*
 *  Original:	MEL	7/87						*
 *									*
 ************************************************************************/
double solstep( del_df, temp, timemeth, compute )
double del_df;
float temp;
int timemeth;
int compute;	/*flag for computation of the split variable*/
{
    double Kt;
    static unsigned int trial_depth=0, trial_count=0;
    int xv = imptosol[XVEL];
    int yv = imptosol[YVEL];
    register int i,j,sj;
    double tsave;
    struct retry_velocity {
        struct nd_str *node;
        double start[2], end[2];
    } *history=NULL;
    int history_count=0;

    /* Bound the existing recursive retry tree. An unresolvable moving mesh
     * is a reported failed run, with a valid checkpoint, never infinite
     * subdivision or a published successful process stage. */
    if (!trial_depth) trial_count=0;
    if (trial_depth>=16 || ++trial_count>128 || del_df<=0 ||
        !isfinite(del_df) || total+del_df==total) {
        fprintf(stderr,"ERROR: diffusion retry budget exhausted at time=%.17g dt=%.17g "
                "depth=%u trials=%u; writing last valid state to step_failure.str\n",
                total,del_df,trial_depth,trial_count);
        if (ig2_write("step_failure.str",0,1.0)<0)
            fprintf(stderr,"ERROR: could not write failed-step checkpoint\n");
        exit(EXIT_FAILURE);
    }
    trial_depth++;
    tsave = total;

    if (n_imp != 0) {
	switch( timemeth ) {
	case TRBDF  :   Kt = timed_trbdf( del_df, temp, compute );
			Kt = (Kt > 1.0)?(DAMP * (Kt - 1.0) + 1.0): Kt;
			if ( (del_df > 1.0) && (Kt < MAXFLOAT / del_df) )
			    Kt = Kt * del_df;
			else if ( del_df <= 1.0 )
			    Kt = Kt * del_df;
			else
			    Kt = MAXFLOAT;
		        break;
	case FORM   :   Kt = timed_trbdf( del_df, temp, FALSE );
		        break;
	}
    }

    /*if Kt is negative split in half and repeat*/
    if ( Kt < 0.0 ) {
	/*back the points up*/
	del_df *= 0.5;

        /* timed_trbdf/do_vel reload velocity history from nd[].sol.
         * Keep it there, not in scratch arrays which the child overwrites.
         * Node objects remain stable when min-fill permutes their indices. */
        if (xv!=-1) {
            history_count=nn;
            history=salloc(struct retry_velocity,history_count);
            for(i=0;i<history_count;i++) {
                history[i].node=nd[i];
                history[i].start[0]=nd[i]->sol[MAXIMP-5];
                history[i].start[1]=nd[i]->sol[MAXIMP-6];
                history[i].end[0]=nd[i]->sol[xv];
                history[i].end[1]=nd[i]->sol[yv];
                nd[i]->sol[xv]=.5*(history[i].start[0]+history[i].end[0]);
                nd[i]->sol[yv]=.5*(history[i].start[1]+history[i].end[1]);
            }
        }

	/*solve first half of step*/
	Kt = solstep( del_df, temp, timemeth, FALSE );

        if (xv!=-1) for(i=0;i<history_count;i++) {
            struct nd_str *node=history[i].node;
            node->sol[MAXIMP-5]=.5*(history[i].start[0]+history[i].end[0]);
            node->sol[MAXIMP-6]=.5*(history[i].start[1]+history[i].end[1]);
            node->sol[xv]=history[i].end[0];
            node->sol[yv]=history[i].end[1];
        }

	total += del_df;

	/*reset the solution space*/
	for( j = 0; j < tm_nsol; j++ ) {
	    sj = tm_sol[j];
	    for(i = 0; i < nn; i++) new[sj][i] = mid[sj][i] = old[sj][i] = nd[i]->sol[sj];
	}

	/*solve second half of step*/
	Kt = solstep( del_df, temp, timemeth, compute );
	total -= del_df;

        if (xv!=-1) for(i=0;i<history_count;i++) {
            struct nd_str *node=history[i].node;
            node->sol[MAXIMP-5]=history[i].start[0];
            node->sol[MAXIMP-6]=history[i].start[1];
            node->sol[xv]=history[i].end[0];
            node->sol[yv]=history[i].end[1];
        }
        free(history);

    }
    total = tsave;
    trial_depth--;
    return( Kt );
}




/************************************************************************
 *									*
 *	dump_data() - Write out stuff at each time step.		*
 *									*
 ************************************************************************/
dump_data( total, dump, movie, cnt )
double total;
int dump;
char *movie;
int *cnt;
{
    char barfname[20];

    *cnt += 1;
    if (*cnt == dump) {
	sprintf(barfname, "s%g", total);
	(void)ig2_write (barfname, /*flip*/ 0, /* scale */ 1.0);
	*cnt = 0;
    }

    if ( movie ) make_movie( movie, total );
}





/************************************************************************
 *									*
 *	time_form()	- Let the user specify a time step.		*
 *									*
 ************************************************************************/
#include "expr.h"

float time_form( expr, total, dt)
    char *expr;
    float total;
    float dt;
{
    float string_to_real();
    double gdt, val;
    char *s, *index();

    /*Calculate the grid-motion-limited time step*/
    if( (s=index( expr, 'g')) != 0 && s[1]=='d' && s[2]=='t') {
	gdt = grid_dt();
	if( gdt != 0)
	    fmacro( "gdt", grid_dt(), "%e");
	else {
	    fprintf( stderr, "*** Sorry, gdt is zero, sticking with timestep = %g\n", dt);
	    return( dt);
	}
    }

    /*Calculate a timestep as an expression using the total time and gdt */
    fmacro( "t", total, "%e");
    fmacro( "dt", dt, "%e");
    val = string_to_real( expr, -1.0);
    if (val < 0) {
	fprintf(stderr, "problems with the formula %s", expr);
	val = dt;	/*constant*/
    }
    umacro( "t");
    umacro( "dt");
    umacro( "gdt");
    return( val);
}


/************************************************************************
 *									*
 *	start_step( temp ) - This routine estimates an initial timestep.*
 *									*
 *  Original:	MEL	10/86						*
 *									*
 ************************************************************************/
double start_step( temp, cont )
float temp;
int cont;
{
    double *oarea;
    register int j;
    register int i;
    register int cimp;
    int is;
    double del;


    /*compute the total area for use in calculations*/
    oarea = salloc(double, (ned>nn)?ned:nn);
    sum_coeff(oarea);

    /*set up old to compute right hand side*/
    for( j = 0; j < n_imp; j++ )
	for( i = 0; i < nn; i++) old[j][i] = nd[i]->sol[j];
    init_rhs(temp, old, olda, oarea);

    /*set up nodes... (this gets the active concentrations initialized)*/
    for( j = 0; j < n_imp; j++ ) {
	switch( (i = soltoimp[j]) ) {
	case Pa   : cimp = P  ; break;
	case Ba   : cimp = B  ; break;
	case Sba  : cimp = Sb ; break;
	case Asa  : cimp = As ; break;
	case iBea : cimp = iBe; break;
	case iMga : cimp = iMg; break;
	case iSea : cimp = iSe; break;
	case iSia : cimp = iSi; break;
	case iSna : cimp = iSn; break;
	case iGea : cimp = iGe; break;
	case iZna : cimp = iZn; break;
	case iCa  : cimp = iC ; break;
	case iGa  : cimp = iG ; break;
	default: cimp = -1; break;
	}
	if (cimp > 0)
	    for( i = 0; i < nn; i++) nd[i]->sol[j] = olda[imptosol[cimp]][i];
    }

    /*this gets the potential initialized..*/
    init_dev(temp);

    /*compute initial values of the steady state variables*/
    steady_state( temp, ss_nsol, ss_sol, new, oarea );

    /*initialize, if needed, the defect concentrations*/
    if ( !damage_read && !cont ) {
	is = imptosol[I];
	get_Cstar( I, nn, temp, noni, new[is], dequ[is] );
	for(i = 0; i < nn; i++) nd[i]->sol[is] = new[is][i];
	is = imptosol[V];
	get_Cstar( V, nn, temp, noni, new[is], dequ[is] );
	for(i = 0; i < nn; i++) nd[i]->sol[is] = new[is][i];
    }
    is = imptosol[T];
    if(is != -1 && !cont )
	for(i = 0; i < nn; i++) nd[i]->sol[is] = Tfrac(nd[i]->mater);

    /*if oxidizing compute some quicky oxide velocities*/
    if (oxidizing) {
	/* Don't use a fancy model if we're just looking for rough calc */
	int save_om = oxide_model;
	if (oxide_model == O_ELASTIC ||
	    oxide_model == O_VISCOUS ||
	    oxide_model == O_VISCOEL)
	    oxide_model = O_VERTICAL;
	oxide_vel( temp, sub_ornt, (gas_type == GAS_WETO2), 1e-32);
	oxide_model = save_om;
    }

    /* Initialization above changes potentials and defect fields. Estimate
     * from the resulting state, not the stale pre-initialization RHS. */
    del=state_loss_step(temp,oarea);
    free(oarea);
    return del;
}

/* A rezone can create a new fast chemical relaxation. Recompute a startup
 * estimate on the actual current mesh, without reinitializing defects,
 * potentials or physical concentrations. This is a depletion predictor,
 * not a proof of nonlinear TR-BDF2 positivity; trial rejection remains. */
static double rezone_loss_step(float temp)
{
    double *area=salloc(double,(ned>nn)?ned:nn),bound;
    upd_diff(0.0,area,FALSE,TRUE);
    bound=state_loss_step(temp,area);
    free(area);return bound;
}

static double state_loss_step(float temp,double *oarea)
{
    int j,i,is;double del,tmp;
    for(j=0;j<n_imp;j++)for(i=0;i<nn;i++)old[j][i]=nd[i]->sol[j];
    init_rhs(temp,old,olda,oarea);
    /* Only generic conservative, non-pseudo transient equations have M*C
     * storage. Algebraic variables and custom storage laws cannot use this
     * estimate. No zero-area or zero-concentration row can force dt=0. */
    for(del = MAXFLOAT, i = 0; i < tm_nsol; i++) {
	is = tm_sol[i];
	if(IS_PSEUDO(soltoimp[is]) || IS_STEADY(soltoimp[is]) ||
           impur[soltoimp[is]].time_val!=time_val)continue;
	for(j = 0; j < nn; j++) {
	    tmp=diffusion_loss_step(old[is][j],oarea[j],newb[is][j]);
	    if(tmp<del)del=tmp;
	}
    }
    return( del );
}


/*-----------------grid_dt----------------------------------------------
 * The grid-motion-limited time step.
 * This should rarely limit the time step, the diffusion generally being
 * more critical.  The exception is when there is no diffusion.
 * Algorithm is very simple, just the oxide grid criterion / maximum vel.
 *----------------------------------------------------------------------*/
double grid_dt()
{
    int ix, iy, i, in;
    int reg1,reg2;
    double vel, maxVel, spacing;

    /*
     * Figure the fastest moving point
     */
    maxVel=0;
    if((ix=imptosol[XVEL]) >= 0 && (iy=imptosol[YVEL]) >= 0) {
	for( i = 0; i < nn; i++) {
	    vel = hypot( nd[ i]->sol[ ix], nd[ i]->sol[ iy]);
	    if( vel > maxVel)
		maxVel = vel;
	}
    }

    if( maxVel == 0) return(MAXFLOAT);


    /*
     * Figure the perpendicular spacing. Usually this is given
     * but sometimes we have to actually figure it out.
     */
    if (grain_growing) spacing  = grain_grid * grain_gdt;
    else spacing  = oxide_grid * oxide_gdt;
    if( spacing <= 0) {

        /* This one isn't supposed to happen */
        if (grain_gdt <= 0) {
            fprintf( stderr, "*** The parameter grain.gdt on the method card is zero, choosing 0.25!\n");
            grain_gdt = 0.25;
        }
        else if (oxide_gdt <= 0) {
            fprintf (stderr, "*** The parameter oxide.gdt on the method card is zero, choosing 0.25!\n");
            oxide_gdt = 0.25;
        }


	/* Look at interface nodes */
	for( i = 0; i < nn; i++)
	    if( nd[ i]->mater == SiO2)
		if( ( in = node_mat( i, Si)) >= 0) {
		    float SilToOx[2], d;
		    assert( 0 == local_normal( in, SiO2, SilToOx));
		    d = rough_perp( nd[ i]->pt, SiO2, SilToOx);
		    if( d > spacing)
			spacing = d;
		}
        /* Look at polyinterface nodes */
        for( i = 0; i < nn; i++)
            if (ask(nd[i],REG_BND) && (pl_interface(i,&reg1,&reg2)) >= 2) {
                float SilToOx[2], d, rough_perp2();
                assert( 0 == local_normal2( i, reg2, SilToOx));
                d = rough_perp( i, reg2, SilToOx);
                if( d > spacing)
                        spacing = d;
            }

        if (spacing <= 0)
            return( MAXFLOAT);  /* Huh?  Grid movement with no oxide? */
        else
            if (grain_growing) spacing *= grain_gdt / 2;
            else  spacing *= oxide_gdt / 2;

    }

    /* Max dt is distance / velocity */
    return( spacing/ maxVel);
}

/*-----------------int_dose---------------------------------------------
 *----------------------------------------------------------------------*/

int_dose()
{
    register int ie, i, j;
    register int n;
    double dose[MAXIMP];

    for(j = 0; j < n_imp; j++) dose[j] = 0.0;

    for(ie = 0; ie < ne; ie++) {

	for(i = 0; i < 3; i++) {
	    n = tri[ie]->nd[i];

	    for(j = 0; j < n_imp; j++)
		dose[j] += nd[n]->sol[j] * tri[ie]->earea[i];
	}
    }

    for(j = 0; j < n_imp; j++)
	printf("integrated dose of imp %d is %e\n", soltoimp[j], dose[j]);
}
