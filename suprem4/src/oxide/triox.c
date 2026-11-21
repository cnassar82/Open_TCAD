/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *************************************************************************/
/*   triox.c                Version 5.1     */
/*   Last Modification : 7/3/91  10:52:40 */

#define X 0
#define Y 1
#define D 2
#define XX 0
#define YY 1
#define XY 2
#define NEL 7			/* number of node in element: fix in FEconvert too */


#include <stdio.h>
#include "suprem_options.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "triangle7.h"


/************************************************************************
 *									*
 *	FEshape - 6 or 7-noded triangle shape functions			*
 *									*
 ************************************************************************/
triox_shape( s,t, Cord, xsj, shp)
    double s,t;		/* natural coordinates*/
    double Cord[NEL][2];	/* x,y coords of all nodes*/
    double *xsj;	/* transformation jacobian determinant */
    double shp[NEL][3];	/* shp[i][2] = shape function, node i;
			   shp[i][0/1] = x,y derivative of shp[i][2] */
{
    int j, k;
    double tmpx, tmpy;
    double xs00, xs01, xs10, xs11, sx00, sx01, sx10, sx11;
    double cord[NEL][2];

    triangle7_basis(s,t,shp);

    /* A finite element must use its actual geometry, including obtuse
     * angles. Moving vertices only for assembly breaks the affine patch
     * test and disagrees with the boundary integrals. */
    for(k=0;k<NEL;k++)for(j=0;j<2;j++)cord[k][j]=Cord[k][j];

    /*The jacobian and its inverse*/
    xs00 = 0; xs01 = 0; xs10 = 0; xs11 = 0;
    for (k = 0; k < NEL; k++) {
	xs00 += (cord[k][0]-cord[0][0]) * shp[k][0];
	xs01 += (cord[k][1]-cord[0][1]) * shp[k][0];
	xs10 += (cord[k][0]-cord[0][0]) * shp[k][1];
	xs11 += (cord[k][1]-cord[0][1]) * shp[k][1];
    }
    *xsj = xs00*xs11 - xs01*xs10;
    if (*xsj < 1e-30) {
	fprintf(stderr, "warning: triangle has non-positive jacobian\n");
    }

    sx00 = xs11/ *xsj;
    sx11 = xs00/ *xsj;
    sx01 = -xs01/ *xsj;
    sx10 = -xs10/ *xsj;

    /*Form global derivatives*/
    for (k=0; k < NEL; k++) {
	tmpx = shp[k][0]*sx00 + shp[k][1]*sx01;
	tmpy = shp[k][0]*sx10 + shp[k][1]*sx11;
	shp[k][0] = tmpx;
	shp[k][1] = tmpy;
    }

}

#include "global.h"		/* for SAT */
#include "constant.h"		
#include "FEgeom.h"		/* for MAXDOF */
#include "FEmath.h"		/* for lambda, noderiv */


/************************************************************************
 *									*
 *	FEstiff - return the mixed stiffness/diffusivity matrix for a   *
 *		  NEL-noded oxidation triangle 				*
 *									*
 ************************************************************************/
/* L2 projection of the dilatational strain onto elementwise P1.
 * This is the condensed pressure space for the P2+bubble velocity element.
 * Keep deviatoric strain fully integrated; do not introduce hourglass modes
 * by reducing the entire stiffness quadrature. Physical coefficients stay
 * unchanged. The formula uses int(bubble)/area=3/20 and the exact P1 mass.
 */
extern int noderiv;
static void triox_strain_basis(double s,double t,double xl[7][2],
                               double shp[7][3],double bx[7][3],double by[7][3])
{
    double dx1=xl[1][0]-xl[0][0],dy1=xl[1][1]-xl[0][1];
    double dx2=xl[2][0]-xl[0][0],dy2=xl[2][1]-xl[0][1];
    double det=dx1*dy2-dy1*dx2,projection[2],gradient[3][2],l[3],alpha,cx,cy;
    int i,k;
    projection[0]=projection[1]=0;
    if(suprem_option("SUPREM4_FE_BBAR_P1")) {
        l[0]=1-s-t;l[1]=s;l[2]=t;
        gradient[1][0]=dy2/det;gradient[1][1]=-dx2/det;
        gradient[2][0]=-dy1/det;gradient[2][1]=dx1/det;
        for(k=0;k<2;k++)gradient[0][k]=-gradient[1][k]-gradient[2][k];
        for(i=0;i<3;i++)for(k=0;k<2;k++)projection[k]-=1.8*l[i]*gradient[i][k];
        projection[0]-=shp[6][0];projection[1]-=shp[6][1];
    }
    for(i=0;i<7;i++) {
        alpha=i<3?1./9:i<6?-4./9:1.;
        cx=.5*alpha*projection[0];cy=.5*alpha*projection[1];
        bx[i][0]=shp[i][0]+cx;bx[i][1]=cx;bx[i][2]=shp[i][1];
        by[i][0]=cy;by[i][1]=shp[i][1]+cy;by[i][2]=shp[i][0];
    }
}
static int triox_projected_stiff(double stiff[21][21],double rhs[7][3],
                                double xl[7][2],double disp[7][3],int (*coeff)(),int mat)
{
    double shp[7][3],bx[7][3],by[7][3],B[14][3],DB[14][3],eps[3],sig[3],m[40];
    double jac,dv,fx,fy,source,value;int q,i,j,k,c,d;
    for(q=0;q<7;q++) {
        triox_shape(triangle7_q[q][0],triangle7_q[q][1],xl,&jac,shp);
        triox_strain_basis(triangle7_q[q][0],triangle7_q[q][1],xl,shp,bx,by);
        dv=jac*triangle7_q[q][2];fx=fy=0;eps[0]=eps[1]=eps[2]=0;
        for(i=0;i<7;i++) {
            for(k=0;k<3;k++){B[2*i][k]=bx[i][k];B[2*i+1][k]=by[i][k];eps[k]+=bx[i][k]*disp[i][0]+by[i][k]*disp[i][1];}
            fx+=shp[i][0]*disp[i][2];fy+=shp[i][1]*disp[i][2];
        }
        (*coeff)(m,mat,eps,sig);
        for(i=0;i<14;i++) {
            value=0;for(k=0;k<3;k++)value+=B[i][k]*sig[k];
            rhs[i/2][i%2]-=dv*value;
            if(!noderiv)for(k=0;k<3;k++) {
                DB[i][k]=0;for(c=0;c<3;c++)DB[i][k]+=m[10+3*k+c]*B[i][c];
            }
        }
        for(i=0;i<7;i++)rhs[i][2]-=dv*m[5]*m[20]*(shp[i][0]*fx+shp[i][1]*fy);
        if(noderiv)continue;
        for(i=0;i<14;i++)for(j=0;j<14;j++) {
            value=0;for(k=0;k<3;k++)value+=B[i][k]*DB[j][k];
            stiff[3*(i/2)+i%2][3*(j/2)+j%2]+=dv*value;
        }
        for(i=0;i<7;i++) {
            source=dv*m[5]*(shp[i][0]*fx+shp[i][1]*fy);
            for(j=0;j<14;j++) {
                value=0;for(k=0;k<3;k++)value+=m[21+k]*B[j][k];
                stiff[3*i+2][3*(j/2)+j%2]+=source*value;
            }
            for(j=0;j<7;j++)stiff[3*i+2][3*j+2]+=dv*m[5]*m[20]*(shp[i][0]*shp[j][0]+shp[i][1]*shp[j][1]);
        }
    }
    return 0;
}

triox_stiff( stiff, wrhs, xl, disp, coeff, mat)
     double stiff[3*NEL][3*NEL]; /* The stiffness matrix */
     double wrhs[NEL][3];	/* extra right hand side stuff */
     double xl[NEL][2];		/* Coordinates of nodes */
     double disp[NEL][3];	/* The current solution/displacement */
     int (*coeff)();		/* Fetches material data */
     int mat;			/* What material this triangle is */
{
    if(suprem_option("SUPREM4_FE_BBAR_P1"))return triox_projected_stiff(stiff,wrhs,xl,disp,coeff,mat);


    register double six, siy, sjx, sjy;
    double db00,db01,db10,db11,db20,db21, ddx,ddy;
    double shp[NEL][3];
    double xsj, dv, dDD, eps[3], sig[3], sxx, syy, sxy, exx, eyy, exy, fx, fy;
    double matco[40], *DT, d00,d01,d02,d10,d11,d12,d20,d21,d22;
    double fr, fxx, fyy, fxy, ddd;
    int l, i, j;

    /*for each integration point compute contribution to stiffness*/
    for (l=0; l < TRIANGLE7_NQ; l++) {

	/* Integrate the cubic-bubble stiffness without hourglass modes. */
	triox_shape(triangle7_q[l][0],triangle7_q[l][1],xl,&xsj,shp);
	dv = triangle7_q[l][2]*xsj;

	/* Current strains and fluxes */
	exx = 0; eyy = 0; exy = 0; fx = 0; fy = 0;
	for (j = 0; j < NEL; j++) {
	    sjx = shp[j][X];
	    sjy = shp[j][Y];
	    exx += sjx * disp[ j][X];
	    eyy += sjy * disp[ j][Y];
	    exy += sjx * disp[ j][Y] + sjy * disp[j][X];
	    fx  += sjx * disp[ j][D];
	    fy  += sjy * disp[ j][D];
	}

	/* Calculate the current stresses and material coeffs*/
	eps[XX] = exx; eps[YY] = eyy; eps[XY] = exy;
	(*coeff) (matco, mat, eps, sig);
	sxx = sig[XX]; syy = sig[YY]; sxy = sig[XY];
	
	dDD = matco[ 5]*dv;
	fx *= dDD;
	fy *= dDD;
	fr  = matco[20];
	fxx = matco[21+XX]; fyy = matco[21+YY]; fxy = matco[21+XY];

	if (!noderiv) {
	    DT = matco +10;
	    d00 = DT[0]*dv; d01 = DT[1]*dv; d02 = DT[2]*dv;
	    d10 = DT[3]*dv; d11 = DT[4]*dv; d12 = DT[5]*dv;
	    d20 = DT[6]*dv; d21 = DT[7]*dv; d22 = DT[8]*dv;
	}

	/*Stiffness and right hand side*/
	for (j = 0; j < NEL; j++) {
	    sjx = shp[j][X]; sjy = shp[j][Y];

	    /*compute rhs*/
	    wrhs[ j][ X] -= dv * (sxx*sjx + sxy*sjy);
	    wrhs[ j][ Y] -= dv * (sxy*sjx + syy*sjy);
	    wrhs[ j][ D] -= fr * (fx* sjx + fy* sjy);

	    if (noderiv) continue;
	    
	    /*for each j node compute db = d*b*/
	    db00 = d00* sjx + d02* sjy;
	    db01 = d01* sjy + d02* sjx;
	    db10 = d10* sjx + d12* sjy;
	    db11 = d11* sjy + d12* sjx;
	    db20 = d20* sjx + d22* sjy;
	    db21 = d21* sjy + d22* sjx;
	    ddx = fr*dDD*sjx;
	    ddy = fr*dDD*sjy;
	    ddd = fx* sjx + fy*sjy;

	    /*for each i node compute s = bTdb*/
	    for (i=0; i < NEL; i++) {
		six = shp[i][X];
		siy = shp[i][Y];
		stiff[ 3*i+X][ 3*j+X] += six*db00 + siy*db20;
		stiff[ 3*i+X][ 3*j+Y] += six*db01 + siy*db21;
		stiff[ 3*i+Y][ 3*j+X] += six*db20 + siy*db10;
		stiff[ 3*i+Y][ 3*j+Y] += six*db21 + siy*db11;
		stiff[ 3*j+D][ 3*i+X] += ddd*(fxx*six + fxy*siy);
		stiff[ 3*j+D][ 3*i+Y] += ddd*(fxy*six + fyy*siy);
		stiff[ 3*j+D][ 3*i+D] += six*ddx  + siy*ddy;
	    }
	}
    }
    return(0);
}


/*
 *These offensive include's reflect the fact that the finite element
 *code doesn't pass around neighbor information which therefore must
 *be extracted out of the data structure.
 */
#include "constant.h"
#include "material.h"
static double simpson[3]={ 1/6., 1/6., 4/6.};
static int edges[3][3] = {{1,2,3},{2,0,4},{0,1,5}};
static double locg[3][2] = {{0.5,0.5},{0,0.5},{0.5,0}};

/* Experimental mixed-order scalar space. T prolongs the three corner
 * concentrations to the seven velocity nodes. Use T^T R and T^T J T,
 * with exact interpolation equations for the redundant scalar unknowns.
 * No concentration is clipped and no reaction coefficient is changed. */
static const double oxidant_prolong[7][3]={
 {1,0,0},{0,1,0},{0,0,1},{0,.5,.5},{.5,0,.5},{.5,.5,0},{1./3,1./3,1./3}
};
static int triox_scalar_map(int i,int *fixity,double *map)
{
 int j,node=i/3,component=i%3,reduced;
 reduced=(component==2 && node>=3 && suprem_option("SUPREM4_FE_OXIDANT_P1")) ||
         (component<2 && node>=3 && node<6 && fixity[i]==4 && suprem_option("SUPREM4_FE_LINEAR_TRACE"));
 for(j=0;j<7;j++)map[j]=reduced?(j<3?oxidant_prolong[node][j]:0):(j==node?1:0);
 return reduced;
}
void triox_project_oxidant(double *disp,double *original,int *fixity)
{
 int i,j;double map[7];
 for(i=0;i<21;i++)if(triox_scalar_map(i,fixity,map)) {
  disp[i]=0;for(j=0;j<7;j++)disp[i]+=map[j]*original[3*j+i%3];
 }
}
void triox_reduce_oxidant(double *stiff,double *rhs,double *original,int *fixity)
{
 double columns[21][21],result[21][21],out[21],map[21][7],constraint;
 int reduced[21],i,j,k;
 memset(columns,0,sizeof(columns));memset(result,0,sizeof(result));memset(out,0,sizeof(out));
 for(i=0;i<21;i++) {
  reduced[i]=triox_scalar_map(i,fixity,map[i]);
  for(k=0;k<7;k++)if(map[i][k])out[3*k+i%3]+=map[i][k]*rhs[i];
 }
 if(!noderiv) {
  for(i=0;i<21;i++)for(j=0;j<21;j++)for(k=0;k<7;k++)if(map[j][k])columns[i][3*k+j%3]+=stiff[21*i+j]*map[j][k];
  for(i=0;i<21;i++)for(j=0;j<21;j++)for(k=0;k<7;k++)if(map[i][k])result[3*k+i%3][j]+=map[i][k]*columns[i][j];
 }
 for(i=0;i<21;i++)if(reduced[i]) {
  constraint=original[i];for(j=0;j<7;j++)constraint-=map[i][j]*original[3*j+i%3];
  out[i]=-constraint;
  if(!noderiv) {
   result[i][i]=1;
   for(j=0;j<7;j++)if(map[i][j])result[i][3*j+i%3]=-map[i][j];
  }
 }
 memcpy(rhs,out,sizeof(out));if(!noderiv)memcpy(stiff,result,sizeof(result));
}
static double *oxidant_projected;
static int oxidant_capacity;

/* garbage storage is like so: */
/* normal [xy] ks [xx yy xy k0 vcv kr] ri ri+1 toc*/
/*         01      2  3  4  5  6   7   8  9    10*/
#define NORMAL 0
#define KS     2
#define KS0 3
#define VCV 4
#define KSR 5
#define RI  8
#define TOC 10			/* thin oxide correction */

/************************************************************************
 *									*
 *	triox_bc - boundary conditions at the silicon interface		*
 *									*
 ************************************************************************/
triox_bc( stiff, wrhs, xl, disp, coeff, mat, fixity, bc, ie, fixme)
     double stiff[3*NEL][3*NEL]; /* The stiffness matrix */
     double wrhs[NEL][3];	/* extra right hand side stuff */
     double xl[NEL][2];		/* Coordinates of nodes */
     double disp[NEL][3];	/* The current solution/displacement */
     int (*coeff)();		/* Fetches material data */
     int mat;			/* What material this triangle is */
     int fixity[NEL][3];	/* The nodal boundary codes */
     double bc[NEL][3];		/* Any nodal boundary float info */
     int ie;			/* because all this is still not enough */
     int fixme;			/* a variable of unsurpassed ugliness */
{
    int l,m,j,k, nbr, in, c, IsItf=0;
    double *ks, nx, ny, len, ks0, vcv, fr, ksv, cpl, *n, toc;
    double matco[40], xsj, dkde0, dkde1, dkde2, dkdu, dkdv, *DT, average_factor;
    static double shpg[NEL][NEL][3], shpx[NEL][3], sigx[3];
    double bx[7][3],by[7][3];
    static const double node_position[7][2]={{0,0},{1,0},{0,1},{.5,.5},{0,.5},{.5,0},{1./3,1./3}};
    static int do_sdep=1;

    /*
     * Surface tension term.
     * Note that STgamma must be set by FEoxbulk so bulk must be
     * done before boundary in FEassmb.
     */
    if( FEelt[ ie]->mat == SiO2) {
	for( j = 0; j < 3; j++) {
	    if( FEelt[ ie]->face[ j] == -1022) {
		double dx, dy;
		/* Two terms to add, one from each half side */
		for( k = 1; k<= 2; k++) {
		    l = (j+k)%3;
		    m=j+3;
		    dx = xl[ l][0] - xl[ m][0];
		    dy = xl[ l][1] - xl[ m][1];
		    if( fixity[ l][ X] != 1) wrhs[ l][ X] -= STgamma*( dx);
		    if( fixity[ l][ Y] != 1) wrhs[ l][ Y] -= STgamma*( dy);
		    if( fixity[ m][ X] != 1) wrhs[ m][ X] += STgamma*( dx);
		    if( fixity[ m][ Y] != 1) wrhs[ m][ Y] += STgamma*( dy);
		}
	    }
	}
    }

    /* Redo static terms? Yech... */
    if (fixme) TrioxBCinit();
    
    /* First clear existing terms.*/
    for( l = 0; l < NEL; l++) {
        if( fixity[ l][X] != 4 && fixity[ l][Y] != 4) continue;
	IsItf++;
	wrhs[ l][X] = 0;
	wrhs[ l][Y] = 0;

	for (m = 0; m < 3*NEL; m++) {
	    stiff[ 3*l+X][ m] = 0;
	    stiff[ 3*l+Y][ m] = 0;
	}
    }

    if (!IsItf) return(0);

    /* Differentiate the actual nodal strain traces. */
    if (!noderiv)
    for (l = 0; l < NEL; l++) {
        triox_shape(node_position[l][0],node_position[l][1],xl,&xsj,shpg[l]);
    }

    /*Loop through the nodes looking for the magic variables*/
    for( l = 0; l < NEL; l++) {
	if( fixity[ l][ X] != 4 && fixity[ l][ Y] != 4) continue;
	in = FEelt[ ie]->nd[ l];
	n =  &FEnd[ in]->garbage[NORMAL];
	ks = &FEnd[ in]->garbage[KS];

	/* The reaction law always uses the oxide constitutive stress of
         * the averaged strain, including at oxide/nitride junctions. */
        (*coeff)(matco, SiO2, FEnd[in]->eps, sigx);
		
	/*
	 * Surface velocity equation: v = ks * C/N1 * normal + toc
	 */
	ks0 = ks[KS0];
	vcv = ks[VCV];
	fr  = ks[KSR];
	ksv = ks0*vcv*fr;
	toc = ks[ TOC];

	if(mat == SiO2) {
	wrhs[ l][ X] = -(disp[l][X] - ksv * n[X] *disp[l][D] - toc * n[X]);
	wrhs[ l][ Y] = -(disp[l][Y] - ksv * n[Y] *disp[l][D] - toc * n[Y]);
	if (noderiv) continue;
	stiff[ 3*l+X][ 3*l+X] = 1;
	stiff[ 3*l+Y][ 3*l+Y] = 1;
	stiff[ 3*l+X][ 3*l+D] =     - ksv * n[X];
	stiff[ 3*l+Y][ 3*l+D] =     - ksv * n[Y];
        }
        if(noderiv)continue;

	if( do_sdep) {
	/* Stress dependence of growth rate */
	/* We know dk/dsig, now get dsig/deps */
        /* Residuals are repeated once per incident oxide element, but
         * ave_stress averages over ALL active elements. Its derivative
         * therefore needs every side with the same averaging denominator.
         * Restricting this term to oxide loses the nitride coupling and
         * overweights oxide contributions at a material junction. */
        assert(FEnd[in]->strain_count > 0);
        average_factor=FEnd[in]->garbage[RI+1]/FEnd[in]->strain_count;
	DT = matco+10;
		
	/* Form dk/eps */
	dkde0 = average_factor * (ks[XX]*DT[0] + ks[YY]*DT[3] + ks[XY]*DT[6]);
	dkde1 = average_factor * (ks[XX]*DT[1] + ks[YY]*DT[4] + ks[XY]*DT[7]);
	dkde2 = average_factor * (ks[XX]*DT[2] + ks[YY]*DT[5] + ks[XY]*DT[8]);

	/* and deps/du */
        /* Use the actual gradient trace of the cubic-bubble basis.
         * Linear extrapolation from three midpoints gives a nonzero
         * corner gradient to an interior bubble whose true trace is zero. */
        for(c=0;c<NEL;c++) {
            shpx[c][X]=shpg[l][c][X];shpx[c][Y]=shpg[l][c][Y];
        }

	
        triox_strain_basis(node_position[l][0],node_position[l][1],xl,shpg[l],bx,by);
	/* and build up drhs/du neighbor by neighbor */
	for (m = 0; m < NEL; m++) {
	    dkdu = dkde0*bx[m][0]+dkde1*bx[m][1]+dkde2*bx[m][2];
	    dkdv = dkde0*by[m][0]+dkde1*by[m][1]+dkde2*by[m][2];
	    stiff[ 3*l+X][ 3*m+X] += -vcv * n[X] * dkdu * disp[l][D];
	    stiff[ 3*l+X][ 3*m+Y] += -vcv * n[X] * dkdv * disp[l][D];
	    stiff[ 3*l+Y][ 3*m+X] += -vcv * n[Y] * dkdu * disp[l][D];
	    stiff[ 3*l+Y][ 3*m+Y] += -vcv * n[Y] * dkdv * disp[l][D];
	    stiff[ 3*l+D][ 3*m+X] +=  FEnd[ in]->garbage[RI] * dkdu;
	    stiff[ 3*l+D][ 3*m+Y] +=  FEnd[ in]->garbage[RI] * dkdv;
	}
	}
    }

    /* Radiative residual and its concentration derivative live on oxide. */
    if(mat != SiO2)return 0;

    /*Integrate the radiative boundary condition*/
    for (j = 0; j < 3; j++) {
	nbr = FEelt[ie]->face[ j];
	if (nbr >= 0 && FEelt[ nbr]->mat == Si) {
	    
	    /*Silicon edges - get the length of the edge*/
	    /*Sloppy because the triangles might be bent, but...*/
	    /*Normal points into oxide*/
	    nx = xl[(j+1)%3][1] - xl[(j+2)%3][1];
	    ny = xl[(j+2)%3][0] - xl[(j+1)%3][0];
	    len = hypot( nx, ny);

	    /*Use a Simpson's rule since we have sigma at the ends*/
	    /*Then take advantage of shp sub i(j)=delta sub ij*/
	    /* for each of the nodes on the edge */
	    for (k = 0; k < 3; k++) {
		l = edges[j][k];
		in = FEelt[ ie]->nd[ l];
		ks = &FEnd[ in]->garbage[KS];
		ks0 = ks[KS0];
		fr =  ks[KSR];
		cpl = ks0* fr* len* simpson[ k];
		stiff[ 3*l+D][ 3*l+D] += cpl;
		wrhs [ l][ D]         -= cpl * disp[ l][ D];
	    }
	}
    }
    return(0);
}
	    

/************************************************************************
 *									*
 *	FEepsg()	-  pointwise stress, strain			*
 *									*
 ************************************************************************/
triox_stress( mat, coeff, xl, ul, where, eps, sig, xsj)
    int mat;
    int (*coeff)();
    double xl[NEL][2];
    double ul[NEL][3];
    double *where;
    double eps[3];
    double sig[3];
    double *xsj;
{
    int j;
    double shp[NEL][3];
    double matco[40],bx[7][3],by[7][3];

    /* Evaluate strains */
    eps[XX] = 0; eps[YY] = 0; eps[XY] = 0;
    triox_shape( where[0], where[1], xl, xsj, shp);
    triox_strain_basis(where[0],where[1],xl,shp,bx,by);
    for (j=0; j < NEL; j++) {
	eps[ XX] += bx[j][0]*ul[j][X]+by[j][0]*ul[j][Y];
	eps[ YY] += bx[j][1]*ul[j][X]+by[j][1]*ul[j][Y];
	eps[ XY] += bx[j][2]*ul[j][X]+by[j][2]*ul[j][Y];
    }

    /* Call the strain to stress routine */
    (*coeff)( matco, mat, eps, sig);

}

/*-----------------TRIOX_NODAL_STRESS---------------------------------*
 * Evaluate the strains and stress at the actual element nodes.
 *--------------------------------------------------------------------*/
triox_nodal_stress( xl, ul, epsn, sign, coeff, mat, xsj)
    double xl[NEL][2];
    double ul[NEL][3];
    double epsn[NEL][3];
    double sign[NEL][3];
    int (*coeff)();
    int mat;
    double *xsj;
{
    static const double position[7][2]={{0,0},{1,0},{0,1},{.5,.5},{0,.5},{.5,0},{1./3,1./3}};
    int j;
    /* The displacement space is cubic, so its strain is quadratic.  A
     * three-point linear extrapolation cannot give its nodal traces. */
    for(j=0;j<NEL;j++) {
        double where[2];where[0]=position[j][0];where[1]=position[j][1];
        triox_stress(mat,coeff,xl,ul,where,epsn[j],sign[j],xsj);
    }
}



/*-----------------TrioxBCinit------------------------------------------
 * Set up nodally defined BC info.
 *----------------------------------------------------------------------*/
TrioxBCinit() 
{
    int in, k, l, nbr, j, ie;
    double matco[40], sig[3], nx, ny, len, cpl, *ks, *normal;
    FEelt_typ *pie; int *ip;

    if(suprem_option("SUPREM4_FE_OXIDANT_P1")) {
        int corner;
        if(oxidant_capacity<FEnn) {
            free(oxidant_projected);oxidant_projected=malloc(FEnn*sizeof(double));
            if(!oxidant_projected)abort();oxidant_capacity=FEnn;
        }
        for(in=0;in<FEnn;in++)oxidant_projected[in]=FEnd[in]->sol[D];
        for(ie=0;ie<FEne;ie++)if(FEdesc[FEelt[ie]->desc].stiff==triox_stiff) {
            for(j=3;j<7;j++) {
                in=FEelt[ie]->nd[j];oxidant_projected[in]=0;
                for(corner=0;corner<3;corner++)oxidant_projected[in]+=oxidant_prolong[j][corner]*FEnd[FEelt[ie]->nd[corner]]->sol[D];
            }
        }
    }

    /* Redundancy here: we only use the strains, not stress */
    ave_stress();
	
    /*
     * Compute reaction rate constants:
     * ks[0..2] is d ks/dstress, ks[3] is ks0, ks[4] is frac/N1
     * ks[5] is the reduction
     */
    for (in = 0; in < FEnn; in++) {
	if (FEnd[ in]->fixity[X] != 4 && FEnd[ in]->fixity[Y] != 4) continue;
	FE_oxbulk( matco, SiO2, FEnd[in]->eps, sig);
	ks = &FEnd[ in]->garbage[KS];
	normal = &FEnd[ in]->garbage[NORMAL];
	ksubs( ks, sig, normal, (double *)FEnd[in]->user );
    }

    /* Compute contributions to radiative integral */
    for (in = 0; in < FEnn; in++)
	FEnd[ in]->garbage[RI] = FEnd[ in]->garbage[RI+1] = 0;

    for (ie = 0; ie < FEne; ie++) {
	pie = FEelt[ ie];
	if (pie->mat != SiO2) continue;

	/*Bump the neighbor count of each node*/
	for (j = 0; j < NEL; j++)
	    FEnd[ pie->nd[ j]]->garbage[ RI+1] += 1;

	/*See if there are any contributions to make from this triangle*/
	for (j = 0; j < 3; j++) {
	    nbr = pie->face[ j];
	    if (nbr >= 0 && FEelt[ nbr]->mat == Si) {

		/*Silicon edges - get the length of the edge*/
		/*Sloppy because the triangles might be bent, but...*/
		/*Normal points into oxide*/
		ip = pie->nd;
		nx = FEnd[ ip[(j+1)%3]]->cord[1] - FEnd[ ip[(j+2)%3]]->cord[1];
		ny = FEnd[ ip[(j+2)%3]]->cord[0] - FEnd[ ip[(j+1)%3]]->cord[0];
		len = hypot( nx, ny);

		/* for each of the nodes on the edge */
		for (k = 0; k < 3; k++) {
		    l = edges[j][k];
		    in = pie->nd[ l];
		    cpl = len * simpson[k];
		    FEnd[in]->garbage[ RI] += cpl * (suprem_option("SUPREM4_FE_OXIDANT_P1")?oxidant_projected[in]:FEnd[in]->sol[D]);
		}
	    }
	}
    }

    /* Divide dependency by number of neighbor triangles */
    for (in = 0; in < FEnn; in++)
	if (FEnd[ in]->garbage[RI+1] != 0)
	    FEnd[ in]->garbage[ RI] /= FEnd[ in]->garbage[ RI+1];
}
	

/*-----------------Get an oxload of this--------------------------------
 * Apply load boundary conditions on a silicon-only problem.
 * Called by FEassmb during re-solution in visco_growth()
 *----------------------------------------------------------------------*/
oxload( stiff, wrhs, xl, disp, coeff, mat, fixity, bc, ie, fixme)
     double stiff[2*7][2*7];	/* The stiffness matrix */
     double wrhs[7][2];		/* extra right hand side stuff */
     double xl[7][2];		/* Coordinates of nodes */
     double disp[7][2];		/* The current solution/displacement */
     int (*coeff)();		/* Fetches material data */
     int mat;			/* What material this triangle is */
     int fixity[7][2];		/* The nodal boundary codes */
     double bc[7][2];		/* Any nodal boundary float info */
     int ie;			/* because all this is still not enough */
     int fixme;			/* a variable of unsurpassed ugliness */
{
    int j, k, l, in, nbr;
    double nx, ny, len, tx, ty;
    
    /* Look for the oxide interface */
    if( FEelt[ ie]->mat != Si) return(0);

    for( j = 0; j < 3; j++) {
	nbr = FEelt[ ie]->face[ j];
	if( nbr >= 0 && FEelt[ nbr]->mat == SiO2) {

	    /* Oh boy we got one */
	    /* Calculate the length of the side */
	    nx = xl[(j+1)%3][1] - xl[(j+2)%3][1];
	    ny = xl[(j+2)%3][0] - xl[(j+1)%3][0];
	    len = hypot( nx, ny);

	    /* load up the nodes on this face */
	    for( k = 0; k < 3; k++) {
		l = edges[j][k];
		in = FEelt[ ie]->nd[ l];
		tx = FEnd[ in]->garbage[ 2+X];
		ty = FEnd[ in]->garbage[ 2+Y];
		wrhs[ l][ X] += simpson[ k] * len * tx;
		wrhs[ l][ Y] += simpson[ k] * len * ty;
	    }
	}
    }
    return(0);
}

	    
/*#################### D E B U G ########################################*/

#ifdef debug

/*-----------------Text_Triox-------------------------------------------
 * Test the triangle shape functions.
 *----------------------------------------------------------------------*/
test_triox()
{
    int i, j, c;
    double rr, ss, sum, sumx, sumy, xl[2], xr[2], derror;
    double xsj, shp[NEL][3], shp2[NEL][3];

    static double garbage[NEL][2] = { {0,0},{1,1},{-1,1}};
    for (i = 3; i < 6; i++) for (j = 0; j < 2; j++)
	garbage[i][j] = 0.5*(garbage[(i+1)%3][j] + garbage[(i+2)%3][j]);
    for (j=0; j<2; j++)
	garbage[6][j] = (garbage[0][j] + garbage[1][j] + garbage[2][j])/3;

    for (rr = 0; rr < 1; rr+= 0.2) {
	for (ss = 0; ss < 1; ss+= 0.2) {
	    if (rr + ss > 1) continue;
	    triox_shape( rr, ss, garbage, &xsj, shp);

	    sum = sumx = sumy = 0;
	    for (i = 0; i < NEL; i++) {
		sum += shp[i][2];
		sumx += shp[i][X];
		sumy += shp[i][Y];
	    }
	    sum -= 1.0;
	    if (fabs(sum) >= 1e-10) printf("[%g %g] sum %g\n",rr,ss,sum);
	    if (fabs(sumx) >= 1e-10) printf("[%g %g] sumx %g\n",rr,ss,sumx);
	    if (fabs(sumy) >= 1e-10) printf("[%g %g] sumy %g\n",rr,ss,sumy);

	    /* See if the derivatives really are */
	    triox_shape( rr+1e-2, ss+1e-2, garbage, &xsj, shp2);
	    for (c = 0; c < 2; c++) {
		xl[c] = 0;
		xr[c] = 0;
		for (j = 0; j < NEL; j++) xl[c] += shp[j][2]*garbage[j][c];
		for (j = 0; j < NEL; j++) xr[c] += shp2[j][2]*garbage[j][c];
	    }
	    for (j = 0; j < NEL; j++) {
		derror = shp2[j][2] - shp[j][2] 
			 - (xr[0]-xl[0])*(shp2[j][0] + shp[j][0])/2
			     - (xr[1]-xl[1])*(shp2[j][1] + shp[j][1])/2;

		if( derror > 1e-4 || derror < -1e-4)
		    printf("derivative error \n");
	    }
	}
    }
}

/*----------------------------------------------------------------------
 * Test the stiffness matrix contributions for boundary elements.
 *----------------------------------------------------------------------*/

test_tbc( ie, iv, vx, vy)
    int ie;			/* Which triangle to work on */
    int iv;			/* Vertex of triangle to tweak */
    double vx, vy;		/* Velocities to tweak */
{
#   define MELDOF 60    
    static double xl[MELDOF], disp[MELDOF], bc[MELDOF];
    static int    eqns[MELDOF], fixity[MELDOF];
    double stiff[MELDOF*MELDOF], wrhs[MELDOF];
    double stiff1[MELDOF*MELDOF], wrhs1[MELDOF];
    int dummy;
    FEdesc_typ *desc; FEelt_typ *pie;
    int i, j, ndf, nel;

    pie = FEelt[ ie];
    desc = &FEdesc[ pie->desc];
    nel = desc->nel;
    ndf = nel*FEdf;

    /*clear out the arrays to be filled*/
    for (i = 0; i < ndf*ndf; i++) stiff[ i] = 0;
    for (i = 0; i < ndf; i++) wrhs[ i] = 0;
    for (i = 0; i < ndf*ndf; i++) stiff1[ i] = 0;
    for (i = 0; i < ndf; i++) wrhs1[ i] = 0;

    /*Change from global to local info*/
    localize( ie, xl, disp, eqns, fixity, bc);

    /*do element boundary conditions*/
    (*(desc->bound)) (stiff, wrhs, xl, disp, desc->coeff, pie->mat, fixity, bc, ie);

    
    disp[ FEdf*iv + X] += vx;
    disp[ FEdf*iv + Y] += vy;
    
    (*(desc->bound)) (stiff1, wrhs1, xl, disp, desc->coeff, pie->mat, fixity, bc, ie);

    for (i = 0; i < nel; i++){
	for (j =0; j < FEdf; j++) {
	    printf("%2d %2d %12.5g %12.5g\n",
		   i, j, wrhs1[FEdf*i+j] - wrhs[FEdf*i+j],
		   stiff[ndf*(FEdf*i+j) + FEdf*iv+X] * vx +
		   stiff[ndf*(FEdf*i+j) + FEdf*iv+Y] * vy);
	}
    }
}
		   

/*-----------------test_jacob-------------------------------------------
 * Test the Jacobian by touching a variable and comparing
 * Jdx to F(x) - F(x+dx).
 *----------------------------------------------------------------------*/

test_jacob( ia, aoff, in, iv, v, il, loff)
    int *ia, aoff;		/* matrix map */
    int in;			/* node to tweak */
    int iv;			/* dof to tweak */
    double v;			/* tweak value */
    int *il, loff;
{
    int neq, i, jn, jv, jf;
    double *a, *a1, *rhs, *rhs1, *deltav, *drhs; double *l;

    if (FEnn*FEdf != ia[0]-1) {
	printf("all confused"); return(0);
    }
    neq = ia[0]-1;
    
    a = salloc( double, ia[neq] + aoff);
    a1 = salloc( double, ia[neq] + aoff);
    rhs=salloc( double, neq);
    rhs1=salloc( double, neq);
    deltav = salloc( double, FEnn*FEdf);
    drhs   = salloc( double, FEnn*FEdf);
    
    FEassmb( ia, aoff, a, rhs);
    
    FEnd[ in]->sol[ iv] += v;
    FEassmb( ia, aoff, a1, rhs1);
    
    FEnd[ in]->sol[ iv] -= v;

    for (i = 0; i < neq; i++) deltav[i] = 0;
    deltav[ FEdf*in +iv] = v;
    mxv( neq, ia, aoff, 0.0, a, deltav, drhs);

    for (i = 0; i < neq; i++) {
	if (rhs1[i] != rhs[i] || drhs[i] != 0) {
	    jn=i/FEdf;
	    jv=i%FEdf;
	    printf("%d %d %c %12.5g %12.5g (%12.5g %12.5g)\n",
		   jn, jv, (jf=FEnd[jn]->fixity[jv])?('0'+jf):' ',
		   rhs[i] - rhs1[i], drhs[i], rhs1[i], rhs[i]);
	}
    }
#ifdef notdef
    /* A test of the condition number */
    l = salloc( double, il[neq] + loff);
    csr_Snumfac( neq, ia, aoff, a, il, loff, l);
    csr_Snumbac( neq, il, loff, l, drhs);
    for (i = 0; i < neq; i++)
	if (drhs[i] != deltav[i])
	    printf( "%d %d %12.5g %12.5g\n", i/FEdf, i%FEdf, drhs[i], deltav[i]);
    free( l);
#endif /*notdef*/
    free( a);
    free( a1);
    free( rhs);
    free( rhs1);
    free( deltav);
    free(drhs);

    return(0);
}

/* The convex debugger is feeble-minded. */
#ifdef notdef
#ifdef CONVEX
extern int *FDia;
extern int FDaoff;
extern int *FDil;
extern int FDloff;
test_jackass()
{
    int in, iv;
    double v;
    FILE *lu = fopen("LameDbx", "r");
    fscanf(lu, "%d %d %lf", &in, &iv, &v);
    fclose(lu);
    test_jacob( FDia, FDaoff, in, iv, v, FDil, FDloff);
}
#endif /*CONVEX*/
#endif /*NOTDEF*/

/*-----------------test_stress------------------------------------------
 * Look at direct and extrapolated stress measurements inside an element.
 *----------------------------------------------------------------------*/
#define r3i 0.57735027

test_stress( n)
    int n;
{
    static double xl[MELDOF], disp[MELDOF], bc[MELDOF];
    static int    eqns[MELDOF], fixity[MELDOF];
    double epsn[MAXSC*MAXNEL], sign[MAXSC*MAXNEL], x, y, xx, yy;
    FEdesc_typ *desc; FEelt_typ *pie;
    int i, in, j, k, nel; short *count, iv1, iv2;
    int ie, iv;
    double where[2], area[3], iarea[3], eps[3], sig[3], f;

    for (ie = 0; ie < FEne; ie++) {
	for (iv = 0; iv < 3; iv++) {

    pie = FEelt[ ie];
    desc = &FEdesc[ pie->desc];
    nel = desc->nel;

    if (pie->mat != SiO2 ||
	pie->face[iv] < 0 ||
	FEelt[pie->face[iv]]->mat != Si||
	n<=1)
	continue;

    iv1=(iv+1)%3;
    iv2=(iv+2)%3;
    localize( ie, xl, disp, eqns, fixity, bc);

#ifdef notdef
    /* The direct evaluation */
    for (f = 0; f < 1.0; f+=1.0/n) {
	area[ iv ] = 0;
	area[ iv1] = f;
	area[ iv2] = 1-f;
	where[0] = area[1];
	where[1] = area[2];
	triox_stress( pie->mat, desc->coeff, xl, disp, where, eps, sig);
	printf("[%8f,%8f] eps(%12.5g, %12.5g, %12.5g)\n", where[0], where[1],
	       eps[XX], eps[YY], eps[XY]);
    }
#endif /*notdef*/
    (*(desc->nodal_stress)) (xl, disp, epsn, sign, desc->coeff, pie->mat);
    /* The extrapolated evaluation */
    for (f = 0.5*(1-r3i); f<1; f += r3i) {
	area[ iv ] = 0;
	area[ iv1] = f;
	area[ iv2] = 1-f;
	x = (area[iv ]*xl[2*iv +X] + 
	     area[iv1]*xl[2*iv1+X] +
	     area[iv2]*xl[2*iv2+X]);
	y = (area[iv ]*xl[2*iv +Y] + 
	     area[iv1]*xl[2*iv1+Y] +
	     area[iv2]*xl[2*iv2+Y]);
	where[0] = area[1];
	where[1] = area[2];
	iarea[ iv ] = area[iv1] + area[iv2] - area[iv ];
	iarea[ iv1] = area[iv2] + area[iv ] - area[iv1];
	iarea[ iv2] = area[iv ] + area[iv1] - area[iv2];
	xx = (iarea[ iv ]*xl[ 2*(3+iv) +X] +
	      iarea[ iv1]*xl[ 2*(3+iv1) +X] +
	      iarea[ iv2]*xl[ 2*(3+iv2) +X]);
	yy = (iarea[ iv ]*xl[ 2*(3+iv) +Y] +
	      iarea[ iv1]*xl[ 2*(3+iv1) +Y] +
	      iarea[ iv2]*xl[ 2*(3+iv2) +Y]);
	for (j = XX; j <= XY; j++)
	    eps[j] = (iarea[ iv ]*epsn[ 3*(3+iv) +j] +
		      iarea[ iv1]*epsn[ 3*(3+iv1)+j] + 
		      iarea[ iv2]*epsn[ 3*(3+iv2)+j]);
	printf("[ %8f, %8f] eps( %12.5g, %12.5g, %12.5g)\n",
	       x, y, 
	       eps[XX], eps[YY], eps[XY]);
    }

}}
    return(0);
}

#ifdef CONVEX
test_pbc()
{
    extern double *FDrhs, ffcos();
    double dx, dy;
    int i, j, k1, k2, k;
    
    for( i =0; i < FEne; i++) {
	for (j = 0; j < 3; j++) {
	    if( FEelt[ i]->face[ j] == -1022) {
		k1 = FEelt[ i]->nd[ (j+1)%3];
		k2 = FEelt[ i]->nd[ (j+2)%3];
		k  = FEelt[ i]->nd[  j+3   ];
		dx = FEnd[ k1]->cord[X] - FEnd[ k2]->cord[X];
		dy = FEnd[ k1]->cord[Y] - FEnd[ k2]->cord[Y];
		printf( "(%g,%g) l(%g,%g,%g) c(%g,%g,%g) r(%g,%g,%g)\n",
		       dx, dy,
		       FDrhs[3*k1+X], FDrhs[3*k1+Y], ffcos(dx,dy,FDrhs[3*k1+X], FDrhs[3*k1+Y]),
		       FDrhs[3*k +X], FDrhs[3*k +Y], ffcos(dx,dy,FDrhs[3*k +X], FDrhs[3*k +Y]),
		       FDrhs[3*k2+X], FDrhs[3*k2+Y], ffcos(dx,dy,FDrhs[3*k2+X], FDrhs[3*k2+Y] ) 
		       );
	    }
	}
    }
}
double ffcos(ax,ay,bx,by)
    double ax,ay,bx,by;
{
    double t1, t2, t3;
    t1 = ax*bx+ay*by;
    t2 = ax*ax+ay*ay;
    t3 = bx*bx+by*by;
    if(t2==0 || t3==0) return(0);
    else return( t1*t1/(t2*t3));
}
/* seem to have lost these suckers */
decomp_(){;}
solve_(){;}
    
#endif CONVEX
#endif /*debug*/
