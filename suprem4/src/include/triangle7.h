#ifndef SUPREM_TRIANGLE7_H
#define SUPREM_TRIANGLE7_H
/* Seven positive Dunavant points, exact through degree five on the
 * reference triangle (weights sum to 1/2). A cubic bubble has quadratic
 * gradients, so its stiffness products require at least degree four.
 * The historical three-midpoint rule introduces extra zero-energy modes.
 */
#define TRIANGLE7_NQ 7
static const double triangle7_q[7][3]={
 {1.0/3,1.0/3,.1125},
 {.470142064105115,.470142064105115,.066197076394253},
 {.059715871789770,.470142064105115,.066197076394253},
 {.470142064105115,.059715871789770,.066197076394253},
 {.101286507323456,.101286507323456,.0629695902724135},
 {.797426985353087,.101286507323456,.0629695902724135},
 {.101286507323456,.797426985353087,.0629695902724135}
};
/* Keep the historical bubble basis normalization (9*l0*l1*l2). This
 * agrees with the old cached midpoint derivatives and elastic element.
 * Define it at every point, rather than only three cached locations. */
static void triangle7_basis(double s,double t,double shp[7][3])
{
    double l0=1-s-t,l1=s,l2=t;int j;
    shp[0][2]=l0*(2*l0-1);shp[1][2]=l1*(2*l1-1);shp[2][2]=l2*(2*l2-1);
    shp[3][2]=4*l1*l2;shp[4][2]=4*l2*l0;shp[5][2]=4*l0*l1;
    shp[0][0]=1-4*l0;shp[0][1]=1-4*l0;
    shp[1][0]=4*l1-1;shp[1][1]=0;
    shp[2][0]=0;shp[2][1]=4*l2-1;
    shp[3][0]=4*t;shp[3][1]=4*s;
    shp[4][0]=-4*t;shp[4][1]=4*(1-s-2*t);
    shp[5][0]=4*(1-2*s-t);shp[5][1]=-4*s;
    shp[6][2]=9*l0*l1*l2;shp[6][0]=9*l2*(l0-l1);shp[6][1]=9*l1*(l0-l2);
    for(j=0;j<3;j++) {
        shp[0][j]+=shp[6][j]/9;shp[1][j]+=shp[6][j]/9;shp[2][j]+=shp[6][j]/9;
        shp[3][j]-=4*shp[6][j]/9;shp[4][j]-=4*shp[6][j]/9;shp[5][j]-=4*shp[6][j]/9;
    }
}
#endif
