/* Consistent sampling of a vertical deformation on neighboring constraints.
 * Experimental: refine BOTH owners of an oxide/poly interface when a more
 * finely sampled free surface would otherwise cross it under the same map.
 * This is not loop deletion: no point, material fragment or chemical inventory
 * is removed. Coordinates are constructed once for both material owners.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "matrix.h"

struct vc_face {int t,f,a,b,type;float xy[2][2];double lo,hi;};
static int face_order(const void *aa,const void *bb)
{
    const struct vc_face *a=aa,*b=bb;
    return a->lo<b->lo ? -1:a->lo>b->lo;
}
static double vc_area(float *a,float *b,float *c)
{return .5*(((double)b[0]-a[0])*((double)c[1]-a[1])-((double)b[1]-a[1])*((double)c[0]-a[0]));}
#include "conservative_transfer.h"
#define vc_density conserved_density
#define vc_parent conserved_parent
#define vc_transfer conserved_transfer

static int vc_reject(int line)
{
    static int reports=0;
    if(getenv("SUPREM4_CONSTRAINT_TRACE") && reports++<24)
        fprintf(stderr,"[vertical-constraints] preparation refused at source line %d\n",line);
    return 0;
}

static int vc_split(struct vc_face *face,int sample)
{
    int owner[2],edge[2],nodes[2][4],region[2],side,i,j,k,t,a,b,newpoint,newnode;
    double oldw[2][4]={{0}},neww[2][4]={{0}},fields[2][4][MAXIMP],values[4];
    double fractions[2][4][MAXIMP];
    double ratio,area,wa,wb,velocity[2],start[2];float c[3]={0},v[3][3];
    struct nd_str *stable[2][3];struct pt_str *newpt;
    int source=node_mat(nd_pt(sample,0),SiO2),xv=imptosol[XVEL],yv=imptosol[YVEL];
    if(source<0 || xv<0 || yv<0)return vc_reject(__LINE__);
    owner[0]=face->t;owner[1]=neigh_fc(face->t,face->f);
    if(owner[1]<0 || dead_tri(owner[1]) || !leaf(tri[owner[0]]) || !leaf(tri[owner[1]]))return vc_reject(__LINE__);
    a=face->a;b=face->b;
    if(xcord(a)==xcord(b))return vc_reject(__LINE__);
    ratio=((double)xcord(sample)-xcord(a))/((double)xcord(b)-xcord(a));
    if(!(ratio>0 && ratio<1))return vc_reject(__LINE__);
    c[0]=xcord(sample);c[1]=(double)ycord(a)+ratio*((double)ycord(b)-ycord(a));
    velocity[0]=pt[sample]->vel[0];velocity[1]=pt[sample]->vel[1];
    start[0]=sol_nd(source,MAXIMP-5);start[1]=sol_nd(source,MAXIMP-6);
    for(side=0;side<2;side++) {
        t=owner[side];region[side]=reg_tri(t);edge[side]=-1;
        for(i=0;i<3;i++) {
            int p=pt_nd(vert_tri(t,(i+1)%3)),q=pt_nd(vert_tri(t,(i+2)%3));
            if((p==a && q==b)||(p==b && q==a)) {
                nodes[side][0]=vert_tri(t,(i+1)%3);nodes[side][1]=vert_tri(t,(i+2)%3);
                nodes[side][2]=vert_tri(t,i);edge[side]=edg_ele(t,i);break;
            }
        }
        if(edge[side]<0)return vc_reject(__LINE__);
        for(i=0;i<3;i++)memcpy(v[i],cord_arr(pt_nd(nodes[side][i])),sizeof(v[i]));
        if(orient2d_float_sign(v[0],c,v[2])<=0 || orient2d_float_sign(c,v[1],v[2])<=0)return vc_reject(__LINE__);
        area=vc_area(v[0],v[1],v[2]);wa=vc_area(v[0],c,v[2])/3;wb=vc_area(c,v[1],v[2])/3;
        for(i=0;i<3;i++) {
            int node=nodes[side][i];stable[side][i]=nd[node];
            for(j=0;j<num_tri_nd(node);j++) {
                int tt=tri_nd(node,j);
                /* A shared material node can span deposition-region labels
                 * (the native-oxide cap closes a hole around the poly gate).
                 * Its field is single-valued across that bookkeeping seam.
                 * Include its entire material support in the inventory. */
                if(dead_tri(tt) || !leaf(tri[tt]) ||
                   mat_reg(reg_tri(tt))!=mat_reg(region[side]))return vc_reject(__LINE__);
                oldw[side][i]+=vc_area(cord_arr(pt_nd(vert_tri(tt,0))),
                    cord_arr(pt_nd(vert_tri(tt,1))),cord_arr(pt_nd(vert_tri(tt,2))))/3;
            }
            if(!(oldw[side][i]>0))return vc_reject(__LINE__);
            neww[side][i]=oldw[side][i]-area/3;
            for(j=0;j<n_imp;j++) {
                fields[side][i][j]=sol_nd(node,j);
                if(!isfinite(fields[side][i][j]))return vc_reject(__LINE__);
            }
        }
        neww[side][0]+=wa;neww[side][1]+=wb;neww[side][2]+=wa+wb;neww[side][3]=wa+wb;
        for(j=0;j<n_imp;j++) {
            double f=pt_nd(nodes[side][0])==a ? ratio:1-ratio;
            fields[side][3][j]=(1-f)*fields[side][0][j]+f*fields[side][1][j];
        }
        for(j=0;j<n_imp;j++) {
            int parent=vc_parent(soltoimp[j]),stream;
            if(parent<0)continue;
            stream=imptosol[parent];if(stream<0)return vc_reject(__LINE__);
            for(i=0;i<4;i++) {
                double total=fields[side][i][stream],active=fields[side][i][j];
                if(total<0 || active<0 || active>total)return vc_reject(__LINE__);
                fractions[side][i][j]=total>0?active/total:0;
            }
        }
        for(j=0;j<n_imp;j++) {
            if(vc_density(soltoimp[j],mat_reg(region[side]))) {
                for(i=0;i<4;i++)values[i]=fields[side][i][j];
                if(!vc_transfer(4,oldw[side],neww[side],values)) {
                    if(getenv("SUPREM4_CONSTRAINT_TRACE")) {
                        static int reports=0;
                        if(reports++<8)fprintf(stderr,"[vertical-constraints] transfer infeasible imp=%d material=%d values=%.17g,%.17g,%.17g,%.17g\n",
                            soltoimp[j],mat_reg(region[side]),fields[side][0][j],fields[side][1][j],fields[side][2][j],fields[side][3][j]);
                    }
                    return vc_reject(__LINE__);
                }
                for(i=0;i<4;i++)fields[side][i][j]=values[i];
            }
        }
        for(j=0;j<n_imp;j++) {
            int parent=vc_parent(soltoimp[j]),stream;
            if(parent<0)continue;
            stream=imptosol[parent];if(stream<0)return vc_reject(__LINE__);
            for(i=0;i<4;i++)fields[side][i][j]=fractions[side][i][j]*fields[side][i][stream];
        }
    }
    /* Existing shared-edge primitive constructs the physical point once and
     * subdivides both owners. Geometry and all new fields were checked above. */
    {
    static int traced=0;
    int trace=getenv("SUPREM4_CONSTRAINT_TRACE") && !traced;
    if(trace)ig2_write("constraint_before.str",0,1.0);
    sp_gredg(edge[0],c,&newpoint);newpt=pt[newpoint];
    for(side=0;side<2;side++) {
        newnode=node_mat(nd_pt(newpoint,0),mat_reg(region[side]));
        if(newnode<0)panic("vertical constraint split lost a material owner");
        for(i=0;i<3;i++)for(j=0;j<n_imp;j++)stable[side][i]->sol[j]=fields[side][i][j];
        for(j=0;j<n_imp;j++)nd[newnode]->sol[j]=fields[side][3][j];
        nd[newnode]->sol[xv]=velocity[0];nd[newnode]->sol[yv]=velocity[1];
        nd[newnode]->sol[MAXIMP-5]=start[0];nd[newnode]->sol[MAXIMP-6]=start[1];
    }
    newpt->vel[0]=velocity[0];newpt->vel[1]=velocity[1];
    if(trace)ig2_write("constraint_split.str",0,1.0);
    bd_connect("synchronized vertical constraints");UpdateSymbolic++;
    if(trace){ig2_write("constraint_after.str",0,1.0);traced=1;}
    }
    return 1;
}

int vertical_constraint_refine(double dt)
{
    struct vc_face *faces;int count=0,t,f,n,i,j,k,point,changed=0,iteration;
    const char *option=getenv("SUPREM4_VERTICAL_CONSTRAINTS");
    if(!option || strcmp(option,"1") || mode!=TWOD)return 0;
    {
    static int entry_reported=0;
    if(getenv("SUPREM4_CONSTRAINT_TRACE") && !entry_reported) {
        fprintf(stderr,"[vertical-constraints] enabled dt=%.17g\n",dt);entry_reported=1;
    }
    }
    for(iteration=0;iteration<64;iteration++) {
        faces=salloc(struct vc_face,3*ne);count=0;
        for(t=0;t<ne;t++) {
            if(dead_tri(t) || !leaf(tri[t]) || mat_reg(reg_tri(t))!=SiO2)continue;
            for(f=0;f<3;f++) {
                struct vc_face *p;
                n=neigh_fc(t,f);
                if(n!=EXPOSED && (n<0 || mat_reg(reg_tri(n))!=Poly))continue;
                p=&faces[count++];p->t=t;p->f=f;p->type=n==EXPOSED?0:1;
                p->a=pt_nd(nd_face_ele(t,0,f));p->b=pt_nd(nd_face_ele(t,1,f));
                for(k=0;k<2;k++) {
                    point=k?p->b:p->a;
                    p->xy[k][0]=(double)xcord(point)+pt[point]->vel[0]*dt;
                    p->xy[k][1]=(double)ycord(point)+pt[point]->vel[1]*dt;
                }
                p->lo=fmin(p->xy[0][0],p->xy[1][0]);p->hi=fmax(p->xy[0][0],p->xy[1][0]);
            }
        }
        qsort(faces,count,sizeof(*faces),face_order);n=0;
        for(i=0;i<count && !n;i++)for(j=i+1;j<count && faces[j].lo<=faces[i].hi && !n;j++) {
            struct vc_face *a=&faces[i],*b=&faces[j],*coarse,*fine;
            if(a->type==b->type || a->a==b->a || a->a==b->b || a->b==b->a || a->b==b->b)continue;
            if(fmax(a->xy[0][1],a->xy[1][1])<fmin(b->xy[0][1],b->xy[1][1]) ||
               fmax(b->xy[0][1],b->xy[1][1])<fmin(a->xy[0][1],a->xy[1][1]))continue;
            if(orient2d_float_sign(a->xy[0],a->xy[1],b->xy[0])*
               orient2d_float_sign(a->xy[0],a->xy[1],b->xy[1])>0 ||
               orient2d_float_sign(b->xy[0],b->xy[1],a->xy[0])*
               orient2d_float_sign(b->xy[0],b->xy[1],a->xy[1])>0)continue;
            coarse=a->type?a:b;fine=a->type?b:a;
            if(pt[coarse->a]->vel[0]!=0 || pt[coarse->b]->vel[0]!=0 ||
               pt[fine->a]->vel[0]!=0 || pt[fine->b]->vel[0]!=0)continue;
            for(k=0;k<2 && !n;k++)n=vc_split(coarse,k?fine->b:fine->a);
        }
        free(faces);if(!n)break;changed+=n;
    }
    if(iteration==64)panic("vertical constraint refinement work limit exceeded");
    if(changed)fprintf(stderr,"[vertical-constraints] inserted=%d dt=%.17g\n",changed,dt);
    return changed;
}
