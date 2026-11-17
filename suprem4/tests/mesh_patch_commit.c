#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "impurity.h"
#include "mesh_patch.h"

#ifdef SUPREM_MESH_PATCH_TEST
#undef malloc
#undef calloc
extern void *malloc(size_t);
extern void *calloc(size_t,size_t);
static int allocation_index, fail_after=-1;
void *mesh_patch_test_malloc(size_t n)
{
    if(fail_after>=0 && allocation_index++==fail_after)return NULL;
    return malloc(n);
}
void *mesh_patch_test_calloc(size_t n,size_t size)
{
    if(fail_after>=0 && allocation_index++==fail_after)return NULL;
    return calloc(n,size);
}
#endif

static unsigned long long hash_bytes(unsigned long long h,const void *p,size_t n)
{
    const unsigned char *s=p;size_t i;
    for(i=0;i<n;i++){h^=s[i];h*=1099511628211ULL;}return h;
}
#define HASH(h,x) h=hash_bytes(h,&(x),sizeof(x))
static unsigned long long semantic_hash(void)
{
    unsigned long long h=1469598103934665603ULL;int i,j;
    HASH(h,np);HASH(h,nn);HASH(h,ne);HASH(h,nreg);HASH(h,need_waste);
    HASH(h,geom_dirty);HASH(h,clkws_dirty);HASH(h,neigh_dirty);HASH(h,edge_dirty);
    HASH(h,nd2edg_dirty);HASH(h,reg_dirty);
    for(i=0;i<np;i++) {
        HASH(h,pt[i]->flags);HASH(h,pt[i]->nn);
        for(j=0;j<MAXDIM;j++){HASH(h,pt[i]->cord[j]);HASH(h,pt[i]->vel[j]);}
        for(j=0;j<pt[i]->nn;j++)HASH(h,pt[i]->nd[j]);
    }
    for(i=0;i<nn;i++) {
        HASH(h,nd[i]->flags);HASH(h,nd[i]->pt);HASH(h,nd[i]->mater);
        for(j=0;j<n_imp;j++)HASH(h,nd[i]->sol[j]);
        HASH(h,nd[i]->ltri.num);HASH(h,nd[i]->ltri.all);
        for(j=0;j<nd[i]->ltri.num;j++)HASH(h,nd[i]->ltri.list[j]);
    }
    for(i=0;i<ne;i++) {
        HASH(h,tri[i]->flags);HASH(h,tri[i]->regnum);
        for(j=0;j<3;j++){HASH(h,tri[i]->nd[j]);HASH(h,tri[i]->nb[j]);}
    }
    HASH(h,reg[0]->tri.num);HASH(h,reg[0]->tri.all);
    for(j=0;j<reg[0]->tri.num;j++)HASH(h,reg[0]->tri.list[j]);
    return h;
}

static void list_add(struct list_str *l,int t)
{
    l->list=realloc(l->list,(l->num+1)*sizeof(int));l->list[l->num++]=t;l->all=l->num;
}
static void fixture(void)
{
    int i,j;float xy[5][2]={{0,0},{-1,-1},{1,-1},{1,1},{-1,1}};
    mode=TWOD;nvrt=3;nedg=3;np=5;nn=9;ne=4;nreg=1;ned=0;n_imp=5;
    soltoimp[0]=B;soltoimp[1]=Ba;soltoimp[2]=P;soltoimp[3]=XVEL;soltoimp[4]=I;
    for(i=0;i<MAXIMP;i++)imptosol[i]=-1;
    imptosol[B]=0;imptosol[Ba]=1;imptosol[P]=2;imptosol[XVEL]=3;imptosol[I]=4;
    pt=calloc(np,sizeof(*pt));nd=calloc(nn,sizeof(*nd));tri=calloc(ne,sizeof(*tri));
    reg[0]=calloc(1,sizeof(*reg[0]));reg[0]->mater=SiO2;
    for(i=0;i<np;i++) {
        pt[i]=calloc(1,sizeof(**pt));pt[i]->cord[0]=xy[i][0];pt[i]->cord[1]=xy[i][1];
        pt[i]->nn=1;pt[i]->nd[0]=i;
        nd[i]=calloc(1,sizeof(**nd));nd[i]->pt=i;nd[i]->mater=SiO2;
        nd[i]->sol[0]=1e18*(i+1);nd[i]->sol[1]=.3*nd[i]->sol[0];
        nd[i]->sol[2]=4e17;nd[i]->sol[3]=i-2;nd[i]->sol[4]=1;
        if(i>0) {
            int gas=i+4;
            set_surf(i);pt[i]->nd[pt[i]->nn++]=gas;
            nd[gas]=calloc(1,sizeof(**nd));nd[gas]->pt=i;nd[gas]->mater=GAS;
        }
    }
    for(i=0;i<ne;i++) {
        int vv[3]={0,i+1,(i+1)%4+1};
        tri[i]=calloc(1,sizeof(**tri));tri[i]->fath=tri[i]->son=NO_TRIANGLE;
        tri[i]->regnum=0;
        for(j=0;j<3;j++){tri[i]->nd[j]=vv[j];list_add(&nd[vv[j]]->ltri,i);}
        tri[i]->nb[0]=EXPOSED;tri[i]->nb[1]=(i+1)%4;tri[i]->nb[2]=(i+3)%4;
        list_add(&reg[0]->tri,i);
    }
}
static void add_external_triangle(void)
{
    int j;int vv[3]={2,1,9};
    pt=realloc(pt,6*sizeof(*pt));nd=realloc(nd,11*sizeof(*nd));tri=realloc(tri,5*sizeof(*tri));
    pt[5]=calloc(1,sizeof(**pt));pt[5]->cord[1]=-2;pt[5]->nn=2;pt[5]->nd[0]=9;pt[5]->nd[1]=10;
    set_surf(5);
    nd[9]=calloc(1,sizeof(**nd));nd[9]->pt=5;nd[9]->mater=SiO2;
    nd[9]->sol[0]=2e18;nd[9]->sol[1]=.3*nd[9]->sol[0];nd[9]->sol[2]=4e17;
    nd[9]->sol[3]=3;nd[9]->sol[4]=1;
    nd[10]=calloc(1,sizeof(**nd));nd[10]->pt=5;nd[10]->mater=GAS;
    tri[4]=calloc(1,sizeof(**tri));tri[4]->fath=tri[4]->son=NO_TRIANGLE;
    for(j=0;j<3;j++){tri[4]->nd[j]=vv[j];list_add(&nd[vv[j]]->ltri,4);tri[4]->nb[j]=EXPOSED;}
    tri[4]->nb[2]=0;tri[0]->nb[0]=4;list_add(&reg[0]->tri,4);
    np=6;nn=11;ne=5;
}
static void boundary_fixture(void)
{
    int i,j;
    float xy[5][2]={{-.99f,0},{-1,0},{-1,-1},{1,-1},{1,0}};
    fixture();
    for(i=0;i<np;i++) {
        pt[i]->cord[0]=xy[i][0];pt[i]->cord[1]=xy[i][1];
        free(nd[i]->ltri.list);memset(&nd[i]->ltri,0,sizeof(struct list_str));
    }
    free(tri[3]);ne=3;reg[0]->tri.num=3;
    for(i=0;i<ne;i++) {
        tri[i]->nd[0]=0;tri[i]->nd[1]=i+1;tri[i]->nd[2]=i+2;
        tri[i]->nb[0]=EXPOSED;
        tri[i]->nb[1]=i==2?BC_OFFSET:i+1;
        tri[i]->nb[2]=i==0?BC_OFFSET:i-1;
        for(j=0;j<3;j++)list_add(&nd[tri[i]->nd[j]]->ltri,i);
    }
}

static double inventory(int stream)
{
    int t,j;double sum=0;
    for(t=0;t<ne;t++)if(!dead_tri(t)) {
        int a=pt_nd(tri[t]->nd[0]),b=pt_nd(tri[t]->nd[1]),c=pt_nd(tri[t]->nd[2]);
        double area=((xcord(b)-xcord(a))*(ycord(c)-ycord(a))-
                     (ycord(b)-ycord(a))*(xcord(c)-xcord(a)))/2.;
        assert(area>0);
        for(j=0;j<3;j++)sum+=area*sol_nd(tri[t]->nd[j],stream)/3.;
    }
    return sum;
}
static void check(void)
{
    int t,j,k,i,gas;
    for(t=0;t<ne;t++)if(!dead_tri(t))for(j=0;j<3;j++) {
        int n=tri[t]->nd[j],found=0;
        assert(!dead_nd(n));
        for(i=0;i<nd[n]->ltri.num;i++)if(nd[n]->ltri.list[i]==t)found++;
        assert(found==1);
        k=tri[t]->nb[j];
        if(k>=0) {
            assert(!dead_tri(k));found=0;
            for(i=0;i<3;i++)if(tri[k]->nb[i]==t)found++;
            assert(found==1);
        }
    }
    for(i=1;i<np;i++) {
        gas=0;for(j=0;j<num_nd(i);j++)if(mat_nd(nd_pt(i,j))==GAS)gas++;
        assert(is_surf(i) && gas==1);
        if(i>=5)continue;
        assert(sol_nd(i,0)>=1e18 && sol_nd(i,0)<=5e18);
        assert(fabs(sol_nd(i,1)/sol_nd(i,0)-.3)<1e-15);
        assert(sol_nd(i,2)==4e17 && sol_nd(i,3)==i-2 && sol_nd(i,4)==1);
    }
}
static void cleanup(void)
{
    int i;for(i=0;i<np;i++)free(pt[i]);
    for(i=0;i<nn;i++){free(nd[i]->ltri.list);free(nd[i]);}
    for(i=0;i<ne;i++)free(tri[i]);
    free(reg[0]->tri.list);free(reg[0]);free(pt);free(nd);free(tri);
}
int main(void)
{
    const char *reason;double before[3];int status,i,variant;
    fixture();
    before[0]=inventory(0);before[1]=inventory(2);before[2]=inventory(4);
    status=mesh_patch_remove_point(0,&reason);
    if(status!=MESH_PATCH_OK)fprintf(stderr,"proposal failed: %d %s\n",status,reason);
    assert(status==MESH_PATCH_OK);assert(dead_pt(0) && dead_nd(0));
    for(i=0;i<3;i++) {
        int stream=i==0?0:i==1?2:4;
        assert(fabs(inventory(stream)-before[i])<=1e-12*before[i]);
    }
    check();cleanup();
    fixture();
    { struct pt_str point=*pt[0];struct nd_str node=*nd[0];struct tri_str triangle=*tri[0];
      set_surf(0);point.flags=pt[0]->flags;
      assert(mesh_patch_remove_point(0,&reason)==MESH_PATCH_PINNED);
      assert(memcmp(&point,pt[0],sizeof point)==0);
      assert(memcmp(&node,nd[0],sizeof node)==0);
      assert(memcmp(&triangle,tri[0],sizeof triangle)==0);
    }
    cleanup();
    fixture();add_external_triangle();
    before[0]=inventory(0);before[1]=inventory(2);before[2]=inventory(4);
    assert(mesh_patch_remove_point(0,&reason)==MESH_PATCH_OK);
    for(i=0;i<3;i++) {
        int stream=i==0?0:i==1?2:4;
        assert(fabs(inventory(stream)-before[i])<=1e-12*before[i]);
    }
    check();assert(!dead_tri(4));cleanup();
    fixture();add_external_triangle();tri[4]->regnum=1;
    { unsigned long long before_state=semantic_hash();
      assert(mesh_patch_remove_point(0,&reason)==MESH_PATCH_PINNED);
      assert(strcmp(reason,"field support spans region boundary")==0);
      assert(before_state==semantic_hash());
    }
    cleanup();
    boundary_fixture();
    before[0]=inventory(0);before[1]=inventory(2);before[2]=inventory(4);
    assert(mesh_patch_remove_point(0,&reason)==MESH_PATCH_OK);
    assert(dead_pt(0) && dead_nd(0));
    {int boundary_edges=0,t,j;
     for(t=0;t<ne;t++)if(!dead_tri(t))for(j=0;j<3;j++)
         if(tri[t]->nb[j]==BC_OFFSET)boundary_edges++;
     assert(boundary_edges==1);}
    for(i=0;i<3;i++) {
        int stream=i==0?0:i==1?2:4;
        assert(fabs(inventory(stream)-before[i])<=1e-12*before[i]);
    }
    check();cleanup();
    boundary_fixture();tri[0]->nb[2]=BC_OFFSET+1;
    {unsigned long long saved=semantic_hash();
     assert(mesh_patch_remove_point(0,&reason)==MESH_PATCH_PINNED);
     assert(saved==semantic_hash());}
    cleanup();
#ifdef SUPREM_MESH_PATCH_TEST
    for(variant=0;variant<2;variant++) {
    for(i=0;i<100;i++) {
        unsigned long long before_state;
        if(variant)boundary_fixture();else fixture();
        before_state=semantic_hash();allocation_index=0;fail_after=i;
        status=mesh_patch_remove_point(0,&reason);fail_after=-1;
        if(status==MESH_PATCH_OK){cleanup();break;}
        assert(status==MESH_PATCH_NOMEM);
        assert(before_state==semantic_hash());cleanup();
    }
    assert(i>15 && i<100);
    printf("PASS injected failure at each of %d preparation allocations leaves semantic state unchanged\n",i);
    }
#endif
    puts("PASS interior/straight-boundary atomic commit, inventories, boundary codes, exposure and refusal");
    return 0;
}
