#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "material_domain.h"
#undef malloc
#undef calloc
extern void *malloc(size_t);
extern void *calloc(size_t,size_t);
static int fail_at=-1,allocation_index,rebuilds;
static jmp_buf failed;
void *domain_test_malloc(size_t n){if(fail_at>=0 && allocation_index++==fail_at)return NULL;return malloc(n);}
void *domain_test_calloc(size_t n,size_t size){if(fail_at>=0 && allocation_index++==fail_at)return NULL;return calloc(n,size);}
void panic(char *message){assert(!strcmp(message,"Out of memory"));longjmp(failed,1);}
int free_skel(int r){assert(0);return 0;}
int dis_1reg(struct reg_str **r)
{
    free((*r)->tri.list);free((*r)->edg.list);free((*r)->nd.list);free((*r)->fac.list);
    free(*r);*r=NULL;return 0;
}
int bd_connect(char *when){assert(nreg==1 && tri[0]->regnum==0 && tri[1]->regnum==0);rebuilds++;return 0;}
static struct pt_str points[4],*pp[4];
static struct nd_str nodes[5],*npointers[5];
static struct tri_str triangles[2],*tp[2];
static void setup(void)
{
    int i,j;const int vv[2][3]={{0,1,2},{0,2,3}};
    memset(points,0,sizeof points);memset(nodes,0,sizeof nodes);memset(triangles,0,sizeof triangles);
    pt=pp;nd=npointers;tri=tp;np=4;nn=4;ne=2;nreg=2;nsreg=0;mode=TWOD;rebuilds=0;
    need_waste=tripts=neigh_dirty=geom_dirty=reg_dirty=0;
    for(i=0;i<4;i++){pt[i]=&points[i];pt[i]->cord[0]=(i==1 || i==2);pt[i]->cord[1]=i/2;}
    for(i=0;i<5;i++) {
        nd[i]=&nodes[i];nd[i]->pt=i%4;nd[i]->mater=SiO2;
        for(j=0;j<MAXIMP;j++)nd[i]->sol[j]=(j+1)*1e16+i;
    }
    for(i=0;i<2;i++) {
        tri[i]=&triangles[i];tri[i]->regnum=i;tri[i]->son=tri[i]->fath=NO_TRIANGLE;
        for(j=0;j<3;j++){tri[i]->nd[j]=vv[i][j];tri[i]->nb[j]=-1;}
        reg[i]=calloc(1,sizeof(struct reg_str));reg[i]->mater=SiO2;reg[i]->sreg=-1;reg[i]->thick=.5;
        reg[i]->tri.list=malloc(sizeof(int));reg[i]->tri.list[0]=i;reg[i]->tri.num=reg[i]->tri.all=1;
    }
    tri[0]->nb[1]=1;tri[1]->nb[2]=0;
}
static void cleanup(void){int i;for(i=0;i<nreg;i++)dis_1reg(&reg[i]);}
int main(void)
{
    struct pt_str saved_points[4];struct nd_str saved_nodes[5];
    struct tri_str saved_triangles[2];int attempt,i;
    setenv("SUPREM4_MATERIAL_DOMAINS","oxide",1);
    for(attempt=0;attempt<30;attempt++) {
        fail_at=-1;setup();memcpy(saved_points,points,sizeof points);memcpy(saved_nodes,nodes,sizeof nodes);
        memcpy(saved_triangles,triangles,sizeof triangles);allocation_index=0;fail_at=attempt;
        if(!setjmp(failed)) {
            assert(material_domain_normalize()==1);assert(nreg==1 && rebuilds==1);
            assert(!memcmp(saved_points,points,sizeof points));assert(!memcmp(saved_nodes,nodes,sizeof nodes));
            for(i=0;i<2;i++){saved_triangles[i].regnum=0;assert(!memcmp(&saved_triangles[i],tri[i],sizeof(*tri[i])));}
            cleanup();break;
        }
        assert(nreg==2 && rebuilds==0 && !reg_dirty && !geom_dirty);
        assert(!memcmp(saved_points,points,sizeof points));assert(!memcmp(saved_nodes,nodes,sizeof nodes));
        assert(!memcmp(saved_triangles,triangles,sizeof triangles));
        for(i=0;i<2;i++)assert(reg[i]->tri.num==1 && reg[i]->tri.list[0]==i);
        cleanup();
    }
    assert(attempt==11);fail_at=-1;setup();
    /* Same material but different nodes on a shared point: keep the seam. */
    nn=5;tri[1]->nd[0]=4;nd[4]->sol[0]=1;
    assert(material_domain_normalize()==0 && nreg==2 && rebuilds==0);cleanup();
    puts("PASS domain commit preserves geometry/fields; discontinuous seam protected; all 11 allocation failures leave state unchanged");
    return 0;
}
