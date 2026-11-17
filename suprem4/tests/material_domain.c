#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "material_domain.h"
#ifdef SUPREM_DOMAIN_TEST
#undef malloc
#undef calloc
extern void *malloc(size_t);
extern void *calloc(size_t,size_t);
static int alloc_at,fail_at=-1;
void *domain_test_malloc(size_t n){if(fail_at>=0 && alloc_at++==fail_at)return NULL;return malloc(n);}
void *domain_test_calloc(size_t n,size_t s){if(fail_at>=0 && alloc_at++==fail_at)return NULL;return calloc(n,s);}
#endif
int main(void)
{
    struct domain_region r[4]={{1,1,0},{1,1,0},{1,1,0},{1,1,0}};
    struct domain_triangle square[2]={{{0,1,2},0},{{0,2,3},1}};
    struct domain_triangle ring[8],bad[3],touch[2]={{{0,1,2},0},{{0,3,4},1}};
    int map[4],n,i,status;
    assert(material_domain_plan(4,2,r,2,square,map,&n)==0 && n==1 && map[0]==map[1]);
    r[1].material=3;
    assert(material_domain_plan(4,2,r,2,square,map,&n)==0 && n==2);
    r[1].material=1;r[1].model_key=2;
    assert(material_domain_plan(4,2,r,2,square,map,&n)==0 && n==2);
    r[1].model_key=1;r[1].protected_region=1;
    assert(material_domain_plan(4,2,r,2,square,map,&n)==0 && n==2);
    r[1].protected_region=0;
    assert(material_domain_plan(5,2,r,2,touch,map,&n)==0 && n==2);
    for(i=0;i<4;i++) {
        int j=(i+1)%4;
        ring[2*i].vertex[0]=i;ring[2*i].vertex[1]=j;ring[2*i].vertex[2]=4+j;
        ring[2*i+1].vertex[0]=i;ring[2*i+1].vertex[1]=4+j;ring[2*i+1].vertex[2]=4+i;
        ring[2*i].region=ring[2*i+1].region=i;
    }
    assert(material_domain_plan(8,4,r,8,ring,map,&n)==0 && n==2);
    /* One semantic region can have holes in a modern domain model. This
     * adapter must preserve cuts until the legacy skeleton supports them. */
    bad[0]=square[0];bad[1]=square[1];bad[2]=square[0];
    map[0]=37;map[1]=38;n=39;
    assert(material_domain_plan(4,2,r,3,bad,map,&n)==-1);
    assert(map[0]==37 && map[1]==38 && n==39);
    for(i=0;i<6;i++) {
        fail_at=i;alloc_at=0;
        status=material_domain_plan(4,2,r,2,square,map,&n);
        assert(status==-2 && map[0]==37 && map[1]==38 && n==39);
    }
    fail_at=-1;
    assert(material_domain_plan(4,2,r,2,square,map,&n)==0 && n==1);
    puts("PASS same-domain quotient, protected physics, point contact, holes, invalid incidence and all 6 allocation failures");
    return 0;
}
