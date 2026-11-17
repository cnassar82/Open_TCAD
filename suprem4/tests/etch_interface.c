/* Neighbor arrays may still reference a killed, unsplit parent. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
extern int test_interface_other_owner(int,int,int);
int main(void)
{
    int i,j;
    int v[4][3]={{0,1,2},{3,4,5},{6,7,8},{9,10,11}};
    int points[12]={0,1,2,0,1,3,0,1,3,0,1,2};
    float xy[4][2]={{0,0},{1,0},{0,1},{0,-1}};
    nvrt=3; nedg=3; np=4; nn=12; ne=4;
    pt=calloc(np,sizeof(*pt));nd=calloc(nn,sizeof(*nd));tri=calloc(ne,sizeof(*tri));
    for(i=0;i<np;i++) {
        pt[i]=calloc(1,sizeof(**pt));pt[i]->cord[0]=xy[i][0];pt[i]->cord[1]=xy[i][1];
    }
    for(i=0;i<nn;i++) {
        int p=points[i];
        nd[i]=calloc(1,sizeof(**nd));nd[i]->pt=p;
        pt[p]->nd[pt[p]->nn++]=i;
        nd[i]->ltri.num=1;nd[i]->ltri.list=malloc(sizeof(int));nd[i]->ltri.list[0]=i/3;
    }
    for(i=0;i<ne;i++) {
        tri[i]=calloc(1,sizeof(**tri));tri[i]->son=NO_TRIANGLE;tri[i]->regnum=i;
        for(j=0;j<3;j++)tri[i]->nd[j]=v[i][j];
    }
    /* 0: target; 1: obsolete parent; 2: live neighbor with separate
     * material nodes; 3: same-side replacement (not a neighbor). */
    set(tri[1],KILL_ME);tri[0]->nb[2]=1;
    assert(test_interface_other_owner(0,2,0)==2);
    set(tri[2],KILL_ME);
    assert(test_interface_other_owner(0,2,0)==-1);
    clr(tri[2],KILL_ME);tri[2]->son=1;
    assert(test_interface_other_owner(0,2,0)==-1);
    for(i=0;i<nn;i++){free(nd[i]->ltri.list);free(nd[i]);}
    for(i=0;i<np;i++)free(pt[i]);
    for(i=0;i<ne;i++)free(tri[i]);
    free(pt);free(nd);free(tri);
    puts("PASS live interface ownership excludes dead parents and same-side replacements");
    return 0;
}
