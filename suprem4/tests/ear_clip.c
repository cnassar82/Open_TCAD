/* Thin contact-etch polygon captured from the full +1% spacing run. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "dbaccess.h"
#include "skel.h"

extern struct LLedge *chop();
void panic(char *message) { fprintf(stderr, "%s\n", message); abort(); }

static double side(int a, int b, int c)
{
    double ax=pt[a]->cord[0], ay=pt[a]->cord[1];
    double bx=pt[b]->cord[0], by=pt[b]->cord[1];
    double cx=pt[c]->cord[0], cy=pt[c]->cord[1];
    return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
}

static void check(const float coordinates[5][2], int valid, const double *angles, int must)
{
    struct pt_str points[5] = {0}, *point_ptrs[5];
    struct nd_str nodes[5] = {0}, *node_ptrs[5];
    struct edg_str edges[5] = {0}, *edge_ptrs[5];
    struct LLedge links[5] = {0}, *ear;
    struct sreg_str region = {0};
    int i, j, rotation, a,b,c, order[5], tmp;
    double total, remainder;
    pt=point_ptrs; nd=node_ptrs; edg=edge_ptrs; mode=TWOD;
    region.len=5;
    for(i=0;i<5;i++) {
        pt[i]=&points[i]; nd[i]=&nodes[i]; edg[i]=&edges[i];
        memcpy(points[i].cord,coordinates[i],sizeof(coordinates[i]));
        points[i].spac=6e-6;
        nodes[i].pt=i; edges[i].nd[0]=i; edges[i].nd[1]=(i+1)%5;
        links[i].edge=i; links[i].iscc=TRUE;
        links[i].prev=&links[(i+4)%5]; links[i].next=&links[(i+1)%5];
        /* Exercise the fallback when cached angles obscure valid ears. */
        links[i].ang=angles ? angles[i] : PI; order[i]=i;
    }
    for(i=1;i<5;i++) for(j=i;j>0 && links[order[j]].ang<links[order[j-1]].ang;j--) {
        tmp=order[j]; order[j]=order[j-1]; order[j-1]=tmp;
    }
    for(i=0;i<5;i++) {
        links[order[i]].gt=&links[order[(i+1)%5]];
        links[order[i]].lt=&links[order[(i+4)%5]];
    }
    region.maxa=&links[order[4]];
    for(rotation=0;rotation<5;rotation++) {
        region.bnd=&links[rotation];
        mgeom=0;
        ear=chop(&region,must);
        if(!valid) { assert(ear==NULL); continue; }
        if (!ear) fprintf(stderr, "No ear: start=(%.9g,%.9g) rotation=%d cached=%d\n",coordinates[0][0],coordinates[0][1],rotation,angles!=NULL);
        assert(ear!=NULL);
        a=nB(ear->prev); b=nB(ear); c=nB(ear->next);
        assert(side(a,b,c)>0);
        for(i=0;i<5;i++) if(i!=a && i!=b && i!=c)
            assert(side(a,b,i)<0 || side(b,c,i)<0 || side(c,a,i)<0);
        total=side(0,1,2)+side(0,2,3)+side(0,3,4);
        remainder=side(a,c,(c+1)%5)+side(a,(c+1)%5,(c+2)%5);
        assert(remainder>0);
        assert(fabs(total-side(a,b,c)-remainder)<fabs(total)*1e-10);
    }
}

int main(void)
{
    const float observed[5][2] = {
        {.012916618026793003,-1.5717132555437274e-5},
        {.012906094081699848,9.318651791545562e-6},
        {.012908217497169971,4.26568249167758e-6},
        {.012910649180412292,-1.5182316701611853e-6},
        {.012913431972265244,-8.138826160575263e-6}
    };
    const float concave[5][2]={{0,0},{3,0},{1,1},{3,3},{0,3}};
    const float crossed[5][2]={{0,0},{3,3},{0,3},{3,0},{1,-1}};
    const float flat[5][2]={{0,0},{1,0},{2,0},{3,0},{4,0}};
    const float duplicate[5][2]={{0,0},{2,0},{2,0},{2,2},{0,2}};
    const double cached_angles[5]={5.2377244576314215e-5,-7.215436538540211e-6,PI,PI,PI};
    check(observed,1,cached_angles,TRUE); check(observed,1,NULL,TRUE); check(concave,1,NULL,TRUE);
    check(crossed,0,NULL,TRUE); check(flat,0,NULL,TRUE); check(duplicate,0,NULL,TRUE);
    {
        const double misleading[5]={1,1,1,1,1};
        check(observed,1,misleading,FALSE);
    }
    puts("PASS: thin polygon ears preserve area and reject invalid boundaries");
    return 0;
}
