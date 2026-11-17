#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
extern int test_removal_polygon_valid(int *,int);
extern int test_removal_ears(int *,int,int *);
void panic(char *message) {fprintf(stderr,"%s\n",message);abort();}
static void check(float xy[][2],int count,int expected)
{
    int nodes[16],i;
    np=nn=count;pt=calloc(count,sizeof(*pt));nd=calloc(count,sizeof(*nd));
    for(i=0;i<count;i++) {
        pt[i]=calloc(1,sizeof(**pt));nd[i]=calloc(1,sizeof(**nd));nodes[i]=i;nd[i]->pt=i;
        pt[i]->cord[0]=xy[i][0];pt[i]->cord[1]=xy[i][1];
    }
    assert(test_removal_polygon_valid(nodes,count)==expected);
    if(expected) {
        int ears[42],nt=test_removal_ears(nodes,count,ears),j;
        double area=0,reference=0,cross;
        assert(nt==count-2);
        for(i=0;i<count;i++)reference+=(double)xy[i][0]*xy[(i+1)%count][1]-(double)xy[i][1]*xy[(i+1)%count][0];
        for(j=0;j<nt;j++) {
            int a=ears[3*j],b=ears[3*j+1],c=ears[3*j+2];
            assert(orient2d_float_sign(pt[a]->cord,pt[b]->cord,pt[c]->cord)>0);
            cross=((double)xy[b][0]-xy[a][0])*(xy[c][1]-xy[a][1])-((double)xy[b][1]-xy[a][1])*(xy[c][0]-xy[a][0]);
            area+=cross;
        }
        assert(fabs(area-reference)<1e-6);
    }
    for(i=0;i<count;i++){free(pt[i]);free(nd[i]);}
    free(pt);free(nd);
}
int main(void)
{
    float concave[5][2]={{0,0},{2,0},{2,2},{1,1},{0,2}};
    float outside[3][2]={{1,0},{.2,.2},{0,1}};
    float zigzag[7][2]={{0,0},{4,0},{4,1},{3,.2},{2,.9},{1,.1},{0,1}};
    float bowtie[4][2]={{0,0},{2,2},{0,2},{2,0}};
    float reverse[5][2]={{0,0},{2,0},{1,0},{2,2},{0,2}};
    float captured[9][2]={
        {.00084136839723214507,2.2692116544931196e-5},
        {.00084338599117472768,2.2823034669272602e-5},
        {.00084344157949090004,2.2822572645964101e-5},
        {.00084349582903087139,2.2813826944911852e-5},
        {.0008435766794718802,2.2800788428867236e-5},
        {.00084369722753763199,2.2781350708100945e-5},
        {.00084387685637921095,2.2752374206902459e-5},
        {.00084414472803473473,2.2709178665536456e-5},
        {.00084672227967530489,2.2939466362004168e-5}};
    check(concave,5,1);check(zigzag,7,1);check(outside,3,0);check(bowtie,4,0);check(reverse,5,0);check(captured,9,0);
    puts("PASS removal cavity simplicity and captured unsafe closing chord");return 0;
}
