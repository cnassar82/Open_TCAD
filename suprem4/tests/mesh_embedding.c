#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mesh_embedding.h"
int main(int argc,char **argv)
{
    /* Both captured LOCOS triangles have positive exact orientation. */
    float captured[5][3]={{.5001741374144331e-4,.14980503692640923e-4},
        {.4999007433070801e-4,.14964174624765292e-4},
        {.5001742465537973e-4,.14980397281760816e-4},
        {.5001794852432795e-4,.14956844097469002e-4},
        {.5004935883334838e-4,.14984182598709594e-4}};
    struct embedding_triangle hit[2]={{{0,1,2},7,14679},{{3,4,0},7,14701}};
    const int side=argc>1?atoi(argv[1]):150;
    int n=(side+1)*(side+1),nt=2*side*side,bad,i,j,k,t=0,v[4],a,b;
    float (*xy)[3]=calloc(n,sizeof(*xy));
    struct embedding_triangle *tri=calloc(nt,sizeof(*tri));
    clock_t start;
    assert(side>0 && side<=2000);
    assert(xy && tri);
    for(i=0;i<2;i++)assert(orient2d_float_sign(captured[hit[i].vertex[0]],captured[hit[i].vertex[1]],captured[hit[i].vertex[2]])>0);
    assert(!mesh_embedding_check((const float (*)[3])captured,2,hit,&bad));
    assert(bad==14679 || bad==14701);
    for(j=0;j<=side;j++)for(i=0;i<=side;i++){xy[j*(side+1)+i][0]=i;xy[j*(side+1)+i][1]=j;}
    for(j=0;j<side;j++)for(i=0;i<side;i++) {
        v[0]=j*(side+1)+i;v[1]=v[0]+1;v[2]=v[1]+side+1;v[3]=v[0]+side+1;
        tri[t].vertex[0]=v[0];tri[t].vertex[1]=v[1];tri[t].vertex[2]=v[2];t++;
        tri[t].vertex[0]=v[0];tri[t].vertex[1]=v[2];tri[t].vertex[2]=v[3];t++;
    }
    for(t=0;t<nt;t++) {
        tri[t].identity=t;
        for(k=0;k<3;k++) {
            a=tri[t].vertex[k];b=tri[t].vertex[(k+1)%3];
            if((xy[a][0]==0 && xy[b][0]==0) || (xy[a][0]==side && xy[b][0]==side) ||
               (xy[a][1]==0 && xy[b][1]==0) || (xy[a][1]==side && xy[b][1]==side))tri[t].boundary|=1<<k;
        }
    }
    start=clock();assert(mesh_embedding_check((const float (*)[3])xy,nt,tri,&bad));
    printf("PASS: captured positive-area overlap rejected; %d-triangle conforming mesh validated in %.3f CPU seconds\n",nt,(double)(clock()-start)/CLOCKS_PER_SEC);
    free(xy);free(tri);return 0;
}
