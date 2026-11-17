/* Retain process construction geometry independently of its mesh embedding.
 * The actual surface snapshot is an exact validity guard, never a tolerance
 * for removing material. Nominal points use cm, as does the mesh database.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "refine.h"
#include "process_surface.h"
static double *nominal,*embedding;
static int count,actual_count,material_owner,read_nominal,read_actual;
void process_surface_clear(void)
{
    free(nominal);free(embedding);nominal=embedding=NULL;
    count=actual_count=read_nominal=read_actual=0;
}
static int complete(void)
{return count>=2 && actual_count>=2 && read_nominal==count && read_actual==actual_count;}
int process_surface_get(int material,int n,const struct line *actual,const double **xy)
{
    int i;
    if(!complete() || material!=material_owner || n!=actual_count)return 0;
    for(i=0;i<n;i++)if(actual[i].x!=embedding[2*i] || actual[i].y!=embedding[2*i+1])return 0;
    /* Material ownership is checked on every sampled exposed point. */
    for(i=0;i<n;i++)if(actual[i].p<0 ||
        node_mat(nd_pt(actual[i].p,0),material_owner)<0)return 0;
    *xy=nominal;return count;
}
void process_surface_record(int material,int n,const double *xy)
{
    struct line *actual=NULL;int i,m;
    process_surface_clear();
    if(mode!=TWOD || n<2)return;
    m=find_surf(&actual);
    if(m<2){free(actual);return;}
    nominal=salloc(double,2*n);embedding=salloc(double,2*m);
    memcpy(nominal,xy,2*n*sizeof(double));
    for(i=0;i<m;i++){embedding[2*i]=actual[i].x;embedding[2*i+1]=actual[i].y;}
    free(actual);count=read_nominal=n;actual_count=read_actual=m;material_owner=material;
}
void process_surface_write(FILE *file,double sx,double sy)
{
    struct line *actual=NULL;const double *xy;int i,n;
    if(!complete() || mode!=TWOD)return;
    n=find_surf(&actual);
    if(!process_surface_get(material_owner,n,actual,&xy)){free(actual);return;}
    free(actual);
    /* Construction records are in cm with hexadecimal binary64 encoding.
     * Avoid a cm->micron->cm round trip destroying an analytic cancellation.
     * Preserve explicit structure scale/flip transformations. */
    sx*=1e-4;sy*=1e-4;
    fprintf(file,"G surface_cm 1 %d %d %d\n",material_owner,count,actual_count);
    for(i=0;i<count;i++)fprintf(file,"G nominal %a %a\n",nominal[2*i]*sx,nominal[2*i+1]*sy);
    for(i=0;i<actual_count;i++)fprintf(file,"G embedding %a %a\n",embedding[2*i]*sx,embedding[2*i+1]*sy);
}
int process_surface_read(const char *line,double sx,double sy)
{
    int version,m,n,k;double x,y;
    sx*=1e4;sy*=1e4;
    if(sscanf(line,"G surface_cm %d %d %d %d",&version,&m,&n,&k)==4) {
        process_surface_clear();
        if(version!=1 || m<0 || m>=MAXMAT || n<2 || k<2 || n>10000000 || k>10000000)return 0;
        nominal=salloc(double,2*n);embedding=salloc(double,2*k);
        count=n;actual_count=k;material_owner=m;return 1;
    }
    if(sscanf(line,"G nominal %le %le",&x,&y)==2) {
        if(!nominal || read_nominal>=count || !isfinite(x*sx) || !isfinite(y*sy))return 0;
        nominal[2*read_nominal]=x*sx;nominal[2*read_nominal+1]=y*sy;read_nominal++;return 1;
    }
    if(sscanf(line,"G embedding %le %le",&x,&y)==2) {
        if(!embedding || read_actual>=actual_count || !isfinite(x*sx) || !isfinite(y*sy))return 0;
        /* The embedding represents stored mesh precision after the same unit
         * conversion as c records; construction values keep binary64. */
        embedding[2*read_actual]=(float)(x*sx);embedding[2*read_actual+1]=(float)(y*sy);read_actual++;return 1;
    }
    return 0;
}
