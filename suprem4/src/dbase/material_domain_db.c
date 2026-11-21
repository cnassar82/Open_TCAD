/* Legacy bridge for a physical material-domain policy. Run only between
 * process commands, never while a region Boolean is iterating its operands.
 * The current opt-in policy joins oxide with a material-global model. Other
 * materials and seams represented by distinct material nodes stay protected.
 */
#include <stdio.h>
#include "suprem_options.h"
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "material_domain.h"
int material_domain_normalize(void)
{
    const char *policy=suprem_option("SUPREM4_MATERIAL_DOMAINS");
    struct domain_region *properties=NULL;
    struct domain_triangle *faces=NULL;
    struct reg_str **replacement=NULL;
    int *map=NULL,i,j,k,r,n,other,a,b,ndomains,status=-2;
    if(!policy || strcmp(policy,"oxide") || mode!=TWOD)return 0;
    if(need_waste || tripts || neigh_dirty)return 0;
    for(i=0;i<ne;i++)if(dead_tri(i) || !leaf(tri[i]) || !root(tri[i]))return 0;
    properties=calloc(nreg,sizeof(*properties));faces=malloc((size_t)ne*sizeof(*faces));
    map=malloc((size_t)nreg*sizeof(int));replacement=calloc(nreg,sizeof(*replacement));
    if(!properties || !faces || !map || !replacement)goto out;
    for(r=0;r<nreg;r++) {
        properties[r].material=mat_reg(r);properties[r].model_key=mat_reg(r);
        properties[r].protected_region=mat_reg(r)!=SiO2;
    }
    for(i=0;i<ne;i++) {
        r=reg_tri(i);faces[i].region=r;
        for(j=0;j<3;j++) {
            faces[i].vertex[j]=pt_nd(vert_tri(i,j));
            other=neigh_fc(i,j);
            if(other<0 || reg_tri(other)==r || mat_reg(reg_tri(other))!=mat_reg(r))continue;
            /* Distinct fields at a same-material seam are a physical barrier.
             * A shared point alone does not authorize erasing that barrier. */
            for(k=1;k<=2;k++) {
                a=vert_tri(i,(j+k)%3);b=-1;
                for(n=0;n<3;n++)if(pt_nd(vert_tri(other,n))==pt_nd(a))b=vert_tri(other,n);
                if(a!=b)properties[r].protected_region=properties[reg_tri(other)].protected_region=1;
            }
        }
    }
    status=material_domain_plan(np,nreg,properties,ne,faces,map,&ndomains);
    if(status || ndomains==nreg)goto out;
    /* Allocate the complete replacement region table before changing labels.
     * Field arrays, point ownership and triangle geometry remain untouched. */
    for(r=0;r<nreg;r++)if(!replacement[map[r]]) {
        replacement[map[r]]=calloc(1,sizeof(struct reg_str));
        if(!replacement[map[r]]){status=-2;goto out;}
        /* Retain non-derived region metadata from the representative.
         * The historical thick member is only copied in this codebase;
         * it is not a regional physics parameter. */
        memcpy(replacement[map[r]],reg[r],sizeof(struct reg_str));
        memset(&replacement[map[r]]->edg,0,sizeof(struct list_str));
        memset(&replacement[map[r]]->tri,0,sizeof(struct list_str));
        memset(&replacement[map[r]]->fac,0,sizeof(struct list_str));
        memset(&replacement[map[r]]->nd,0,sizeof(struct list_str));
        replacement[map[r]]->sreg=-1;
    }
    fprintf(stderr,"[material-domain] policy=oxide regions=%d domains=%d mapping=",nreg,ndomains);
    for(r=0;r<nreg;r++)fprintf(stderr,"%s%d:%d",r?",":"",r,map[r]);
    fprintf(stderr,"\n");
    /* Publish the quotient; all allocations and topology tests have finished. */
    for(i=0;i<ne;i++)tri[i]->regnum=map[tri[i]->regnum];
    while(nsreg>0)free_skel(0);
    for(r=0;r<nreg;r++) {
        reg[r]->sreg=-1;dis_1reg(&reg[r]);
    }
    for(r=0;r<ndomains;r++){reg[r]=replacement[r];replacement[r]=NULL;}
    nreg=ndomains;reg_dirty=geom_dirty=TRUE;
    status=1;
out:
    if(replacement)for(r=0;r<nreg;r++)free(replacement[r]);
    free(properties);free(faces);free(map);free(replacement);
    if(status<0)panic(status==-2?"Out of memory":"Invalid material-domain incidence");
    if(status==1)bd_connect("canonical oxide material domains");
    return status;
}
