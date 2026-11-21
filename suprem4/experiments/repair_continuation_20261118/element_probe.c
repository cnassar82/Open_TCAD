#include <stdio.h>
#include <string.h>
int noderiv=0;
float oxide_obfix=1;
int dcopy(double *a,double *b,int n){memcpy(b,a,n*sizeof(double));return 0;}
int tri7_stiff();int triox_stiff();
int oxide(double *m,int mat,double *e,double *s) {
 memset(m,0,40*sizeof(double));m[0]=3;m[1]=1;m[2]=1;
 m[5]=1;m[10]=3;m[11]=1;m[13]=1;m[14]=3;m[18]=1;m[20]=1;
 s[0]=3*e[0]+e[1];s[1]=e[0]+3*e[1];s[2]=e[2];return 0;
}
int elastic(double *m,int mat,double *e,double *s) {
 memset(m,0,20*sizeof(double));m[0]=3;m[1]=1;m[2]=1;
 s[0]=3*e[0]+e[1];s[1]=e[0]+3*e[1];s[2]=e[2];return 0;
}
int main(void){
 double xl[7][2]={{0,0},{1,0},{0,1},{.5,.5},{0,.5},{.5,0},{1./3,1./3}};
 double u[7][2]={{0}},k[14][14]={{0}},r[7][2]={{0}},v[7][3]={{0}},a[21][21]={{0}},b[7][3]={{0}};
 int i,j;
 tri7_stiff(k,r,xl,u,elastic,1);triox_stiff(a,b,xl,v,oxide,1);
 for(i=0;i<14;i++)for(j=0;j<14;j++)printf("%.17g%c",k[i][j],j==13?'\n':' ');
 for(i=0;i<21;i++)for(j=0;j<21;j++)printf("%.17g%c",a[i][j],j==20?'\n':' ');
 return 0;
}
