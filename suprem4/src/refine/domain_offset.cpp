/* Domain-based conformal deposition using Clipper2 1.5.4 (Boost license).
 * C ABI: all coordinates are cm. Caller owns the malloc-allocated result.
 * Original interface vertices remain in the C database without modification.
 */
#include "../third_party/clipper2/include/clipper2/clipper.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <vector>
using namespace Clipper2Lib;

// Certify simplification against every vertex of each removed subchain.
// A local deletion tolerance alone does not bound accumulated deletions.
// Each chord joins retained vertices, so this also bounds the continuous
// polygonal subchain's distance from that finite chord.
static double simplification_error(const Path64 &before,const Path64 &after) {
 if(before.empty() || after.size()<3)return INFINITY;
 size_t cursor=0;
 while(cursor<before.size() && before[cursor]!=after[0])++cursor;
 if(cursor==before.size())return INFINITY;
 double worst=0;
 for(size_t j=0;j<after.size();++j) {
  auto a=after[j],b=after[(j+1)%after.size()];
  double dx=double(b.x)-a.x,dy=double(b.y)-a.y,den=dx*dx+dy*dy;
  if(den==0)return INFINITY;
  size_t walked=0;
  do {
   auto p=before[cursor];
   double x=double(p.x)-a.x,y=double(p.y)-a.y;
   double t=std::max(0.,std::min(1.,(x*dx+y*dy)/den));
   worst=std::max(worst,std::hypot(x-t*dx,y-t*dy));
   cursor=(cursor+1)%before.size();
   if(++walked>before.size())return INFINITY;
  }while(before[cursor]!=b);
 }
 return worst;
}
static int offset_domain(int count,const double *xy,double thickness,
                         double space,double **result,bool nominal) {
 *result=nullptr;
 try {
  if(count<2 || count>INT_MAX/2 || !xy || !(thickness>0) || !(space>0) ||
     !std::isfinite(thickness) || !std::isfinite(space) || thickness>100)return -1;
  double extent=thickness;
  for(int i=0;i<2*count;i++) {
   if(!std::isfinite(xy[i]) || std::fabs(xy[i])>100)return -1;
   extent=std::max(extent,std::fabs(xy[i])+thickness);
  }
  int exponent=0;std::frexp(extent,&exponent);
  // Construct the outer boundary on a dyadic grid exactly representable in
  // the database's binary32 coordinates. Do not round a fine-grid polygon
  // afterward: that can merge vertices and reverse short edges.
  double inverse=std::ldexp(1.0,std::max(-149,exponent-(nominal?48:24))),scale=1.0/inverse;
  Path64 source;
  for(int i=0;i<count;i++) {
   source.emplace_back(std::llround(xy[2*i]*scale),std::llround(xy[2*i+1]*scale));
  }
  const auto grid_left=source.front().x,grid_right=source.back().x;
  double left=xy[0],right=xy[2*count-2];
  if(!(left<right))return -1;
  double bottom=source.front().y*inverse,top=bottom;
  for(const auto &p:source){bottom=std::max(bottom,p.y*inverse);top=std::min(top,p.y*inverse);}
  double delta=thickness;
  bottom+=10*delta+right-left;
  Path64 original=source;
  original.emplace_back(grid_right,int64_t(std::llround(bottom*scale)));original.emplace_back(grid_left,int64_t(std::llround(bottom*scale)));
  auto occupied=Union(Paths64{original},FillRule::NonZero);
  if(occupied.size()!=1 || Area(occupied[0])<=0)return -2;
  double arc_error=std::min(space*space/(8*thickness),thickness*1e-4);
  // Offsetting a finely subdivided surface can generate redundant one-grid
  // edges, including tiny backwards joins. They are discretization artifacts
  // of the new surface, not features of the original material interface.
  // Allocate the existing chord-error budget between arc construction and
  // simplification; never simplify or replace the occupied input domain.
  auto grown=InflatePaths(occupied,delta*scale,JoinType::Round,EndType::Polygon,2,.5*arc_error*scale);
  auto simplified=SimplifyPaths(grown,nominal?2.:.25*arc_error*scale);
  double simplify_error=0;
  for(size_t i=0;i<grown.size();++i)
   simplify_error=std::max(simplify_error,simplification_error(grown[i],simplified[i]));
  if(simplify_error>.5*arc_error*scale)return -8;
  grown=std::move(simplified);
  int64_t bottom_grid=std::llround(bottom*scale),top_grid=std::llround((top-delta-thickness)*scale);
  Paths64 window{{{grid_left,top_grid},{grid_right,top_grid},{grid_right,bottom_grid},{grid_left,bottom_grid}}};
  auto domain=Intersect(grown,window,FillRule::NonZero);
  auto coat=Difference(domain,occupied,FillRule::NonZero);
  if(domain.size()!=1 || coat.size()!=1 || Area(coat[0])<=0)return -3;
  auto &ring=domain[0];
  // The clipping slab has exact side identities on the declared integer grid.
  int a=-1,b=-1;
  for(size_t i=0;i<ring.size();i++) {
   if(ring[i].x==grid_left && (a<0 || ring[i].y<ring[a].y))a=i;
   if(ring[i].x==grid_right && (b<0 || ring[i].y<ring[b].y))b=i;
  }
  if(a<0 || b<0 || a==b)return -4;
  std::vector<Point64> paths[2];
  for(int direction=0;direction<2;direction++) {
   int i=a,increment=direction?-1:1;
   paths[direction].push_back(ring[i]);
   for(size_t j=0;i!=b && j<ring.size();j++) {
    i=(i+increment+ring.size())%ring.size();paths[direction].push_back(ring[i]);
   }
  }
  auto maxy=[](const std::vector<Point64> &p){double y=p[0].y;for(auto &v:p)y=std::max(y,double(v.y));return y;};
  int chosen=maxy(paths[0])<maxy(paths[1])?0:1;
  const auto &outer=paths[chosen];
  if(outer.size()>size_t(INT_MAX/2))return -6;
  if(maxy(outer)>=(bottom-delta)*scale)return -5;
  *result=static_cast<double*>(std::malloc(2*outer.size()*sizeof(double)));
  if(!*result)return -6;
  for(size_t i=0;i<outer.size();i++) {
   double x=outer[i].x*inverse;
   if(outer[i].x==grid_left)x=left;
   if(outer[i].x==grid_right)x=right;
   double y=outer[i].y*inverse;
   if(nominal) {
    // Recover the analytic supporting plane of a horizontal generating face.
    // This is integer feature identity, not a distance-based material deletion.
    // It keeps (y - thickness) + thickness identical to the source plane.
    for(int j=0;j<count-1;j++) {
     if(xy[2*j+1]!=xy[2*j+3] || source[j].y!=source[j+1].y)continue;
     if(outer[i].x<std::min(source[j].x,source[j+1].x) ||
        outer[i].x>std::max(source[j].x,source[j+1].x))continue;
     if(outer[i].y==source[j].y-std::llround(delta*scale)) {
      y=xy[2*j+1]-thickness;break;
     }
    }
   }
   (*result)[2*i]=x;(*result)[2*i+1]=y;
  }
  fprintf(stderr,"[domain-offset%s] source=%d outer=%zu coat_area_um2=%.17g grid_um=%.17g arc_error_um=%.17g simplify_error_um=%.17g\n",
          nominal?"-construction":"",count,outer.size(),Area(coat[0])*inverse*inverse*1e8,inverse*1e4,arc_error*1e4,
          simplify_error*inverse*1e4);
  return outer.size();
 } catch(...) {std::free(*result);*result=nullptr;return -7;}
}
extern "C" int suprem_domain_offset(int count,const double *xy,double thickness,
                                    double space,double **result) {
 return offset_domain(count,xy,thickness,space,result,false);
}
extern "C" int suprem_nominal_offset(int count,const double *xy,double thickness,
                                     double space,double **result) {
 return offset_domain(count,xy,thickness,space,result,true);
}
