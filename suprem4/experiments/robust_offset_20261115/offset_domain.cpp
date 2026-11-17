// Standalone domain-construction experiment, not wired into the simulator.
#include <clipper2/clipper.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace Clipper2Lib;
int main(int argc,char **argv) {
 try {
  double thickness=argc>1?std::stod(argv[1]):.5;
  PathD surface;double x,y;
  while(std::cin>>x>>y)surface.emplace_back(x,y);
  if(surface.size()<2 || thickness<=0)throw std::runtime_error("invalid input");
  double left=surface.front().x,right=surface.back().x;
  double bottom=surface.front().y,top=bottom;
  for(const auto &p:surface){bottom=std::max(bottom,p.y);top=std::min(top,p.y);}
  bottom+=10*thickness;
  PathD original=surface;
  original.emplace_back(right,bottom);original.emplace_back(left,bottom);
  PathsD occupied=Union(PathsD{original},FillRule::NonZero,8);
  if(occupied.size()!=1)throw std::runtime_error("source is not one simple occupied domain");
  // Arc chord error is declared separately from the integer grid resolution.
  auto grown=InflatePaths(occupied,thickness,JoinType::Round,EndType::Polygon,2,8,1e-5);
  PathsD window{{{left,top-thickness-1},{right,top-thickness-1},{right,bottom},{left,bottom}}};
  auto coat=Difference(Intersect(grown,window,FillRule::NonZero,8),occupied,FillRule::NonZero,8);
  double area=0;int holes=0;
  for(const auto &path:coat){double a=Area(path);area+=a;if(a<0)holes++;}
  std::cout<<std::setprecision(17)<<"{\"grid_um\":1e-8,\"arc_error_um\":1e-5,\"area_um2\":"<<area<<",\"holes\":"<<holes<<",\"paths\":[";
  bool first=true;
  for(const auto &path:coat){if(!first)std::cout<<',';first=false;std::cout<<'[';bool pfirst=true;
   for(const auto &p:path){if(!pfirst)std::cout<<',';pfirst=false;std::cout<<'['<<p.x<<','<<p.y<<']';}std::cout<<']';}
  std::cout<<"]}\n";
 }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
