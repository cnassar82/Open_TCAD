// Analytic coating domains, independent of the legacy offset implementation.
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
extern "C" int suprem_domain_offset(int,const double*,double,double,double**);
extern "C" int suprem_nominal_offset(int,const double*,double,double,double**);
static double area(const std::vector<double> &xy) {
 double a=0;size_t n=xy.size()/2;
 for(size_t i=0;i<n;i++){size_t j=(i+1)%n;a+=xy[2*i]*xy[2*j+1]-xy[2*j]*xy[2*i+1];}
 return std::fabs(a)*.5;
}
static void check(std::vector<double> surface,double expected,double shift) {
 for(size_t i=0;i<surface.size();i+=2)surface[i]=float(surface[i]+shift);
 for(size_t i=1;i<surface.size();i+=2)surface[i]=float(surface[i]);
 auto saved=surface;double *output=nullptr;
 int n=suprem_domain_offset(surface.size()/2,surface.data(),.5e-4,.02e-4,&output);
 assert(n>=2 && output && surface==saved);
 assert(output[0]==surface[0] && output[2*n-2]==surface[surface.size()-2]);
 for(int i=0;i<2*n;i++)assert(double(float(output[i]))==output[i]);
 auto ring=surface;
 for(int i=n-1;i>=0;i--){ring.push_back(output[2*i]);ring.push_back(output[2*i+1]);}
 for(double v:ring)assert(std::isfinite(v));
 assert(std::fabs(area(ring)*1e8-expected)<3e-5);
 std::free(output);
}
int main(int argc,char **argv) {
 for(double height:{0.,.1e-4,-.427e-4})for(double width:{4e-4,12e-4,100e-4})
  for(double micron:{.37,.5,1.35}) {
   double surface[]={0,double(float(height)),double(float(width)),double(float(height))};
   double thickness=double(float(micron))*1e-4,*construction=nullptr;
   int n=suprem_nominal_offset(2,surface,thickness,.04e-4,&construction);
   assert(n==2 && construction);
   for(int i=0;i<n;i++)assert(float(construction[2*i+1]+thickness)==float(height));
   std::free(construction);
  }
 for(double shift:{0.,-.000137,.000631}) {
  check({0,0,4e-4,0},2,shift);
  check({0,0,2e-4,0,2e-4,1e-4,4e-4,1e-4},2.5+(std::acos(-1.)/4-1)*.25,shift);
 }
 double *output=nullptr;
 // Unsupported topologies are reported before touching the caller's domain.
 double reversed[]={1,0,0,0};
 assert(suprem_domain_offset(2,reversed,1,.1,&output)<0 && !output);
 assert(argc==2);
 std::ifstream fixture(argv[1]);assert(fixture);
 std::vector<double> surface;std::string line;
 while(std::getline(fixture,line)) {
  if(line.empty() || line[0]=='#')continue;
  std::istringstream row(line);double x,y;assert(row>>x>>y);
  surface.push_back(x);surface.push_back(y);
 }
 assert(surface.size()==286);
 auto saved=surface;
 int count=suprem_domain_offset(surface.size()/2,surface.data(),double(float(.15e-4)),.04e-4,&output);
 assert(count>=2 && output && surface==saved);
 // Dilation of this shallow, single-valued surface is single-valued too.
 // The old construction inserted a one-grid backwards horizontal edge.
 for(int i=1;i<count;i++)assert(output[2*i]>=output[2*(i-1)]);
 std::free(output);
 puts("PASS domain offset: planar/step analytic areas, translations, side constraints and failure contract");
}
