import sys,random,math,json,hashlib
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,'suprem4/tests')
from mesh_patch_oracle import propose_interior_star
# Independent coverage calculation: exact convex polygon clipping, not the
# proposal's incidence/ear/area/overlap predicates.
def det(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def area(poly):
 if len(poly)<3:return F(0)
 return abs(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(poly,poly[1:]+poly[:1])))/2
def clip(poly,tri):
 for a,b in zip(tri,tri[1:]+tri[:1]):
  out=[]
  if not poly:return []
  for s,e in zip(poly,poly[1:]+poly[:1]):
   ds,de=det(a,b,s),det(a,b,e)
   if (ds<0)!=(de<0):
    f=ds/(ds-de);out.append(tuple(s[k]+f*(e[k]-s[k]) for k in range(2)))
   if de>=0:out.append(e)
  poly=out
 return poly
rng=random.Random(20261115)
cases=0;outputs=0
for n in range(4,13):
 for case in range(4):
  pts={0:(F(0),F(0))}
  for i in range(n):
   angle=2*math.pi*(i+rng.uniform(-.2,.2))/n;r=rng.randint(100,1000)
   pts[i+1]=(F(round(r*math.cos(angle))),F(round(r*math.sin(angle))))
  old=[(0,i+1,(i+1)%n+1) for i in range(n)]
  assert all(det(*(pts[v] for v in t))>0 for t in old)
  for perm in range(2):
   ids=list(pts);rng.shuffle(ids);rename=dict(zip(pts,ids))
   p={rename[k]:v for k,v in pts.items()};t=[tuple(rename[v] for v in tr) for tr in old];rng.shuffle(t)
   before=repr((p,t));result=propose_interior_star(rename[0],p,t)
   assert before==repr((p,t))
   assert result.status=='OK',result
   oldpolys=[[p[v] for v in tr] for tr in t]
   newpolys=[[p[v] for v in tr] for tr in result.triangles]
   for new in newpolys:
    assert sum(area(clip(new[:],old)) for old in oldpolys)==area(new)
   for oldpoly in oldpolys:
    assert sum(area(clip(oldpoly[:],new)) for new in newpolys)==area(oldpoly)
   cases+=1;outputs+=len(newpolys)
print(json.dumps({'cases':cases,'candidate_triangles':outputs,'exact_bidirectional_coverage':'PASS','input_unchanged':'PASS','seed':20261115,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}))
