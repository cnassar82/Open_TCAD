"""Read-only feasibility test: triangulations valid before AND after motion.
Optional outputs contain geometry only; they are not process results.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"tests"))
from partition_audit import orientation,audit_partition,_pairs,_bounds
source=Path(sys.argv[1]);dt=float(sys.argv[2]);coords={};nodes={};regions={};faces={};species=[];outer={}
for l in source.read_text().splitlines():
 s=l.split()
 if not s:continue
 if s[0]=='c':coords[int(s[1])]=np.array(list(map(float,s[2:4])))*1e-4
 elif s[0]=='r':regions[int(s[1])]=int(s[2])
 elif s[0]=='s':species=list(map(int,s[2:]))
 elif s[0]=='n':nodes.setdefault(int(s[1])+1,{})[int(s[2])]=dict(zip(species,map(float,s[3:])))
 elif s[0]=='t':
  vs=list(map(int,s[3:6]));faces[int(s[1])]=[int(s[2]),vs]
  for j,n in enumerate(map(int,s[6:9])):
   if n<0:outer[tuple(sorted((vs[(j+1)%3],vs[(j+2)%3])))]=n
velocity={}
for i,ns in nodes.items():
 m=3 if 3 in ns else 4 if 4 in ns else next(iter(ns));velocity[i]=np.array([ns[m].get(8,0),ns[m].get(9,0)])
future={i:(c+velocity[i]*dt).astype(np.float32).astype(float) for i,c in coords.items()}
def area(v,xy):
 a,b,c=[xy[i] for i in v];u=b-a;w=c-a;return float(u[0]*w[1]-u[1]*w[0])
def side(v,xy):return orientation(*(tuple(map(float,xy[i])) for i in v))
def signed(v):return v if side(v,coords)>0 else [v[0],v[2],v[1]]
def edges(v):return [tuple(sorted((v[i],v[(i+1)%3]))) for i in range(3)]
incidence={}
for t,(r,v) in faces.items():
 for e in edges(v):incidence.setdefault(e,set()).add(t)
alignment=[]
if '--align' in sys.argv:
 for iteration in range(100):
  records=[]
  for e,owners in incidence.items():
   if len(owners)==1 or len({faces[t][0] for t in owners})>1:
    records.append(_bounds([tuple(map(float,future[i])) for i in e],e))
  split=None
  for first,second in _pairs(records,True):
   a,b=first[5];c,d=second[5]
   if orientation(a,b,c)*orientation(a,b,d)>=0 or orientation(c,d,a)*orientation(c,d,b)>=0:continue
   for coarser,finer in [(first[4],second[4]),(second[4],first[4])]:
    owners=incidence[coarser]
    if len(owners)!=2:continue
    materials={regions[faces[t][0]] for t in owners}
    # A silicon reaction interface has different kinematics. This prototype
    # only studies the shared vertical overburden map on oxide/poly faces.
    if materials!={1,4}:continue
    aa,bb=sorted(coarser,key=lambda i:coords[i][0])
    for point in finer:
     if coords[aa][0]<coords[point][0]<coords[bb][0]:split=(aa,bb,point);break
    if split:break
   if split:break
  if not split:break
  a,b,sample=split;t=(coords[sample][0]-coords[a][0])/(coords[b][0]-coords[a][0]);point=max(coords)+1
  xy=(coords[a]+t*(coords[b]-coords[a])).astype(np.float32).astype(float)
  xy[0]=coords[sample][0]
  coords[point]=xy;velocity[point]=velocity[sample].copy();future[point]=(xy+velocity[point]*dt).astype(np.float32).astype(float)
  owners=list(incidence[tuple(sorted((a,b)))]);replacement=[]
  for owner in owners:
   r,v=faces[owner];opp=next(i for i in v if i not in [a,b]);left=signed([opp,a,point]);right=signed([opp,point,b])
   if min(side(left,coords),side(right,coords))<=0:raise ValueError('rounded interface subdivision invalid at input coordinates')
   replacement.append((owner,r,left,right))
  for owner,r,left,right in replacement:
   for e in edges(faces[owner][1]):
    incidence[e].remove(owner)
    if not incidence[e]:del incidence[e]
   faces[owner][1]=left;extra=max(faces)+1;faces[extra]=[r,right]
   for tid in [owner,extra]:
    for e in edges(faces[tid][1]):incidence.setdefault(e,set()).add(tid)
  projected=coords[a]+t*(coords[b]-coords[a]);alignment.append({'edge':[a,b],'sample':sample,'new_point':point,'construction_error_cm':float(np.linalg.norm(xy-projected))})
 else:raise ValueError('alignment iteration bound')
bad={t for t,(r,v) in faces.items() if side(v,future)<=0};initial=len(bad);plan=[]
while bad:
 changed=False
 for t in list(bad):
  r,v=faces[t]
  for e in edges(v):
   owners=incidence[e]
   if len(owners)!=2:continue
   u=next(i for i in owners if i!=t)
   if faces[u][0]!=r:continue
   a=next(i for i in v if i not in e);b,c=e;d=next(i for i in faces[u][1] if i not in e)
   if a==d or tuple(sorted((a,d))) in incidence:continue
   left=signed([a,b,d]);right=signed([a,d,c])
   if min(side(left,coords),side(right,coords),side(left,future),side(right,future))<=0:continue
   for q in [t,u]:
    for edge in edges(faces[q][1]):
     incidence[edge].remove(q)
     if not incidence[edge]:del incidence[edge]
   faces[t][1]=left;faces[u][1]=right
   for q in [t,u]:
    for edge in edges(faces[q][1]):incidence.setdefault(edge,set()).add(q)
   bad.discard(t);bad.discard(u);plan.append([t,u,list(e),[a,d]])
   changed=True;break
  if changed:break
 if not changed:break
report={'source':str(source),'dt':dt,'initial_bad':initial,'remaining_bad':len(bad),'flips':len(plan),'unresolved':sorted(bad)[:20],'plan':plan,'alignment':alignment,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'passed_geometry':False}
if not bad and len(sys.argv)>3:
 out=Path(sys.argv[3]);outer={}
 for line in source.read_text().splitlines():
  s=line.split()
  if not s or s[0]!='t':continue
  vs=list(map(int,s[3:6]));ns=list(map(int,s[6:9]))
  for j,n in enumerate(ns):
   if n<0:outer[tuple(sorted((vs[(j+1)%3],vs[(j+2)%3])))]=n
 lines=['# DIAGNOSTIC geometry only: no chemical fields, not process acceptance']
 lines += [f'c {i} {x*1e4:.17g} {y*1e4:.17g} 0' for i,(x,y) in future.items()]
 lines += [f'r {r} {m}' for r,m in regions.items()]
 for t,(r,v) in faces.items():
  ns=[]
  for j in range(3):
   edge=tuple(sorted((v[(j+1)%3],v[(j+2)%3])));owners=incidence[edge]-{t}
   ns.append(next(iter(owners)) if owners else outer[edge])
  lines.append(' '.join(map(str,['t',t,r,*v,*ns,-1,-1])))
 out.write_text('\n'.join(lines)+'\n')
 try:report['partition_audit']=audit_partition(out);report['passed_geometry']=True
 except ValueError as error:report['passed_geometry']=False;report['partition_error']=str(error)
print(json.dumps(report,indent=2))
