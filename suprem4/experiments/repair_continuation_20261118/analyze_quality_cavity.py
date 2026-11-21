#!/usr/bin/env python3
"""Enumerate all legal triangulations of a small captured removal cavity."""
import argparse,hashlib,json,math,struct,sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path[:0]=[str(ROOT/'decks/subcmos2014'),str(ROOT/'suprem4/tests')]
from validate import read
from mesh_patch_oracle import propose_interior_star,orient,inside_or_on
p=argparse.ArgumentParser(description=__doc__);p.add_argument('structure',type=Path);p.add_argument('--point',type=int,required=True,help='1-based exported point ID');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
coords,tris,_=read(a.structure);old=[vs for mat,vs in tris if a.point in vs];materials={mat for mat,vs in tris if a.point in vs}
if len(materials)!=1:raise ValueError('single-material diagnostic required')
f32=lambda x:struct.unpack('f',struct.pack('f',x*1e-4))[0]
points={v:tuple(Fraction(f32(x)) for x in coords[v]) for t in old for v in t}
oracle=propose_interior_star(a.point,points,old)
if oracle.status!='OK':raise ValueError(oracle.reason)
order=oracle.boundary
if len(order)>12:raise ValueError('exhaustive diagnostic limited to 12 boundary vertices')
def quality(t):
 ps=[points[v] for v in t];den=sum((u[0]-v[0])**2+(u[1]-v[1])**2 for u,v in zip(ps,ps[1:]+ps[:1]))
 return 2*math.sqrt(3)*float(orient(points,*t)/den)
def ears(poly):
 for i,b in enumerate(poly):
  t=(poly[i-1],b,poly[(i+1)%len(poly)])
  if orient(points,*t)>0 and not any(inside_or_on(points,v,*t) for v in poly if v not in t):yield i,t
@lru_cache(None)
def enumerate_triangles(poly):
 if len(poly)==3:return ((tuple(poly),),) if orient(points,*poly)>0 else ()
 result=set()
 for i,t in ears(poly):
  for rest in enumerate_triangles(poly[:i]+poly[i+1:]):result.add(tuple(sorted((tuple(sorted(t)),)+tuple(tuple(sorted(x)) for x in rest))))
 # Canonicalize triangle identity above; orient each recovered face CCW.
 oriented=[]
 for tri in result:
  oriented.append(tuple(x if orient(points,*x)>0 else (x[0],x[2],x[1]) for x in tri))
 return tuple(oriented)
all_t=enumerate_triangles(tuple(order));best=max(all_t,key=lambda ts:min(map(quality,ts)))
work=list(order);greedy=[]
while len(work)>3:
 i,t=max(ears(work),key=lambda pair:quality(pair[1]));greedy.append(t);work.pop(i)
greedy.append(tuple(work))
r={'accepted_full_cmos':False,'structure_sha256':hashlib.sha256(a.structure.read_bytes()).hexdigest(),'point':a.point,'material':next(iter(materials)),'boundary_points':list(order),'old_triangles':old,'triangulations':len(all_t),'old_min_quality':min(map(quality,old)),'greedy_min_quality':min(map(quality,greedy)),'optimal_min_quality':min(map(quality,best)),'greedy_triangles':greedy,'optimal_triangles':best,'coordinates_cm':{str(v):[float(x) for x in xy] for v,xy in points.items()}}
a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
