#!/usr/bin/env python3
"""The same positive poly/Si distance can contain oxide or an invalid void."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'decks/subcmos2014'))
from validate import gate_dielectric

def mesh(oxide_top):
 coords={};tris=[]
 def rectangle(x0,x1,y0,y1,mat):
  n=len(coords);coords.update({n:(x0,y0),n+1:(x1,y0),n+2:(x1,y1),n+3:(x0,y1)})
  tris.extend([(mat,(n,n+1,n+2)),(mat,(n,n+2,n+3))])
 rectangle(0,10,.1,2,3);rectangle(5.5,6.5,-.4,.08,4)
 if oxide_top is not None:rectangle(5.5,6.5,.08,oxide_top,1)
 return coords,tris,{}
assert len(gate_dielectric(mesh(.1)))==5
for top in (None,.09):
 try:gate_dielectric(mesh(top))
 except ValueError as e:assert 'missing oxide material' in str(e)
 else:raise AssertionError('positive empty gap accepted as a gate dielectric')
coords={};tris=[];identities={}
def rectangle(x0,x1,y0,y1,mat):
 vs=[]
 for point in ((x0,y0),(x1,y0),(x1,y1),(x0,y1)):
  if point not in identities:
   index=len(coords);identities[point]=index;coords[index]=point
  vs.append(identities[point])
 tris.extend([(mat,tuple(vs[:3])),(mat,(vs[0],vs[2],vs[3]))])
rectangle(0,10,.1,2,3)
for left,right in ((5.5,5.51),(5.51,5.52),(5.52,6.5)):
 rectangle(left,right,-.4,.08,4)
 if left!=5.51:rectangle(left,right,.08,.1,1)
try:gate_dielectric((coords,tris,{}))
except ValueError as e:assert 'gate underside edge' in str(e)
else:raise AssertionError('gap between the five probes accepted')
print('PASS gate oxide occupancy: intact film accepted; empty, partial, and between-probe gaps rejected')
