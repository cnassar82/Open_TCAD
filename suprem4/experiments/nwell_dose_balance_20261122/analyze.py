"""Phosphorus dose per material along vertical cuts (cm^-2), n-well dose balance."""
import importlib.util,json,sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('v','/home/chris/Open_TCAD/decks/subcmos2014/validate.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
here=Path(__file__).resolve().parent
NAMES={1:'oxide',2:'nitride',3:'silicon'}
def dose(mesh,x,mat):
 coords,tris,nodes=mesh
 sub=[(mm,vs) for mm,vs in tris if mm==mat and min(coords[q][0] for q in vs)<=x+1e-9 and max(coords[q][0] for q in vs)>=x-1e-9]
 mesh=(coords,sub,nodes);tot=0
 for a,b in v.intervals(mesh,x,mat):
  n=max(20,int((b-a)/5e-4));h=(b-a)/n
  for i in range(n):tot+=v.concentration(mesh,x,a+(i+.5)*h,mat)['phosphorus_cm3']*h*1e-4
 return tot
out={}
for case in ('through_nitride','after_strip'):
 for st in ('implanted','stripped','driven'):
  m=v.read(here/case/f'{st}.str')
  for x in (5.0,1.3):
   out[f'{case}/{st}/x={x}']={NAMES[k]:dose(m,x,k) for k in NAMES}
m=v.read(Path('/home/chris/Open_TCAD/decks/advcmos2014/results/full_wide_20261121/pmos/wells.str'))
out['full_flow_wells.str/x=5.0']={NAMES[k]:dose(m,5.0,k) for k in NAMES}
json.dump(out,open(here/'dose_balance.json','w'),indent=1)
for k,d in out.items():print(f'{k:40s}',' '.join(f'{n} {d[n]:.3e}' for n in d),f'total {sum(d.values()):.3e}')
