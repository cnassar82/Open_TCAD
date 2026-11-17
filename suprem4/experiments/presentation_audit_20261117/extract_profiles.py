import sys,json,hashlib,csv
from pathlib import Path
import numpy as np
root=Path('/home/chris/Open_TCAD');sys.path.insert(0,str(root/'decks/subcmos2014'))
from validate import read,intervals
source=root/'decks/subcmos2014/results/directional_final_fresh_20261117'
out=Path(__file__).resolve().parent;profiles=out/'profiles';profiles.mkdir(exist_ok=True)
stages=['initial','pad','well_mask','wells','locos_mask','locos','threshold','gate_oxide','poly_doped','gate','ldd','spacers','source_drain','annealed','contacts','metal1']
locations={'body':1.25,'field_stop':2.5,'source':4.25,'source_extension':5.25,'channel':6.,'drain_extension':6.75,'drain':7.75}
streams=[5,3,23,21]
def cut(mesh,x,mat=3):
 coords,tris,nodes=mesh;vs=np.array([v for m,v in tris if m==mat],dtype=int)
 if not len(vs):return None
 points=np.array([[coords[v] for v in tri] for tri in vs]);values=np.array([[[nodes[v,mat].get(s,0. if s in (3,5) or {23:5,21:3}.get(s) not in nodes[v,mat] else np.nan) for s in streams] for v in tri] for tri in vs])
 yy=[];vv=[]
 for a,b in [(0,1),(1,2),(2,0)]:
  xa,xb=points[:,a,0],points[:,b,0];sel=(np.minimum(xa,xb)<=x)&(np.maximum(xa,xb)>=x)&(xa!=xb)
  alpha=(x-xa[sel])/(xb[sel]-xa[sel]);yy.extend(points[sel,a,1]+alpha*(points[sel,b,1]-points[sel,a,1]));vv.extend(values[sel,a]+alpha[:,None]*(values[sel,b]-values[sel,a]))
  sel=(xa==xb)&(xa==x)
  for k in (a,b):yy.extend(points[sel,k,1]);vv.extend(values[sel,k])
 if not yy:return None
 yy=np.array(yy);vv=np.array(vv);order=np.argsort(yy,kind='stable');yy=yy[order];vv=vv[order]
 yy,index=np.unique(yy,return_index=True);vv=vv[index]
 return yy,vv

def roots(y,net):
 result=[]
 for i in range(len(y)-1):
  if net[i]==0:result.append(float(y[i]))
  elif net[i]*net[i+1]<0:result.append(float(y[i]-net[i]*(y[i+1]-y[i])/(net[i+1]-net[i])))
 return result
summary={};all_profiles={}
for dev in ('nmos','pmos'):
 summary[dev]={}
 for stage in stages:
  path=source/dev/(stage+'.str');mesh=read(path);result={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'cuts':{}}
  for label,x in locations.items():
   y,v=cut(mesh,x);depth=y-y[0];net=v[:,1]-v[:,0]
   profile={'y_um':y.tolist(),'depth_um':depth.tolist(),'B_cm3':v[:,0].tolist(),'P_cm3':v[:,1].tolist(),'B_active_cm3':v[:,2].tolist(),'P_active_cm3':v[:,3].tolist()}
   key=f'{dev}/{stage}/{label}';all_profiles[key]=profile
   f=profiles/f'{dev}_{stage}_{label}.csv'
   with f.open('w') as fh:
    w=csv.writer(fh);w.writerow(['depth_below_local_Si_um','absolute_y_um','B_chemical_cm3','P_chemical_cm3','B_active_cm3','P_active_cm3','net_chemical_donor_cm3','net_active_donor_cm3'])
    for d,abs_y,row in zip(depth,y,v):w.writerow([d,abs_y,*row,row[1]-row[0],row[3]-row[2]])
   dose=lambda j:float(np.trapezoid(v[:,j],y)*1e-4)
   sample=lambda d,j:float(np.interp(d,depth,v[:,j]))
   result['cuts'][label]={'x_um':x,'surface_y_um':float(y[0]),'B_sheet_cm2':dose(0),'P_sheet_cm2':dose(1),'B_excess_above_1e15_sheet_cm2':float(np.trapezoid(v[:,0]-1e15,y)*1e-4),'peak_B_cm3':float(v[:,0].max()),'peak_P_cm3':float(v[:,1].max()),'chemical_junction_depths_um':roots(depth,net),'active_junction_depths_um':roots(depth,v[:,3]-v[:,2]),'B_2e15_crossing_depths_um':roots(depth,v[:,0]-2e15),'samples':{str(d):{'B':sample(d,0),'P':sample(d,1),'Ba':sample(d,2),'Pa':sample(d,3),'net_chemical':sample(d,1)-sample(d,0),'net_active':sample(d,3)-sample(d,2)} for d in (0,.01,.05,.1,.5,1,2,3,5)}}
  result['film_spans_um']={str(x):{str(m):intervals(mesh,x,m) for m in (1,2,3,4,6)} for x in (1.25,2.5,4.25,6.,7.75)}
  if any(m==4 for m,vs in mesh[1]):
   y,v=cut(mesh,6,4);result['poly_channel']={'thickness_um':float(y[-1]-y[0]),'B_sheet_cm2':float(np.trapezoid(v[:,0],y)*1e-4),'P_sheet_cm2':float(np.trapezoid(v[:,1],y)*1e-4),'mid_B_cm3':float(np.interp((y[0]+y[-1])/2,y,v[:,0])),'mid_P_cm3':float(np.interp((y[0]+y[-1])/2,y,v[:,1])),'mid_B_active_cm3':float(np.interp((y[0]+y[-1])/2,y,v[:,2])),'mid_P_active_cm3':float(np.interp((y[0]+y[-1])/2,y,v[:,3]))}
  summary[dev][stage]=result
  print(dev,stage,flush=True)
(out/'profile_metrics.json').write_text(json.dumps(summary,indent=2,allow_nan=True)+'\n')
(out/'profile_arrays.json').write_text(json.dumps(all_profiles,allow_nan=True)+'\n')
