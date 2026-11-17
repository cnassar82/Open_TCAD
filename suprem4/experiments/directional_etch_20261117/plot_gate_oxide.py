import sys,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
root=Path('/home/chris/Open_TCAD');sys.path.insert(0,str(root/'decks/subcmos2014'))
from validate import read,intervals,gate_dielectric
source=root/'decks/subcmos2014/results/directional_final_fresh_20261117'
out=Path(__file__).resolve().parent
colors={1:'#86bfdf',3:'#d79265',4:'#d7bd65',6:'#b1b8c0'}
fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
for ax,dev in zip(axes,('nmos','pmos')):
 mesh=read(source/dev/'metal1.str');checks=gate_dielectric(mesh)
 coords,tris,_=mesh;surface=intervals(mesh,6,3)[0][0]
 patches=[];facecolors=[]
 for mat,vs in tris:
  ps=[coords[v] for v in vs]
  if max(p[0] for p in ps)<5.45 or min(p[0] for p in ps)>6.55:continue
  ys=[(p[1]-surface)*1000 for p in ps]
  if max(ys)<-80 or min(ys)>50:continue
  patches.append([(p[0],y) for p,y in zip(ps,ys)]);facecolors.append(colors.get(mat,'#cccccc'))
 ax.add_collection(PolyCollection(patches,facecolors=facecolors,edgecolors='none',antialiaseds=False))
 ax.set_xlim(5.45,6.55);ax.set_ylim(50,-80)
 ax.set_xlabel('Lateral position (um)');ax.set_ylabel('Depth relative to channel surface (nm)')
 ax.set_title(f'{dev.upper()}: center gate oxide {checks[2]["oxide_thickness_um"]*1000:.2f} nm')
 ax.text(6,-40,'Polysilicon gate',ha='center',va='center',fontsize=10)
 ax.text(6,-8,'Gate oxide',ha='center',va='center',fontsize=9)
 ax.text(6,28,'Silicon channel',ha='center',va='center',fontsize=10)
fig.suptitle('Corrected final structures after contact etch and first metal',fontsize=12)
fig.savefig(out/'final_gate_oxide.png',dpi=180)
fig.savefig(out/'final_gate_oxide.svg')
print(out/'final_gate_oxide.png')
