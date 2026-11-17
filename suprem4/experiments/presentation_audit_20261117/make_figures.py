import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import matplotlib.tri as mtri
out=Path(__file__).resolve().parent;root=Path('/home/chris/Open_TCAD');sys.path.insert(0,str(root/'decks/subcmos2014'))
from validate import read
r=json.loads((out/'profile_metrics.json').read_text());a=json.loads((out/'profile_arrays.json').read_text())
figdir=out/'figures';figdir.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 fig.savefig(figdir/(name+'.png'),dpi=170);fig.savefig(figdir/(name+'.svg'));plt.close(fig)
def profile(dev,stage,loc):
 p=a[f'{dev}/{stage}/{loc}'];return np.array(p['depth_um']),np.array(p['B_cm3']),np.array(p['P_cm3'])
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for ax,dev,target,jtarget in zip(axes,('nmos','pmos'),(1e16,3e16),(3.,2.5)):
 x,b,p=profile(dev,'wells','channel');ax.semilogy(x,b,label='B chemical',color='#1769aa')
 if max(p)>0:ax.semilogy(x,p,label='P chemical',color='#b54e28')
 ax.axhline(target,color='black',ls='--',lw=1,label='Slide nominal well doping')
 ax.axvline(jtarget,color='#7c397d',ls='--',lw=1,label='Slide well-depth target')
 ax.axhline(1e15,color='#777',ls=':',lw=1,label='Starting substrate')
 if dev=='pmos':ax.axvline(r[dev]['wells']['cuts']['channel']['chemical_junction_depths_um'][-1],color='#b54e28',ls=':',label='Simulated PN junction')
 else:ax.axvline(r[dev]['wells']['cuts']['channel']['B_2e15_crossing_depths_um'][-1],color='#1769aa',ls=':',label='B = 2 x background marker')
 ax.set(xlim=(0,5),ylim=(1e13,3e17),xlabel='Depth below local silicon surface (um)',ylabel='Concentration (cm^-3)',title=dev.upper()+' after well drive');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Well profiles: nominal slide values do not define a measured profile')
save(fig,'wells_vs_targets')
fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
for ax,dev in zip(axes,('nmos','pmos')):
 for stage,color in [('wells','#777777'),('threshold','#1769aa'),('gate_oxide','#41976b'),('metal1','#bf5e25')]:
  x,b,p=profile(dev,stage,'channel');ax.semilogy(x*1000,b,color=color,label=stage+' B')
  if max(p)>0:ax.semilogy(x*1000,p,color=color,ls='--',label=stage+' P')
 ax.set(xlim=(0,400),ylim=(1e14,1e18),xlabel='Depth below local silicon surface (nm)',ylabel='Chemical concentration (cm^-3)',title=dev.upper()+' channel, x=6 um');ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2)
fig.suptitle('Channel evolution: solid boron, dashed phosphorus')
save(fig,'channel_evolution')
fig,axes=plt.subplots(1,2,figsize=(10,4.1),layout='constrained')
for ax,dev in zip(axes,('nmos','pmos')):
 for loc,color in [('source','#1769aa'),('source_extension','#bf5e25'),('channel','#777777')]:
  x,b,p=profile(dev,'metal1',loc);ax.semilogy(x,b,color=color,label=loc+' B');ax.semilogy(x,p,color=color,ls='--',label=loc+' P')
 ax.set(xlim=(0,1.2),ylim=(1e14,3e20),xlabel='Depth below local silicon surface (um)',ylabel='Chemical concentration (cm^-3)',title=dev.upper()+' final junction profiles');ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2)
fig.suptitle('Source, spacer extension and channel are different sampling locations')
save(fig,'final_junction_profiles')
helper=out/'extract_profiles.py';ns={'__file__':str(helper)};exec(helper.read_text().split('summary={};all_profiles={}')[0],ns)
fig,axes=plt.subplots(1,2,figsize=(10,4.1),layout='constrained')
for ax,dev in zip(axes,('nmos','pmos')):
 for label,path,color in [('Original 2D mesh',out/'poly_implant'/dev/'after.str','#b54e28'),('2D spacing=0.02 um',out/'poly_refinement'/f'{dev}_h0.02'/'after.str','#7f4b98')]:
  y,v=ns['cut'](read(path),6,4);ax.semilogy((y-y[0])*1000,v[:,1],'.-',label=label,color=color,ms=3)
 ref=np.genfromtxt(out/'poly_1d_reference/250/profile.csv',delimiter=',',names=True);ax.semilogy(ref['depth_um']*1000,ref['P_cm3'],color='#1769aa',label='Resolved planar 1D control')
 ax.axvline(66.3,color='black',ls=':',lw=1,label='Table projected range: 66.3 nm')
 ax.set(xlim=(0,250),ylim=(1e15,1e22),xlabel='Depth below poly top (nm)',ylabel='As-implanted phosphorus (cm^-3)',title=dev.upper()+' gate center');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Poly implant is under-resolved in 2D: no internal sample near the peak')
save(fig,'poly_implant_resolution')
fig,axes=plt.subplots(1,2,figsize=(11,4.1),layout='constrained');levels=[-1e20,-1e19,-1e18,-1e17,-1e16,-1e15,0,1e15,1e16,1e17,1e18,1e19,1e20]
for ax,dev in zip(axes,('nmos','pmos')):
 coords,tris,nodes=read(root/'decks/subcmos2014/results/directional_final_fresh_20261117'/dev/'metal1.str')
 ids=sorted({v for mat,vs in tris if mat==3 for v in vs});lookup={v:i for i,v in enumerate(ids)};pts=np.array([coords[v] for v in ids]);tri=np.array([[lookup[v] for v in vs] for mat,vs in tris if mat==3]);net=np.array([nodes[v,3].get(21,0)-nodes[v,3].get(23,0) for v in ids]);grid=mtri.Triangulation(pts[:,0],pts[:,1],tri)
 cf=ax.tricontourf(grid,net,levels=levels,norm=SymLogNorm(linthresh=1e15,vmin=-1e20,vmax=1e20),cmap='coolwarm',extend='both');ax.tricontour(grid,net,levels=[0],colors='black',linewidths=.6)
 ax.set(xlim=(0,10),ylim=(1.6,-.1),xlabel='Lateral position (um)',ylabel='Absolute depth coordinate (um)',title=dev.upper()+' isolated cut');ax.set_facecolor('#eeeeee')
 for x in (1.25,4.25,5.25,6,6.75,7.75):ax.axvline(x,lw=.45,color='#555',alpha=.35)
cb=fig.colorbar(cf,ax=axes,label='Active net donor concentration (cm^-3); blue=p, red=n',ticks=[-1e20,-1e18,-1e16,0,1e16,1e18,1e20],shrink=.85)
cb.set_ticklabels([r'$-10^{20}$',r'$-10^{18}$',r'$-10^{16}$','0',r'$10^{16}$',r'$10^{18}$',r'$10^{20}$'])
fig.suptitle('Final silicon doping and active junction locations (black); oxide/poly omitted')
save(fig,'final_doping_maps')
print(figdir)
