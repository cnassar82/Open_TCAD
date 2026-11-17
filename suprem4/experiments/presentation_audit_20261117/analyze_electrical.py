from pathlib import Path
import re,json,math,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent;original=out.parent/'directional_etch_20261117/electrical'
def parse(path,kind):
 cols={};rows=[]
 for line in path.read_text().splitlines():
  if m:=re.match(r'#\s+(\d+)\s+(.+)\s+\[([^]]+)\]',line):cols[m[2].strip()]=int(m[1])-1
  elif line and not line.startswith('#'):rows.append(list(map(float,line.split())))
 rows=np.array(rows);assert len(rows)==67 and np.isfinite(rows).all()
 indices=lambda p:[j for k,j in cols.items() if k.startswith(p)]
 g=indices('Vapp(Gate[');dv=indices('Vapp(Drain[')[0];di=indices('I(Drain[')[0];si=indices('I(Source[')[0];bi=indices('I(Body[')[0];gi=indices('I(Gate[')
 assert all(np.max(abs(rows[:,j]-rows[:,g[0]]))<1e-12 for j in g)
 x=abs(rows[:,g[0]] if kind=='transfer' else rows[:,dv]);current=abs(rows[:,di]);assert np.allclose(x,np.arange(67)*.05,atol=1e-12)
 kcl=np.max(abs(rows[:,di]+rows[:,si]+rows[:,bi]+rows[:,gi].sum(axis=1)))
 return x,current,{'points':67,'max_KCL_residual_A':float(kcl),'drain_first_V':float(rows[0,dv]),'drain_last_V':float(rows[-1,dv]),'gate_first_V':float(rows[0,g[0]]),'gate_last_V':float(rows[-1,g[0]])}
def crossing(x,y,target):
 if target<min(y) or target>max(y):return None
 for i in range(len(x)-1):
  if y[i]<=target<=y[i+1]:return float(x[i]+(target-y[i])/(y[i+1]-y[i])*(x[i+1]-x[i]))
 return None
results={};fig,axes=plt.subplots(2,2,figsize=(10,7.4),layout='constrained')
for row,dev in enumerate(('nmos','pmos')):
 sign=1 if dev=='nmos' else -1
 x,low,record=parse(original/dev/'transfer.dat','transfer');gm=np.diff(low)/np.diff(x);j=int(np.argmax(gm));vint=float(x[j]-low[j]/gm[j]);vt=sign*(vint-.025)
 lowcc=crossing(x,low,1e-7)
 highpath=out/'electrical'/dev/'transfer_vd3p3';execution=json.loads((highpath/'execution.json').read_text());assert execution['returncode']==0
 assert 'DIVERGED' not in (highpath/'stdout.log').read_text() and 'ERROR' not in (highpath/'stdout.log').read_text()
 _,high,highrecord=parse(highpath/'transfer.dat','transfer');assert abs(highrecord['drain_first_V']-sign*3.3)<1e-12 and abs(highrecord['drain_last_V']-sign*3.3)<1e-12
 highcc=crossing(x,high,1e-7)
 # Four-point subthreshold fit: report its actual current/bias window.
 windows=[]
 for k in range(len(x)-3):
  xx=x[k:k+4];yy=low[k:k+4]
  if min(yy)>=1e-10 and max(yy)<=1e-7 and yy[-1]/yy[0]>=10**.5:
   slope,intercept=np.polyfit(xx,np.log10(yy),1)
   if slope>0:windows.append({'SS_mV_per_dec':float(1000/slope),'gate_drive_range_V':xx[[0,-1]].tolist(),'current_range_A':yy[[0,-1]].tolist()})
 ss=min(windows,key=lambda w:w['SS_mV_per_dec']) if windows else None
 result={'threshold_method':'largest adjacent-point d|Id|/d|Vg| tangent; signed (x_intercept - |Vd|/2), low Vd=50mV; 50mV gate sampling','Vt_maxgm_V':vt,'Vt_target_V':sign*.75,'Vt_error_V':vt-sign*.75,'gm_max_A_per_V_um':float(gm[j]),'gm_peak_gate_interval_V':x[j:j+2].tolist(),'Vt_constant_100nA_V':None if lowcc is None else sign*lowcc,'Vt_constant_100nA_high_drain_V':None if highcc is None else sign*highcc,'DIBL_constant_100nA_mV_per_V':None if lowcc is None or highcc is None else (lowcc-highcc)/3.25*1000,'SS_best_four_point_low_drain':ss,'low_drain':{**record,'Ioff_A_per_um':float(low[0]),'Ion_A_per_um':float(low[-1])},'high_drain':{**highrecord,'Ioff_A_per_um':float(high[0]),'Ion_A_per_um':float(high[-1]),'execution':execution},'output_curves':{}}
 ax=axes[row,0];ax.semilogy(x,low*1e6,label='|Vds|=0.05 V');ax.semilogy(x,high*1e6,label='|Vds|=3.3 V');ax.axvline(.75,ls='--',color='black',lw=1,label='Slide |Vt| target');ax.axvline(abs(vt),ls=':',color='#bb4c29',lw=1,label='Max-gm estimate (low drain)');ax.set(title=dev.upper()+' transfer',xlabel='Gate drive |Vgs| (V)',ylabel='|Id| (uA / um)',xlim=(0,3.3));ax.grid(alpha=.2);ax.legend(fontsize=8)
 for gate in (.75,1.5,3.3):
  folder=out/'electrical'/dev/f'output_vg{gate:g}'.replace('.','p');ex=json.loads((folder/'execution.json').read_text());assert ex['returncode']==0
  text=(folder/'stdout.log').read_text();assert 'ERROR' not in text and 'DIVERGED' not in text
  xd,yd,rc=parse(folder/'output.dat','output');assert abs(rc['gate_first_V']-sign*gate)<1e-12 and abs(rc['gate_last_V']-sign*gate)<1e-12
  result['output_curves'][str(gate)]={**rc,'Id_at_3p3V_drain_A_per_um':float(yd[-1]),'execution':ex}
  axes[row,1].plot(xd,yd*1e6,label=f'|Vgs|={gate:g} V')
  if gate==3.3:result['Ion_path_relative_difference']=float(abs(yd[-1]-high[-1])/high[-1])
 axes[row,1].set(title=dev.upper()+' output',xlabel='Drain bias |Vds| (V)',ylabel='|Id| (uA / um)',xlim=(0,3.3));axes[row,1].grid(alpha=.2);axes[row,1].legend(fontsize=8)
 results[dev]=result
fig.suptitle('Actual corrected-structure electrical results; no fit to the slide thresholds\n300 K, 1 um width, ideal n+ gate, Lucent mobility, QF=0',fontsize=12)
fig.savefig(out/'figures/electrical_vs_targets.png',dpi=170);fig.savefig(out/'figures/electrical_vs_targets.svg')
(out/'electrical_metrics.json').write_text(json.dumps(results,indent=2)+'\n')
records=[json.loads(p.read_text()) for p in (out/'electrical').glob('*/*/execution.json')];assert len(records)==8 and all(r['returncode']==0 for r in records)
(out/'electrical/executions_verified.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(results,indent=2))
