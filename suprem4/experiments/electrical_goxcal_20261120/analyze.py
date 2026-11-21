"""Metrics with the presentation_audit_20261117 definitions: max-gm linear
extrapolation at |Vd|=50 mV (x_int - Vd/2), constant current 100 nA/um, SS
best 4-point fit in 1e-10..1e-7 A, DIBL from constant-current shift."""
import re,json,numpy as np
from pathlib import Path
here=Path(__file__).resolve().parent
def parse(p):
 cols={};rows=[]
 for line in p.read_text().splitlines():
  if m:=re.match(r'#\s+(\d+)\s+(.+)\s+\[([^]]+)\]',line):cols[m[2].strip()]=int(m[1])-1
  elif line and not line.startswith('#'):rows.append(list(map(float,line.split())))
 r=np.array(rows);ix=lambda s:[j for k,j in cols.items() if k.startswith(s)]
 g=ix('Vapp(Gate[')[0];di=ix('I(Drain[')[0];si=ix('I(Source[')[0];bi=ix('I(Body[')[0];gi=ix('I(Gate[')
 kcl=float(np.max(abs(r[:,di]+r[:,si]+r[:,bi]+r[:,gi].sum(axis=1))))
 return abs(r[:,g]),abs(r[:,di]),kcl,r[:,g],r[:,di]
def cross(x,y,t):
 for i in range(len(x)-1):
  if y[i]<=t<=y[i+1]:return float(x[i]+(t-y[i])/(y[i+1]-y[i])*(x[i+1]-x[i]))
out={}
for dev in ('nmos','pmos'):
 s=1 if dev=='nmos' else -1
 x,lo,k1,_,_=parse(here/dev/'transfer_vd0p05/transfer.dat');_,hi,k2,_,_=parse(here/dev/'transfer_vd3p3/transfer.dat')
 _,neg,k3,gneg,dneg=parse(here/dev/'transfer_neg_vd0p05/transfer.dat')
 gm=np.diff(lo)/np.diff(x);j=int(np.argmax(gm));vt=s*(x[j]-lo[j]/gm[j]-.025)
 cl=cross(x,lo,1e-7);ch=cross(x,hi,1e-7)
 w=[]
 for k in range(len(x)-3):
  yy=lo[k:k+4]
  if min(yy)>=1e-10 and max(yy)<=1e-7 and yy[-1]/yy[0]>=10**.5:
   sl=np.polyfit(x[k:k+4],np.log10(yy),1)[0]
   if sl>0:w.append(1000/sl)
 out[dev]={'Vt_maxgm_V':round(vt,3),'Vt_target_V':s*.75,'Vt_cc_100nA_V':None if cl is None else round(s*cl,3),
  'Vt_cc_100nA_Vd3p3_V':None if ch is None else round(s*ch,3),'DIBL_mV_per_V':None if cl is None or ch is None else round((cl-ch)/3.25*1000,1),
  'SS_mV_dec':round(min(w),1) if w else None,'Id_Vg0_Vd0p05_A':float(lo[0]),'Id_Vg0_Vd3p3_A':float(hi[0]),
  'Ion_Vg3p3_Vd3p3_A':float(hi[-1]),'Ion_Vg3p3_Vd0p05_A':float(lo[-1]),'gm_max_A_per_V':float(gm[j]),
  'Id_at_Vg_-1V_reverse_A':float(abs(dneg[-1])),'Vg_reverse_last_V':float(gneg[-1]),'max_KCL_A':max(k1,k2,k3)}
 # dump a few points for inspection
 out[dev]['Id_lowVd_samples']={f'{v:.2f}':float(i) for v,i in zip(x[::6],lo[::6])}
old=json.loads(Path('/home/chris/Open_TCAD/suprem4/experiments/presentation_audit_20261117/electrical_metrics.json').read_text())
for dev in out:
 o=old[dev];out[dev]['reference_20261117']={'Vt_maxgm_V':o['Vt_maxgm_V'],'Vt_cc_100nA_V':o['Vt_constant_100nA_V'],'DIBL':o['DIBL_constant_100nA_mV_per_V'],'SS':(o['SS_best_four_point_low_drain'] or {}).get('SS_mV_per_dec'),'Ioff':o['low_drain']['Ioff_A_per_um'],'Ion':o['high_drain']['Ion_A_per_um']}
(here/'metrics.json').write_text(json.dumps(out,indent=1)+'\n');print(json.dumps(out,indent=1))
