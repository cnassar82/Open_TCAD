#!/usr/bin/env python3
"""Report diagnostic Genius transfer sweeps against PDF nominal targets.

No threshold is inferred from invalid missing-dielectric structures.
"""
import argparse,json,re,math,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validate import read,intervals,concentration,gate_dielectric
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--process-results',type=Path,required=True);p.add_argument('--electrical-results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);report={'pdf':'rit_process/SubCmos2014_0p5um.pdf','target_page':11,'devices':{}};fig,axes=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
for ax,dev in zip(axes,('nmos','pmos')):
 d=a.electrical_results/dev;cols={};rows=[]
 for line in (d/'transfer.dat').read_text().splitlines():
  if m:=re.match(r'#\s+(\d+)\s+(.+)\s+\[([^]]+)\]',line):cols[m[2].strip()]=int(m[1])-1
  elif line and not line.startswith('#'):rows.append(list(map(float,line.split())))
 def column(prefix):return [v for k,v in cols.items() if k.startswith(prefix)]
 gv=column('Vapp(Gate[');dv=column('Vapp(Drain[')[0];di=column('I(Drain[')[0];si=column('I(Source[')[0];bi=column('I(Body[')[0];gi=column('I(Gate[')
 assert len(rows)==67 and len(gv)>0
 sign=1 if dev=='nmos' else -1
 assert all(abs(row[dv]-sign*.05)<1e-12 and all(abs(row[j]-row[gv[0]])<1e-12 for j in gv) for row in rows)
 assert abs(rows[-1][gv[0]]-sign*3.3)<1e-12
 assert all(all(math.isfinite(v) for v in row) for row in rows)
 execution=json.loads((d/'transfer_execution.json').read_text());log=(d/'transfer_stdout.log').read_text();assert execution['returncode']==0 and 'DIVERGED' not in log and 'ERROR' not in log
 vg=[r[gv[0]] for r in rows];idrain=[r[di] for r in rows]
 mesh=read(a.process_results/dev/'metal1.str');gateox=read(a.process_results/dev/'gate_oxide.str');well=read(a.process_results/dev/'wells.str');annealed=read(a.process_results/dev/'annealed.str')
 oxide_before=intervals(gateox,6,1);tox=oxide_before[-1][1]-oxide_before[-1][0]
 field=intervals(annealed,2.5,1);field_thick=sum(b-c for c,b in field)
 surface=intervals(well,6,3)[0][0];well05=concentration(well,6,surface+.5)['net_donor_cm3']
 try:occupancy={'passed':True,'probes':gate_dielectric(mesh)}
 except ValueError as error:occupancy={'passed':False,'error':str(error)}
 r={'execution':execution,'converged_gate_bias_points':len(rows),'drain_bias_V':sign*.05,'width_um':1,'max_terminal_kcl_error_A':max(abs(row[di]+row[si]+row[bi]+sum(row[j] for j in gi)) for row in rows),'threshold_target_V':sign*.75,'extracted_threshold_V':None,'threshold_status':'Not meaningful: final gate dielectric missing','gate_oxide_before_contact_nm':tox*1000,'gate_oxide_target_nm':15,'pre_contact_oxide_error_percent':100*(tox/.015-1),'final_gate_oxide':occupancy,'field_oxide_probe_x_um':2.5,'field_oxide_pre_ILD_nm':field_thick*1000,'well_net_doping_probe_depth_um':.5,'well_net_doping_cm3':well05,'well_nominal_target_cm3':1e16 if dev=='nmos' else 3e16,'Id_at_Vg0_A':idrain[0],'Id_at_max_gate_drive_A':idrain[-1],'max_abs_Id_A':max(map(abs,idrain)),'electrically_validated':False}
 report['devices'][dev]=r
 (a.output/(dev+'_transfer.csv')).write_text('Vgs_V,Vds_V,Id_A,width_um\n'+''.join(f'{v:.12g},{sign*.05:.12g},{i:.16g},1\n' for v,i in zip(vg,idrain)))
 if dev=='nmos':ax.plot([abs(v) for v in vg],[i*1e15 for i in idrain],'.-',lw=.8,ms=3,color='#236b8e');ax.set_ylabel('Drain current (fA / um)');ax.set_ylim(-1,1);ax.set_title('NMOS: numerical-floor current')
 else:ax.plot([abs(v) for v in vg],[abs(i)*1e6 for i in idrain],lw=2,color='#a34173');ax.set_ylabel('|Drain current| (uA / um)');ax.set_ylim(0,1.05);ax.set_title('PMOS: essentially no gate modulation')
 ax.axvline(.75,ls='--',lw=1,color='black',label='PDF |VT| target: 0.75 V');ax.set_xlabel('Gate drive |Vgs| (V)');ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper left');ax.set_xlim(0,3.3)
fig.suptitle('Diagnostic sweeps: final gate oxide is missing\nThese curves do not validate transistor thresholds',fontsize=12)
fig.savefig(a.output/'diagnostic_transfer.png',dpi=170);fig.savefig(a.output/'diagnostic_transfer.svg')
(a.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
