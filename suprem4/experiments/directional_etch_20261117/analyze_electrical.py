import re,json,math,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent/'electrical'
manifest=json.loads((base/'export_manifest.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(9,3.7),layout='constrained');summary={}
for dev,color in [('nmos','#236b8e'),('pmos','#a34173')]:
 d=base/dev;cols={};rows=[]
 for line in (d/'transfer.dat').read_text().splitlines():
  if m:=re.match(r'#\s+(\d+)\s+(.+)\s+\[([^]]+)\]',line):cols[m[2].strip()]=int(m[1])-1
  elif line and not line.startswith('#'):rows.append(list(map(float,line.split())))
 indices=lambda p:[j for k,j in cols.items() if k.startswith(p)]
 gates=indices('Vapp(Gate[');drain=indices('Vapp(Drain[')[0];di=indices('I(Drain[')[0]
 si=indices('I(Source[')[0];bi=indices('I(Body[')[0];gi=indices('I(Gate[')
 sign=1 if dev=='nmos' else -1
 assert len(rows)==67 and gates and all(all(math.isfinite(v) for v in row) for row in rows)
 assert all(abs(row[drain]-sign*.05)<1e-12 and all(abs(row[g]-sign*i*.05)<1e-12 for g in gates) for i,row in enumerate(rows))
 log=(d/'transfer_stdout.log').read_text();execution=json.loads((d/'transfer_execution.json').read_text())
 assert execution['returncode']==0 and not re.search(r'ERROR|DIVERGED',log)
 original=Path(manifest[dev]['input_structure']);assert hashlib.sha256(original.read_bytes()).hexdigest()==manifest[dev]['input_sha256']
 assert manifest[dev]['gate_oxide_material_validation']['passed'] and not manifest[dev]['diagnostic_only']
 vg=[abs(row[gates[0]]) for row in rows];current=[abs(row[di]) for row in rows]
 kcl=max(abs(row[di]+row[si]+row[bi]+sum(row[g] for g in gi)) for row in rows)
 summary[dev]={'execution':execution,'converged_gate_bias_points':len(rows),'drain_bias_V':sign*.05,'width_um':1,'gate_oxide_material_validation':manifest[dev]['gate_oxide_material_validation'],'input_structure_sha256':manifest[dev]['input_sha256'],'Id_magnitude_at_Vg0_A':current[0],'Id_magnitude_at_3p3V_drive_A':current[-1],'on_off_ratio':current[-1]/current[0],'max_terminal_kcl_error_A':kcl,'threshold_extracted':False,'PDF_threshold_calibrated':False,'model_assumptions':manifest[dev]['assumptions']}
 (d/'transfer.csv').write_text('Vgs_V,Vds_V,Id_A,width_um\n'+''.join(f'{row[gates[0]]:.12g},{row[drain]:.12g},{row[di]:.16g},1\n' for row in rows))
 for ax in axes:ax.plot(vg,[i*1e6 for i in current],label=dev.upper(),color=color,lw=2);ax.set_xlabel('Gate drive |Vgs| (V)');ax.set_ylabel('|Id| (uA / um)');ax.grid(alpha=.2);ax.legend();ax.set_xlim(0,3.3)
 axes[1].set_yscale('log')
axes[0].set_title('Linear current');axes[1].set_title('Log current')
fig.suptitle('Corrected final RIT cuts: |Vds| = 50 mV, 300 K\nUnfitted DD baseline; PDF threshold calibration not established',fontsize=11)
fig.savefig(base/'corrected_transfer.png',dpi=180);fig.savefig(base/'corrected_transfer.svg')
(base/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
