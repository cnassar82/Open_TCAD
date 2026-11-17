"""Assemble the dated audit's reproducible comparison and provenance."""
from pathlib import Path
import hashlib,json,csv,re
root=Path('/home/chris/Open_TCAD');out=Path(__file__).resolve().parent
m=json.loads((out/'profile_metrics.json').read_text());e=json.loads((out/'electrical_metrics.json').read_text())
process=json.loads((root/'decks/subcmos2014/process.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def thickness(dev,stage,x,material):
 return sum(b-a for a,b in m[dev][stage]['film_spans_um'][str(float(x))][str(material)])*1000
rows=[]
def row(dev,name,stage,target,actual,unit,page,status='deviation',note=''):
 rows.append(dict(device=dev,quantity=name,checkpoint=stage,reference=target,actual=actual,unit=unit,pdf_pages=page,status=status,relative_error_percent=None if target is None or target==0 else 100*(actual-target)/abs(target),note=note))
for dev in m:
 row(dev,'initial substrate B','initial',1e15,m[dev]['initial']['cuts']['channel']['samples']['0']['B'],'cm^-3','11,18','matches within storage precision')
 row(dev,'pad oxide','pad',50,thickness(dev,'pad',6,1),'nm','23,24')
 row(dev,'well mask oxide after nominal 25nm strip','well_mask',25 if dev=='nmos' else 275,thickness(dev,'well_mask',6,1),'nm','41,43,48')
 row(dev,'LOCOS nitride','locos_mask',150,thickness(dev,'locos_mask',6,2),'nm','60','matches geometric deposition')
 row(dev,'field oxide after pad strip','locos',600,thickness(dev,'locos',2.5,1),'nm','69')
 row(dev,'Kooi channel oxide','threshold',100,thickness(dev,'threshold',6,1),'nm','71,73')
 row(dev,'final gate oxide','metal1',15,thickness(dev,'metal1',6,1),'nm','11,80–82')
 row(dev,'final thermal field oxide (before ILD)','annealed',500,thickness(dev,'annealed',2.5,1),'nm','97,133','deviation; reference has inconsistent labels','p11/p85 use 600nm; p78/p97/p133 use 500nm')
 row(dev,'added ILD at field','metal1',500,thickness(dev,'metal1',2.5,1)-thickness(dev,'annealed',2.5,1),'nm','117,118','matches geometric deposition','p17 summary says 400nm')
 row(dev,'final poly thickness','metal1',500,m[dev]['metal1']['poly_channel']['thickness_um']*1000,'nm','85','matches geometric deposition')
 row(dev,'low-drain maximum-gm threshold','metal1',e[dev]['Vt_target_V'],e[dev]['Vt_maxgm_V'],'V','11','deviation','50mV drain and gate sampling; presentation gives no extraction convention')
 row(dev,'3.3V off-current','metal1',None,e[dev]['high_drain']['Ioff_A_per_um'],'A/um','11','characterization; no current target in PDF')
 row(dev,'3.3V on-current','metal1',None,e[dev]['high_drain']['Ion_A_per_um'],'A/um','11','characterization; no current target in PDF')
row('pmos','N-well PN depth after drive','wells',2.5,m['pmos']['wells']['cuts']['channel']['chemical_junction_depths_um'][-1],'um','11')
row('nmos','P-well B=2x background marker','wells',3.,m['nmos']['wells']['cuts']['channel']['B_2e15_crossing_depths_um'][-1],'um','11','not a like-for-like junction comparison','P-well in p-substrate has no bottom PN junction; marker is an explicit diagnostic definition')
poly={}
ns={'__file__':str(out/'extract_profiles.py')};exec((out/'extract_profiles.py').read_text().split('summary={};all_profiles={}')[0],ns)
import numpy as np
for dev in m:
 y,v=ns['cut'](ns['read'](out/'poly_implant'/dev/'after.str'),6,4)
 sheet=float(np.trapezoid(v[:,1],y)*1e-4)
 poly[dev]={'gate_center_sheet_cm2':sheet,'fraction_of_requested_2e16':sheet/2e16}
 row(dev,'as-implanted gate-center poly P sheet','extra pre-activation control',2e16,sheet,'cm^-2','87','fails spatial resolution control')
files=[root/'rit_process/SubCmos2014_0p5um.pdf',root/'suprem4/build/native64_rxq2pnrb/src/suprem4',root/'suprem4/experiments/rit_electrical_20261116/genius_runtime/bin/genius.LINUX']
files += [root/f'decks/subcmos2014/{d}.in' for d in m]
files += [root/'decks/subcmos2014/process.json']
files += [root/f'decks/subcmos2014/results/directional_final_fresh_20261117/{d}/metal1.str' for d in m]
files += [root/f'suprem4/experiments/directional_etch_20261117/electrical/{d}/device.tif' for d in m]
checks={}
for dev in m:
 checks[dev]=[{'name':r['name'],'pdf_page':r['pdf_page'],'deck_total_minutes':sum(s[0] for s in r['segments'])} for r in process['devices'][dev]['recipes']]
inputs={}
for dev in m:
 inputs[dev]=re.findall(r'^implant .+$',(root/f'decks/subcmos2014/{dev}.in').read_text(),re.M)
result={'date':'2026-11-17','presentation_agreement_verified':False,'poly_implant_resolution_qualified':False,'scope':'current corrected full process cuts through first metal; not full two-metal presentation process','comparison':rows,'poly_implant_control':poly,'electrical':e,'input_implant_commands':inputs,'recipes':checks,'profile_csv_count':len(list((out/'profiles').glob('*.csv'))),'provenance_sha256':{str(p):sha(p) for p in files}}
assert result['profile_csv_count']==224
(out/'presentation_comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
with (out/'presentation_comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps({'comparisons':len(rows),'profiles':result['profile_csv_count'],'poly':poly},indent=2))
