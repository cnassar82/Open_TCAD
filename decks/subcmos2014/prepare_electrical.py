#!/usr/bin/env python3
"""Export the four-contact RIT cuts for an explicitly assumed Genius DD model.

Preserves every mesh vertex/triangle and exports active dopants. Connected n+
poly and gate metal become one ideal gate conductor; poly depletion is omitted.
No electrical parameter is fitted to the PDF threshold target.
"""
import argparse,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'suprem_viewer'))
from make_genius_contacts import parse_suprem,region_extents,ContactSpec,write_tif,build_edges
from validate import read,gate_dielectric
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--process-results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--allow-invalid-gate-diagnostic',action='store_true',help='Explicitly allow a missing-dielectric structure for failure diagnosis only')
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);manifest={}
for device in ('nmos','pmos'):
 d=a.output/device;d.mkdir();original=a.process_results/device/'metal1.str'
 try:gate_check={'passed':True,'probes':gate_dielectric(read(original))}
 except ValueError as error:
  gate_check={'passed':False,'error':str(error)}
  if not a.allow_invalid_gate_diagnostic:p.error(f'{device}: {error}; invalid electrical input')
 nodes,tris,regions,codes,rows=parse_suprem(original,auto_name_electrodes=False);bounds=region_extents(nodes,tris)
 metal=sorted([r for r in regions.values() if r.mat_id==6],key=lambda r:bounds[r.idx][0]);assert len(metal)==4
 names=['Body','Source','Gate','Drain'];contacts=[ContactSpec(name,'gate' if name=='Gate' else 'ohmic','region',region=r.idx,work_function=4.17) for name,r in zip(names,metal)]
 poly=[r.idx for r in regions.values() if r.mat_id==4];assert len(poly)==1
 _,_,owners=build_edges(tris);gate=contacts[2].region
 shared=sum(owner_set=={poly[0],gate} for owner_set in owners.values());assert shared>0,'poly and gate metal do not make contact'
 for tri in tris:
  if tri.region==poly[0]:tri.region=gate
 del regions[poly[0]]
 # Medici region IDs are dense one-based indices, unlike arbitrary mesh tags.
 contact_bounds={c.name:bounds[c.region] for c in contacts}
 original_contact_ids={c.name:c.region for c in contacts}
 renumber={old:i+1 for i,old in enumerate(sorted(regions))}
 for tri in tris:tri.region=renumber[tri.region]
 for region in regions.values():region.idx=renumber[region.idx]
 regions={region.idx:region for region in regions.values()}
 for contact in contacts:contact.region=renumber[contact.region]
 counts=write_tif(d/'device.tif',nodes,tris,regions,codes,rows,1e-6,contacts)
 silicon=[r.name for r in regions.values() if r.mat_id==3];assert len(silicon)==1
 deck='GLOBAL TExternal=300 DopingScale=1e18 Z.Width=1.0\nIMPORT TIFFile=device.tif\n'
 for c in contacts:deck+=f'CONTACT Type={c.genius_type} ID={c.name} '+('WorkFunction=4.17 ' if c.kind=='gate' else '')+'Res=0 Cap=0 Ind=0\n'
 # Override Genius's nonzero default explicitly; this is an unfitted baseline.
 _,_,electrical_owners=build_edges(tris)
 interface_pairs={tuple(sorted(pair)) for pair in electrical_owners.values() if len(pair)==2 and {regions[r].mat_id for r in pair}=={1,3}}
 for pair in sorted(interface_pairs):
  label='_to_'.join(sorted(regions[r].name for r in pair))
  deck+=f'BOUNDARY Type=InsulatorInterface ID={label} QF=0\n'
 deck+=f'PMI Region={silicon[0]} Type=Mobility Model=Lucent\nMODEL Region={silicon[0]} Fermi=true\n'
 deck+='METHOD Type=DDML1 NS=Basic LS=MUMPS MaxIteration=80 Damping=Potential\nSOLVE Type=EQUILIBRIUM\n'
 (d/'equilibrium.in').write_text(deck+'EXPORT CGNSFile=equilibrium.cgns VTKFile=equilibrium.vtu\nEND\n')
 sign=1 if device=='nmos' else -1
 # A drain ramp reaches the low drain bias before the gate transfer scan.
 transfer=deck+f'SOLVE Type=DC Vscan=Drain Vstart=0 Vstep={sign*.01:g} Vstop={sign*.05:g} out.prefix=drain_prep\n'
 transfer+=f'SOLVE Type=DC Vscan=Gate Vstart=0 Vstep={sign*.05:g} Vstop={sign*3.3:g} out.prefix=transfer\n'
 transfer+='EXPORT CGNSFile=transfer_final.cgns VTKFile=transfer_final.vtu\nEND\n';(d/'transfer.in').write_text(transfer)
 manifest[device]={'gate_oxide_material_validation':gate_check,'diagnostic_only':not gate_check['passed'],'input_structure':str(original.resolve()),'input_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'tif_sha256':hashlib.sha256((d/'device.tif').read_bytes()).hexdigest(),'counts':counts,'contacts':{c.name:{'region':c.region,'original_region':original_contact_ids[c.name],'bounds_um':contact_bounds[c.name]} for c in contacts},'merged_poly_region':poly[0],'poly_metal_shared_edges':shared,'gate_work_function_eV':4.17,'width_um':1.,'temperature_K':300.,'fixed_interface_charge_cm2':0.,'drain_bias_V':sign*.05,'active_dopants':True,'assumptions':['Ideal n+ gate; 4.17 eV consistent with installed Si affinity model','Active B/P profiles; fully ionized fixed dopants','Lucent mobility; Fermi statistics','No fitted interface charge or interface traps','One micron device width for current normalization']}
(a.output/'export_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(a.output.resolve())
