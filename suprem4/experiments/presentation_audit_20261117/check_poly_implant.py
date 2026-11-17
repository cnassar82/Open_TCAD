from pathlib import Path
import os,subprocess,json,sys
root=Path('/home/chris/Open_TCAD');out=Path(__file__).resolve().parent/'poly_implant';out.mkdir(exist_ok=True)
binary=root/'suprem4/build/native64_rxq2pnrb/src/suprem4';data=binary.parent.parent/'data'
env={**os.environ,'SUPREM4_PROCESS_SURFACE':'1','SUPREM4_DEPOSITION_GEOMETRY':'domain','SUPREM4_MATERIAL_DOMAINS':'oxide','SUPREM4_MESH_PATCH_REMOVAL':'1'}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(data/n)
for dev in ('nmos','pmos'):
 d=out/dev;d.mkdir(exist_ok=False);source=root/'decks/subcmos2014/results/directional_final_fresh_20261117'/dev/'gate_oxide.str'
 deck=f'option quiet\nmethod fermi vertical grid.oxide=0.04\nstructure infile={source}\ndeposit poly thick=0.5 divisions=10\nstructure outfile=before.str\nimplant phosphorus dose=2e16 energy=50 pearson\nstructure outfile=after.str\nquit\n';(d/'input.in').write_text(deck)
 with (d/'run.log').open('w') as log:r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
 print(dev,r.returncode,flush=True);assert r.returncode==0
