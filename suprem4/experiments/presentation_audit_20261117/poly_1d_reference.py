from pathlib import Path
import os,subprocess,json,numpy as np
root=Path('/home/chris/Open_TCAD');out=Path(__file__).resolve().parent/'poly_1d_reference';out.mkdir(exist_ok=True)
binary=root/'suprem4/build/native64_rxq2pnrb/src/suprem4';env=os.environ.copy()
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:env[k]=str(binary.parent.parent/'data'/n)
results=[]
for div in (10,50,250):
 d=out/str(div);d.mkdir(exist_ok=False)
 deck=f'''option quiet
mode one.dim
line x loc=0 spac=0.01 tag=top
line x loc=2 spac=0.1 tag=bottom
region silicon xlo=top xhi=bottom
boundary exposed xlo=top xhi=top
boundary backside xlo=bottom xhi=bottom
init boron conc=1e15 orient=100
deposit oxide thick=0.015 divisions=15
deposit poly thick=0.5 divisions={div}
implant phosphorus dose=2e16 energy=50 pearson
structure outfile=after.str
quit
''';(d/'input.in').write_text(deck)
 with (d/'run.log').open('w') as log:r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60)
 assert r.returncode==0
 coords={};species=[];values=[]
 for line in (d/'after.str').read_text().splitlines():
  v=line.split()
  if v[0]=='c':coords[int(v[1])]=float(v[2])
  elif v[0]=='s':species=list(map(int,v[2:]))
  elif v[0]=='n' and int(v[2])==4:values.append((coords[int(v[1])+1],float(v[3+species.index(3)])))
 values.sort();y=np.array([p[0] for p in values]);c=np.array([p[1] for p in values]);depth=y-y[0]
 record={'divisions':div,'returncode':r.returncode,'P_sheet_cm2':float(np.trapezoid(c,y)*1e-4),'peak_P_cm3':float(c.max()),'peak_depth_um':float(depth[c.argmax()])}
 np.savetxt(d/'profile.csv',np.column_stack((depth,c)),delimiter=',',header='depth_um,P_cm3',comments='')
 results.append(record);print(record,flush=True)
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
