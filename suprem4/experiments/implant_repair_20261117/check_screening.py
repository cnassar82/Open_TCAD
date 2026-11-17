from pathlib import Path
import os,subprocess,json,math,statistics,numpy as np,sys
root=Path('/home/chris/Open_TCAD');binary=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(exist_ok=False)
N=statistics.NormalDist();dose=1e13;oxrp=.0623;oxsigma=.0238;sirp=.0663;sisigma=.0284
env=os.environ.copy()
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),('IMPFILE','sup4gs.imp')]:env[k]=str(binary.parent.parent/'data'/n)
records=[]
for thick in (0,.025,.05,.1,.15):
 d=out/str(thick);d.mkdir();deck='''option quiet
mode one.dim
line x loc=0 spac=0.002 tag=top
line x loc=1 spac=0.01 tag=bottom
region silicon xlo=top xhi=bottom
boundary exposed xlo=top xhi=top
boundary backside xlo=bottom xhi=bottom
init boron conc=1e15 orient=100
'''
 if thick:deck+=f'deposit oxide thick={thick} divisions={round(thick/.001)}\n'
 deck+='implant phosphorus dose=1e13 energy=50 gauss\nstructure outfile=after.str\nquit\n';(d/'input.in').write_text(deck)
 with (d/'run.log').open('w') as log:r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60)
 assert r.returncode==0
 coord={};species=[];values=[]
 for line in (d/'after.str').read_text().splitlines():
  a=line.split()
  if a[0]=='c':coord[int(a[1])]=float(a[2])
  elif a[0]=='s':species=list(map(int,a[2:]))
  elif a[0]=='n' and int(a[2])==3:values.append((coord[int(a[1])+1],float(a[3+species.index(3)])))
 values.sort();z=np.array([v[0] for v in values]);c=np.array([v[1] for v in values])
 retained=1 if thick==0 else (1-N.cdf((thick-oxrp)/oxsigma))/N.cdf(oxrp/oxsigma)
 entered=1-retained;entry=sirp+sisigma*N.inv_cdf(N.cdf(-sirp/sisigma)+entered*N.cdf(sirp/sisigma))
 expected=dose*np.exp(-.5*((z+entry-sirp)/sisigma)**2)/(sisigma*1e-4*math.sqrt(2*math.pi)*N.cdf(sirp/sisigma))
 mask=expected>expected.max()*1e-3;err=float(np.max(abs(c[mask]/expected[mask]-1)))
 rec=dict(oxide_um=thick,expected_retained_fraction=retained,actual_Si_sheet_cm2=float(np.trapezoid(c,z)*1e-4),expected_Si_sheet_cm2=dose*retained,max_relative_profile_error=err);print(rec,flush=True);records.append(rec)
 np.savetxt(d/'profile.csv',np.column_stack((z,c,expected)),delimiter=',',header='Si_depth_um,actual_P_cm3,expected_P_cm3',comments='')
(out/'results.json').write_text(json.dumps(records,indent=2)+'\n')
