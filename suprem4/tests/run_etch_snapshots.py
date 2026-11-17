#!/usr/bin/env python3
"""Stress a captured strip/spacer operation with thickness and origin changes.

This tests local geometric robustness; it does not establish full-flow convergence.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_oxidation_perturbations import audit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from str_compare import parse_str,summarize
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--input',type=Path,required=True)
p.add_argument('--depth',type=float,required=True)
p.add_argument('--output',type=Path)
a=p.parse_args();binary=a.binary.resolve();data=binary.parent.parent/'data'
source=a.input.resolve();audit(source)
work=a.output or Path(tempfile.mkdtemp(prefix='suprem4_etch_snapshots_'))
work.mkdir(parents=True,exist_ok=True);print(f'workdir: {work}',flush=True)
env=os.environ.copy()
for key,name in [('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc')]:env[key]=str(data/name)
def run(case):
 factor,shift=case;folder=work/f'depth{factor:g}_origin{shift:g}';folder.mkdir()
 lines=[]
 for l in source.read_text().splitlines():
  s=l.split()
  if s and s[0]=='c':s[2]=f'{float(s[2])+shift:.17g}';l=' '.join(s)
  lines.append(l)
 (folder/'input.str').write_text('\n'.join(lines)+'\n')
 text=f'structure infile=input.str\netch oxide dry thick={a.depth*factor:.17g}\nstructure outfile=result.str\nquit\n'
 r=subprocess.run([str(binary)],input=text,text=True,cwd=folder,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
 (folder/'suprem4.log').write_text(r.stdout)
 result={'depth_factor':factor,'origin_shift_um':shift,'returncode':r.returncode,'passed':False}
 try:
  if r.returncode or 'panic' in r.stdout.lower():raise ValueError('simulator failed')
  result['audit']=audit(folder/'result.str');result['summary']=summarize(parse_str(folder/'result.str'));result['passed']=True
 except (ValueError,OSError,KeyError) as e:result['error']=str(e)
 print(('PASS' if result['passed'] else 'FAIL')+f' depth={factor:g} origin={shift:g}: {result.get("error","")}',flush=True)
 return result
cases=[(f,x) for f in (.99,1,1.01) for x in (-.01,0,.01)]
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,cases))
(work/'results.json').write_text(json.dumps({'input':str(source),'binary':str(binary),'depth_um':a.depth,'results':results},indent=2)+'\n')
raise SystemExit(not all(r['passed'] for r in results))
