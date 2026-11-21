#!/usr/bin/env python3
"""Dopant transport through a thin gate oxide must not depend on how many
element layers the oxide was created with.

Gate stack of the RIT process: 16.9 nm oxide, 0.5 um poly, P 2e16 cm^-2 at
50 keV, 117-minute N2 activation. The oxide is (a) deposited as ONE element
layer, (b) grown thermally (one to two irregular layers), and (c) deposited
in 16 layers. The reference (r) builds the same stack from line/region
statements: a rectangular grid of right triangles, whose control volumes
are exact Voronoi cells under any area formula, with thin-film refinement
off. Thin-film conformity must give (a) and (b) at least 8 dielectric
edges between silicon and poly, and (a), (b), (c) must hold P in silicon
within 2x of (r) and match its surface polarity. (Under the legacy |cot|
control volumes, flat oxide triangles inflated node volumes ~3x and the
layered cases fell ~100x below (r).) A control runs (a) with conformity
disabled and must show the one-layer defect, so the test cannot pass
vacuously.
"""
import argparse, os, subprocess, sys
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'decks/subcmos2014'));sys.path.insert(0,str(ROOT/'suprem4/tools'))
from validate import read
import str_compare

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
data=binary.parent.parent/'data'
base={**{k:v for k,v in os.environ.items() if not k.startswith('SUPREM4_')}}
for k,n in [('KEYFILE','suprem.uk'),('SUP4KEYFILE','suprem.uk'),('SUP4MODELRC','modelrc'),
            ('IMPFILE','sup4gs.imp'),('SUP4IMPDATA','sup4gs.imp')]:base[k]=str(data/n)

HEAD='''option quiet
method fermi vertical grid.oxide=0.04
line x loc=0 spac=0.5 tag=left
line x loc=1 spac=0.5 tag=right
line y loc=0 spac=0.005 tag=top
line y loc=0.3 spac=0.02
line y loc=2 spac=0.2 tag=bottom
region silicon xlo=left xhi=right ylo=top yhi=bottom
boundary exposed xlo=left xhi=right ylo=top yhi=top
boundary backside xlo=left xhi=right ylo=bottom yhi=bottom
initialize boron conc=1e16 orient=100
'''
GROW='''diffuse time=12 temp=800 nitrogen
diffuse time=15 temp=800 nitrogen continue
diffuse time=5 temp=825 dryo2
diffuse time=5 temp=875 dryo2 continue
diffuse time=45 temp=900 dryo2 continue
diffuse time=10 temp=875 nitrogen continue
diffuse time=10 temp=825 nitrogen continue
diffuse time=15 temp=800 nitrogen continue
'''
TAIL='''deposit poly thick=0.5 divisions=10 space=0.01
implant phosphorus dose=2e+16 energy=50 pearson
'''+''.join(f'diffuse time={t} temp={T} nitrogen{" continue" if i else ""}\n' for i,(t,T) in enumerate(
    [(12,800),(15,800),(5,825),(5,875),(5,925),(5,975),(15,1000),(10,975),(10,925),(10,875),(10,825),(15,800)]))+'''structure outfile=final.str
quit
'''
RECT=HEAD.replace('line y loc=0 spac=0.005 tag=top',
    'line y loc=-0.5169 spac=0.01 tag=ptop\nline y loc=-0.0169 spac=0.002 tag=pbot\n'
    'line y loc=0 spac=0.002 tag=top').replace('line x loc=1 spac=0.5','line x loc=1 spac=0.1').replace(
    'line x loc=0 spac=0.5','line x loc=0 spac=0.1').replace('boundary exposed xlo=left xhi=right ylo=top yhi=top',
    'region poly xlo=left xhi=right ylo=ptop yhi=pbot\nregion oxide xlo=left xhi=right ylo=pbot yhi=top\n'
    'boundary exposed xlo=left xhi=right ylo=ptop yhi=ptop')+TAIL.split('pearson\n',1)[0].split('\n',1)[1]+'pearson\n'+TAIL.split('pearson\n',1)[1]
cases={'rect_reference':(RECT,{'SUPREM4_THIN_FILM_LAYERS':'0'}),
       'one_layer':(HEAD+'deposit oxide thick=0.0169 divisions=1\n'+TAIL,{}),
       'grown':(HEAD+GROW+TAIL,{}),
       'reference_16':(HEAD+'deposit oxide thick=0.0169 divisions=16\n'+TAIL,{'SUPREM4_THIN_FILM_LAYERS':'0'}),
       'control_disabled':(HEAD+'deposit oxide thick=0.0169 divisions=1\n'+TAIL,{'SUPREM4_THIN_FILM_LAYERS':'0'})}

def run(name):
    deck,extra=cases[name];d=out/name;d.mkdir()
    (d/'input.in').write_text(deck)
    r=subprocess.run([str(binary)],input=deck,text=True,cwd=d,env={**base,**extra},
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=3600)
    (d/'run.log').write_text(r.stdout)
    if r.returncode or 'panic' in r.stdout:raise SystemExit(f'{name}: simulator failed ({r.returncode})')
    return name

def inventory(path):
    coords,tris,nodes=read(path);xs=[c[0] for c in coords.values()];W=(max(xs)-min(xs))*1e-4;tot=0
    for m,vs in tris:
        if m!=3:continue
        (ax,ay),(bx,by),(cx,cy)=(coords[v] for v in vs)
        tot+=abs((bx-ax)*(cy-ay)-(by-ay)*(cx-ax))/2e8*sum(nodes[v,3].get(3,0.) for v in vs)/3
    return tot/W

def surface_net(path):
    """Minimum (B - P) over silicon nodes at the Si/oxide interface."""
    coords,tris,nodes=read(path);mats=defaultdict(set)
    for m,vs in tris:
        for v in vs:mats[v].add(m)
    return min(nodes[v,3].get(5,0.)-nodes[v,3].get(3,0.) for v in mats if 3 in mats[v] and 1 in mats[v])

def min_hops(path):
    d=str_compare.parse_str(path);R=d['regions'];mats=defaultdict(set);adj=defaultdict(set)
    for t in d['triangles']:
        for v in t['verts']:mats[v].add(R[t['region']])
        if R[t['region']]==1:
            v=t['verts']
            for i in range(3):adj[v[i]].add(v[(i+1)%3]);adj[v[(i+1)%3]].add(v[i])
    start=[v for v in adj if 3 in mats[v]];dist={v:0 for v in start};q=deque(start)
    while q:
        u=q.popleft()
        for w in adj[u]:
            if w not in dist:dist[w]=dist[u]+1;q.append(w)
    return min(dist[v] for v in adj if 4 in mats[v] and v in dist)

with ThreadPoolExecutor(4) as pool:list(pool.map(run,cases))
res={n:{'P_in_Si_cm2':inventory(out/n/'final.str'),'min_hops':min_hops(out/n/'final.str'),
        'min_surface_B_minus_P':surface_net(out/n/'final.str')} for n in cases}
R=res['rect_reference'];ref=R['P_in_Si_cm2'];failed=False
for n,r in res.items():print(n,r,flush=True)
for n in ('one_layer','grown','reference_16'):
    r=res[n];ok=(r['min_hops']>=8 and ref/2<=r['P_in_Si_cm2']<=2*ref and
                 (r['min_surface_B_minus_P']>0)==(R['min_surface_B_minus_P']>0))
    failed|=not ok;print(('PASS ' if ok else 'FAIL ')+n,flush=True)
if res['control_disabled']['P_in_Si_cm2']<=5*ref:
    failed=True;print('FAIL control: disabled conformity did not reproduce the defect',flush=True)
else:print('PASS control: disabled conformity reproduces the one-layer defect',flush=True)
sys.exit(1 if failed else 0)
