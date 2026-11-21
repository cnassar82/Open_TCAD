#!/usr/bin/env python3
"""Archive isolated source differences; not an instruction to apply them."""
import difflib,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
base=ROOT/'suprem4/src';candidate=HERE/'work/src'
paths=sorted({p.relative_to(base) for p in base.rglob('*') if p.suffix in ('.c','.h')}|{p.relative_to(candidate) for p in candidate.rglob('*') if p.suffix in ('.c','.h')})
parts=[];records=[]
sha=lambda b:hashlib.sha256(b).hexdigest()
for rel in paths:
 a=base/rel;b=candidate/rel
 av=a.read_bytes() if a.exists() else b'';bv=b.read_bytes() if b.exists() else b''
 if av==bv:continue
 parts.extend(difflib.unified_diff(av.decode().splitlines(True),bv.decode().splitlines(True),fromfile='a/src/'+str(rel) if a.exists() else '/dev/null',tofile='b/src/'+str(rel) if b.exists() else '/dev/null'))
 records.append({'path':str(rel),'main_sha256':sha(av) if a.exists() else None,'experimental_sha256':sha(bv) if b.exists() else None})
patch=''.join(parts).encode();(HERE/'experimental_solver.patch').write_bytes(patch)
r={'status':'experimental, not accepted','base':'exact main working-tree file hashes below','files':records,'patch_sha256':sha(patch),'note':'Archive only. Includes experimental trust-region, matrix-capture, normal-layer and import-quality code, plus differences in diagnostics. Do not blanket-apply this overlay as a qualified repair. Frozen build manifests identify tested variants. No general launcher promotion.'}
(HERE/'experimental_solver_provenance.json').write_text(json.dumps(r,indent=2)+'\n');print(len(records),'source differences archived')
