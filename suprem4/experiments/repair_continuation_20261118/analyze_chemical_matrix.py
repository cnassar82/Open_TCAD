#!/usr/bin/env python3
"""Compare the untouched coupled chemical Jacobian with an independent LU."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);args=p.parse_args()
with (args.folder/'chemical_matrix.bin').open('rb') as f:
    magic,nv,ns=np.fromfile(f,np.int32,3);assert magic==0x43484d31
    species=np.fromfile(f,np.int32,ns)
    rhs=np.fromfile(f,np.float64,nv*ns)
    current=np.fromfile(f,np.float64,nv*ns)
    rows=[];cols=[];values=[]
    for bi in range(ns):
        for bj in range(ns):
            kind,length,offset,index_length=np.fromfile(f,np.int32,4)
            index=np.fromfile(f,np.int32,index_length)
            data=np.fromfile(f,np.float64,length)
            if not kind:continue
            ids=np.arange(nv)
            rows.append(bi*nv+ids);cols.append(bj*nv+ids);values.append(data[:nv])
            if kind!=1:
                ids=np.repeat(np.arange(nv),np.diff(index[:nv+1]))
                neighbors=index[index[0]:index[nv]];entries=np.arange(index[0],index[nv])
                rows.extend((bi*nv+neighbors,bi*nv+ids))
                cols.extend((bj*nv+ids,bj*nv+neighbors))
                values.extend((data[entries],data[entries+offset]))
    assert not f.read(1),'unexpected trailing matrix data'
a=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(nv*ns,nv*ns)).tocsc()
a.sum_duplicates();a.eliminate_zeros()
assert np.isfinite(a.data).all() and np.isfinite(rhs).all() and np.isfinite(current).all()
characteristic=np.maximum(abs(current),1.)
for bi,imp in enumerate(species):
    if imp==14:characteristic[bi*nv:(bi+1)*nv]=1.  # electrostatic potential (Psi)
ac=a@diags(characteristic)
maximum=abs(ac).max(axis=1).toarray().ravel();assert (maximum>0).all(),'zero Jacobian row'
row=1/maximum;scaled=(diags(row)@ac).tocsc();b=rhs*row
print(f'Captured {nv} nodes, species {species.tolist()}, {a.shape[0]} equations, {a.nnz} nonzeros',flush=True)
factor=splu(scaled)
solution=factor.solve(b)
# Residual correction uses the same captured matrix, without a diagonal floor.
for _ in range(3):
    residual=b-scaled@solution
    solution+=factor.solve(residual)
delta=solution*characteristic
legacy=np.fromfile(args.folder/'chemical_legacy_update.bin',np.float64)
assert legacy.shape==current.shape
def errors(y):
    residual=scaled@y-b
    denominator=abs(scaled)@abs(y)+abs(b)
    component=np.divide(abs(residual),denominator,out=np.zeros_like(residual),where=denominator!=0)
    normwise=np.max(abs(residual))/(np.max(abs(scaled).sum(axis=1))*np.max(abs(y))+np.max(abs(b)))
    return {'normwise_backward_error':float(normwise),'componentwise_backward_error':float(max(component)),
            'scaled_relative_residual':float(np.linalg.norm(residual)/np.linalg.norm(b))}
report={'nodes':int(nv),'species':species.tolist(),'equations':a.shape[0],'nonzeros':int(a.nnz),
        'legacy':errors(legacy/characteristic),'independent_pivoted':errors(solution),'updates':[]}
for bi,imp in enumerate(species):
    sl=slice(bi*nv,(bi+1)*nv);c=current[sl];old=legacy[sl];new=delta[sl]
    worst=int(np.argmin(c+old))
    record={'species':int(imp),'legacy_negative_count':int(np.count_nonzero(c+old<0)),
            'pivoted_negative_count':int(np.count_nonzero(c+new<0)),
            'worst_legacy_node':worst,'current':float(c[worst]),
            'legacy_update':float(old[worst]),'pivoted_update_same_node':float(new[worst])}
    if imp==3:
        record['captured_tail_node_31586']={'current':float(c[31586]),'legacy_update':float(old[31586]),'pivoted_update':float(new[31586])}
    report['updates'].append(record)
(args.folder/'linear_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
