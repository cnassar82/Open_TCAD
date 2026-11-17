from pathlib import Path
import numpy as np,json,time
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import splu
p=Path(__file__).resolve().parent/'fe_trace/pmos/FE_matrix.bin';raw=p.read_bytes();n,off,ni=map(int,np.frombuffer(raw[:12],np.int32));ia=np.frombuffer(raw[12:12+4*ni],np.int32);a=np.frombuffer(raw[12+4*ni:12+4*ni+8*(ni+off)],np.float64);rhs=np.frombuffer(raw[12+4*ni+8*(ni+off):],np.float64)
cols=np.repeat(np.arange(n),np.diff(ia[:n+1]));rows=ia[n+1:];A=coo_matrix((np.concatenate((a[:n],a[n+1:ni],a[n+1+off:ni+off])),(np.concatenate((np.arange(n),rows,cols)),np.concatenate((np.arange(n),cols,rows)))),shape=(n,n)).tocsc()
scale=1/np.asarray(abs(A).max(axis=1).toarray()).ravel();B=diags(scale)@A;cs=1/np.asarray(abs(B).max(axis=0).toarray()).ravel();B=(B@diags(cs)).tocsc();start=time.monotonic();lu=splu(B);x=cs*lu.solve(scale*rhs);res=A@x-rhs;den=abs(A)@abs(x)+abs(rhs);backward=float(np.max(abs(res)/np.maximum(den,np.finfo(float).tiny)))
scaled_normwise=float(np.max(abs(B@(x/cs)-scale*rhs))/(np.max(np.asarray(abs(B).sum(axis=1)).ravel())*np.max(abs(x/cs))+np.max(abs(scale*rhs))))
record=dict(scaled_normwise_backward_error=scaled_normwise,n=n,nnz=A.nnz,all_matrix_finite=bool(np.isfinite(A.data).all()),diagonal_min=float(abs(a[:n]).min()),diagonal_max=float(abs(a[:n]).max()),pivoted_solution_finite=bool(np.isfinite(x).all()),componentwise_backward_error=backward,elapsed_seconds=time.monotonic()-start)
print(record,flush=True);(p.parent/'pivoted_oracle.json').write_text(json.dumps(record,indent=2)+'\n')
