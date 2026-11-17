"""Extract the unique exposed boundary chain without coordinate snapping."""
import gzip
from pathlib import Path
import sys
p=Path(sys.argv[1]);raw=gzip.decompress(p.read_bytes()).decode() if p.suffix=='.gz' else p.read_text()
coordinates={};adj={};edges=set()
for line in raw.splitlines():
 v=line.split()
 if not v:continue
 if v[0]=='c':coordinates[int(v[1])]=tuple(map(float,v[2:4]))
 if v[0]=='t':
  vs=list(map(int,v[3:6]))
  for i,n in enumerate(map(int,v[6:9])):
   if n!=-1022:continue
   a,b=vs[(i+1)%3],vs[(i+2)%3];edge=tuple(sorted((a,b)))
   if edge in edges:raise ValueError('duplicate exposed edge')
   edges.add(edge);adj.setdefault(a,[]).append(b);adj.setdefault(b,[]).append(a)
ends=[a for a,b in adj.items() if len(b)==1]
if len(ends)!=2 or any(len(b)>2 for b in adj.values()):raise ValueError('exposed boundary is not a single open chain')
start=min(ends,key=lambda a:coordinates[a][0]);previous=None;current=start;visited=set()
while current is not None:
 if current in visited:raise ValueError('loop in exposed chain')
 visited.add(current);print(*coordinates[current])
 following=[n for n in adj[current] if n!=previous]
 previous,current=current,following[0] if following else None
if len(visited)!=len(adj):raise ValueError('disconnected exposed boundary')
