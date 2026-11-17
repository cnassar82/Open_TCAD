#!/usr/bin/env python3
"""Diagnostic only: quotient adjacent same-material disk regions.
No point, field, triangle, neighbor, or material domain changes. Tests whether
artificial region constraints are responsible for immovable surface slivers.
"""
import argparse,hashlib,json
from collections import defaultdict
from pathlib import Path

def coalesce(source,target):
    lines=source.read_text().splitlines();materials={};tris={};edgeowners=defaultdict(list)
    for line in lines:
        s=line.split()
        if s and s[0]=='r':materials[int(s[1])]=int(s[2])
        elif s and s[0]=='t':
            t,r=int(s[1]),int(s[2]);v=tuple(map(int,s[3:6]));tris[t]=(r,v)
            for j in range(3):edgeowners[tuple(sorted((v[j],v[(j+1)%3])))].append((t,v[j],v[(j+1)%3]))
    assert all(len(owners)<=2 for owners in edgeowners.values())
    parent={r:r for r in materials};members={r:{r} for r in materials}
    def group(r):
        while parent[r]!=r:r=parent[r]
        return r
    def disk(regions):
        directed=[]
        for owners in edgeowners.values():
            inside=[x for x in owners if tris[x[0]][0] in regions]
            if len(inside)==1:directed.append(inside[0][1:])
        nxt={};incoming=set()
        for a,b in directed:
            if a in nxt or b in incoming:return False
            nxt[a]=b;incoming.add(b)
        if not nxt or set(nxt)!=incoming:return False
        start=next(iter(nxt));v=start;seen=set()
        while v not in seen:seen.add(v);v=nxt[v]
        return v==start and len(seen)==len(nxt)
    candidates=[]
    for e,owners in edgeowners.items():
        if len(owners)!=2:continue
        a,b=(tris[o[0]][0] for o in owners)
        if a!=b and materials[a]==materials[b]:candidates.append(tuple(sorted((a,b))))
    merges=[];blocked=[]
    for a,b in sorted(set(candidates)):
        a,b=group(a),group(b)
        if a==b:continue
        if disk(members[a]|members[b]):
            merges.append([sorted(members[a]),sorted(members[b])]);parent[b]=a;members[a]|=members.pop(b)
        else:blocked.append([a,b])
    ids={r:i+1 for i,r in enumerate(sorted(members))};mapping={r:ids[group(r)] for r in materials}
    output=[]
    for line in lines:
        s=line.split()
        if s and s[0]=='r':
            r=int(s[1])
            if r!=group(r):continue
            s[1]=str(mapping[r]);line=' '.join(s)
        elif s and s[0]=='t':s[2]=str(mapping[int(s[2])]);line=' '.join(s)
        output.append(line)
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text('\n'.join(output)+'\n')
    # Comparison ignores only the region labels, retaining material identity.
    def semantic(text):
        mats={int(s[1]):int(s[2]) for l in text.splitlines() if (s:=l.split()) and s[0]=='r'}
        out=[]
        for l in text.splitlines():
            s=l.split()
            if s and s[0]=='r':continue
            if s and s[0]=='t':s[2]=str(mats[int(s[2])])
            out.append(s)
        return out
    assert semantic(source.read_text())==semantic(target.read_text())
    report={'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'region_mapping':mapping,'merges':merges,'rejected_non_disk_unions':blocked,
            'regions_before':len(materials),'regions_after':len(members),
            'material_geometry_and_fields_identical':True,'diagnostic_only':True}
    target.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('target',type=Path)
    a=p.parse_args();r=coalesce(a.source,a.target);print(json.dumps(r,indent=2))
