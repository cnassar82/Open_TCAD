"""Standalone exact 2D cavity proposal oracle; no simulator mutation.

The supported domain is one interior point in one region, with a simple
triangular star. This is a reference for a future production proposal kernel,
not a simulator replacement or a field-transfer implementation.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Proposal:
    status: str
    reason: str
    triangles: tuple = ()
    boundary: tuple = ()


def orient(points, a, b, c):
    (ax, ay), (bx, by), (cx, cy) = (points[p] for p in (a,b,c))
    return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax)


def inside_or_on(points, p, a, b, c):
    return all(orient(points,*edge,p) >= 0 for edge in ((a,b),(b,c),(c,a)))


def _cross(points,a,b,c,d):
    return (orient(points,a,b,c)*orient(points,a,b,d)<0 and
            orient(points,c,d,a)*orient(points,c,d,b)<0)


def _interiors_overlap(points, a, b):
    return not any(all(orient(points,u,v,w)<=0 for w in b)
                   for u,v in zip(a,a[1:]+a[:1])) and not any(
                   all(orient(points,u,v,w)<=0 for w in a)
                   for u,v in zip(b,b[1:]+b[:1]))


def _on_segment(points,p,a,b):
    if orient(points,a,b,p)!=0:return False
    x,y=points[p];ax,ay=points[a];bx,by=points[b]
    return min(ax,bx)<=x<=max(ax,bx) and min(ay,by)<=y<=max(ay,by)


def _bad_boundary_pair(points,a,b,c,d):
    if _cross(points,a,b,c,d):return True
    if all(orient(points,a,b,p)==0 for p in (c,d)):
        axis=0 if points[a][0]!=points[b][0] else 1
        lo=max(min(points[a][axis],points[b][axis]),min(points[c][axis],points[d][axis]))
        hi=min(max(points[a][axis],points[b][axis]),max(points[c][axis],points[d][axis]))
        if hi>lo:return True
    return any(_on_segment(points,p,u,v) and p not in (u,v)
               for p,u,v in ((a,c,d),(b,c,d),(c,a,b),(d,a,b)))


def propose_interior_star(center, points, old_triangles):
    """Return a complete read-only triangulation or an explicit refusal."""
    points = {p:tuple(Fraction(v) for v in xy) for p,xy in points.items()}
    if not old_triangles or center not in points:
        return Proposal('INVALID','missing center or triangles')
    links = {}
    incoming = Counter()
    old_area2 = Fraction(0)
    for tri in old_triangles:
        if len(tri)!=3 or len(set(tri))!=3 or center not in tri or any(p not in points for p in tri):
            return Proposal('INVALID','triangle identity or center')
        area = orient(points,*tri)
        if area<=0:return Proposal('INVALID','nonpositive old triangle')
        old_area2 += area
        i=tri.index(center);a,b=tri[(i+1)%3],tri[(i+2)%3]
        if a in links:return Proposal('INVALID','nonmanifold outgoing ray')
        links[a]=b;incoming[b]+=1
        if incoming[b]>1:return Proposal('INVALID','nonmanifold incoming ray')
    for i,a in enumerate(old_triangles):
        for b in old_triangles[i+1:]:
            if _interiors_overlap(points,a,b):
                return Proposal('INVALID','overlapping old triangles')
    if set(links)!=set(incoming):
        return Proposal('UNSUPPORTED','open boundary fan')
    first=min(links);order=[first];seen={first}
    while True:
        nxt=links[order[-1]]
        if nxt==first:break
        if nxt in seen or len(order)>=len(links):
            return Proposal('INVALID','disconnected or repeated fan')
        order.append(nxt);seen.add(nxt)
    if len(order)!=len(links):return Proposal('INVALID','disconnected fan')
    if len(order)<3:return Proposal('INVALID','fewer than three boundary points')
    for i in range(len(order)):
        a,b=order[i],order[(i+1)%len(order)]
        for j in range(i+2,len(order)):
            if i==0 and j==len(order)-1:continue
            c,d=order[j],order[(j+1)%len(order)]
            if _bad_boundary_pair(points,a,b,c,d):
                return Proposal('INVALID','crossing or contacting old boundary')
    polygon_area2=sum((points[order[i]][0]*points[order[(i+1)%len(order)]][1]-
                       points[order[i]][1]*points[order[(i+1)%len(order)]][0])
                      for i in range(len(order)))
    if polygon_area2!=old_area2:return Proposal('INVALID','old star does not match boundary domain')
    work=list(order);out=[]
    while len(work)>3:
        cut=False
        for i in range(len(work)):
            a,b,c=work[i-1],work[i],work[(i+1)%len(work)]
            if orient(points,a,b,c)<=0:continue
            if any(inside_or_on(points,p,a,b,c) for p in work if p not in (a,b,c)):
                continue
            out.append((a,b,c));work.pop(i);cut=True;break
        if not cut:return Proposal('UNSUPPORTED','no admissible ear')
    if orient(points,*work)<=0:return Proposal('INVALID','nonpositive final ear')
    out.append(tuple(work))
    if sum(orient(points,*tri) for tri in out)!=old_area2:
        return Proposal('INVALID','candidate area differs')
    for i,a in enumerate(out):
        for b in out[i+1:]:
            if _interiors_overlap(points,a,b):
                return Proposal('INVALID','overlapping candidate triangles')
    edges=Counter(tuple(sorted((tri[i],tri[(i+1)%3]))) for tri in out for i in range(3))
    boundary={tuple(sorted((order[i],order[(i+1)%len(order)]))) for i in range(len(order))}
    if {e for e,c in edges.items() if c==1}!=boundary or any(c>2 for c in edges.values()):
        return Proposal('INVALID','candidate boundary differs')
    return Proposal('OK','exact interior cavity',tuple(out),tuple(order))
