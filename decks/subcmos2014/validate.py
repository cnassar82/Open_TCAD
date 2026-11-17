"""Independent dopant polarity and geometry probes for the generated process cuts."""
import math
def read(path):
    coords={}; regions={}; triangles=[]; nodes={}; species=[]
    for line in path.read_text().splitlines():
        p=line.split()
        if not p:continue
        if p[0]=='c':coords[int(p[1])]=tuple(map(float,p[2:4]))
        elif p[0]=='r':regions[int(p[1])]=int(p[2])
        elif p[0]=='t':triangles.append((regions[int(p[2])],tuple(map(int,p[3:6]))))
        elif p[0]=='s':species=list(map(int,p[2:]))
        elif p[0]=='n':nodes[int(p[1])+1,int(p[2])]=dict(zip(species,map(float,p[3:])))
    return coords,triangles,nodes


def intervals(mesh,x,material):
    coords,triangles,_=mesh; spans=[]
    for mat,vs in triangles:
        if mat!=material:continue
        pts=[coords[v] for v in vs];ys=[]
        for a,b in zip(pts,pts[1:]+pts[:1]):
            if a[0]==b[0]:
                if abs(a[0]-x)<1e-10:ys.extend((a[1],b[1]))
            elif min(a[0],b[0])<=x<=max(a[0],b[0]):
                ys.append(a[1]+(b[1]-a[1])*(x-a[0])/(b[0]-a[0]))
        if ys:spans.append((min(ys),max(ys)))
    merged=[]
    for a,b in sorted(spans):
        if merged and a<=merged[-1][1]+1e-7:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
        else:merged.append((a,b))
    return merged


def concentration(mesh,x,y,material=3):
    coords,triangles,nodes=mesh
    for mat,vs in triangles:
        if mat!=material:continue
        a,b,c=[coords[v] for v in vs]
        det=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        u=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/det
        v=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/det
        w=1-u-v
        if min(u,v,w)>=-1e-8:
            values=[nodes[vertex,mat] for vertex in vs]
            def interp(s):return sum(t*d.get(s,0) for t,d in zip((u,v,w),values))
            return {'boron_cm3':interp(5),'phosphorus_cm3':interp(3),
                    'net_donor_cm3':interp(3)-interp(5)}
    raise ValueError(f'probe outside material {material}: ({x},{y})')


def gate_dielectric(mesh):
    """Require oxide material, rather than just a positive poly/silicon gap."""
    coords,tris,_=mesh
    poly=[coords[v] for mat,vs in tris if mat==4 for v in vs]
    if not poly:raise ValueError('missing poly gate')
    left,right=min(p[0] for p in poly),max(p[0] for p in poly)
    result=[]
    for fraction in (.1,.25,.5,.75,.9):
        x=left+(right-left)*fraction
        ps=intervals(mesh,x,4);ss=intervals(mesh,x,3)
        if not ps or not ss:raise ValueError(f'missing gate stack at x={x} um')
        bottom,top=ps[0][1],ss[0][0]
        tolerance=8*math.ulp(max(abs(bottom),abs(top),1e-12))
        if top<=bottom:raise ValueError(f'no positive gate dielectric at x={x} um')
        if not any(a<=bottom+tolerance and b>=top-tolerance for a,b in intervals(mesh,x,1)):
            raise ValueError(f'missing oxide material in gate gap at x={x} um: poly bottom={bottom}, silicon top={top}')
        result.append({'x_um':x,'oxide_thickness_um':top-bottom})
    # Stored triangles are positively oriented. A boundary edge directed
    # toward decreasing x has a downward-facing outward normal. Require an
    # oxide neighbor along every such gate edge, including between probes.
    # Use exact represented coordinates, allowing distinct material nodes at
    # the same interface point without accepting a finite gap.
    edges={}
    for mat,vs in tris:
        for a,b in zip(vs,vs[1:]+vs[:1]):
            key=tuple(sorted((coords[a],coords[b])))
            edges.setdefault(key,[]).append(mat)
    underside_edges=0
    for mat,vs in tris:
        if mat!=4:continue
        for a,b in zip(vs,vs[1:]+vs[:1]):
            owners=edges[tuple(sorted((coords[a],coords[b])))]
            if coords[b][0]<coords[a][0] and owners.count(4)==1:
                underside_edges+=1
                if sorted(owners)!=[1,4]:
                    raise ValueError(f'missing oxide material along gate underside edge {coords[a]} to {coords[b]}: owners={owners}')
    if not underside_edges:raise ValueError('no downward-facing gate interface edges')
    return result


def validate(path,device,stage):
    mesh=read(path);result={}
    def surface(x):return intervals(mesh,x,3)[0][0]
    if stage=='wells':
        probe=concentration(mesh,6,surface(6)+.5)
        expected=1 if device=='pmos' else -1
        if probe['net_donor_cm3']*expected<=0:raise ValueError(f'wrong well polarity: {probe}')
        result['well_0p5um_below_surface']=probe
        result['deep_substrate']=concentration(mesh,6,7)
    if stage in ('gate','annealed','contacts','metal1'):
        coords,tris,_=mesh
        poly=[coords[v] for mat,vs in tris if mat==4 for v in vs]
        if not poly:raise ValueError('missing poly gate')
        width=max(p[0] for p in poly)-min(p[0] for p in poly)
        if abs(width-1)>1e-4:raise ValueError(f'gate width {width} um, expected 1 um')
        result['gate_length_um']=width
        spans=intervals(mesh,6,4)
        probe=concentration(mesh,6,sum(spans[0])/2,4)
        if probe['net_donor_cm3']<=0:raise ValueError('poly gate is not n-type')
        result['poly_midpoint']=probe
        # y increases into the substrate: oxide separates poly bottom from silicon.
        result['gate_dielectric_um']=surface(6)-spans[0][1]
        if result['gate_dielectric_um']<=0:raise ValueError('missing gate dielectric')
        result['gate_oxide_material_probes']=gate_dielectric(mesh)
    if stage in ('source_drain','annealed','metal1'):
        expected=1 if device=='nmos' else -1
        for name,x,sign in [('source',4.25,expected),('drain',7.75,expected),('body',1.25,-expected)]:
            probe=concentration(mesh,x,surface(x)+.05)
            if sign*probe['net_donor_cm3']<=0:raise ValueError(f'wrong {name} polarity: {probe}')
            result[name+'_50nm_below_surface']=probe
    if stage=='metal1':
        # Count connected metal islands by shared vertices.
        _,tris,_=mesh;parent={}
        def root(v):
            while parent[v]!=v:parent[v]=parent[parent[v]];v=parent[v]
            return v
        for mat,vs in tris:
            if mat!=6:continue
            for v in vs:parent.setdefault(v,v)
            for v in vs[1:]:parent[root(v)]=root(vs[0])
        count=len({root(v) for v in parent})
        if count!=4:raise ValueError(f'{count} metal islands, expected body/source/gate/drain')
        result['metal_islands']=count
        for name,x,mat in [('body',1.25,3),('source',4.25,3),('gate',6,4),('drain',7.75,3)]:
            metal=intervals(mesh,x,6);target=intervals(mesh,x,mat)
            if not metal or not target:raise ValueError(f'missing {name} contact')
            gap=target[0][0]-metal[-1][1]
            if abs(gap)>1e-5:raise ValueError(f'{name} contact gap {gap} um')
            result[name+'_contact_gap_um']=gap
    return result
