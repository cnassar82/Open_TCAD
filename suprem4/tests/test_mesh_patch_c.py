"""Exercise the production proposal kernel against independent geometry checks."""
import ctypes as C
import json
import math
import os
from functools import lru_cache
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest
from fractions import Fraction as F
from mesh_patch_oracle import propose_interior_star
from partition_audit import audit_partition

ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(os.environ.get('SUPREM_TEST_SOURCE',ROOT/'src')).resolve()
class Point(C.Structure):_fields_=[('xy',C.c_float*2)]
class Triangle(C.Structure):_fields_=[('vertex',C.c_int*3),('region',C.c_int)]
class Result(C.Structure):
    _fields_=[('count',C.c_int),('triangle',C.POINTER(Triangle)),
              ('oldq',C.c_double),('newq',C.c_double),('reason',C.c_char_p)]
def orient(p,t):
    a,b,c=[tuple(F(v) for v in p[i]) for i in t]
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])

class KernelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.work=Path(cls.tmp.name)
        so=cls.work/'kernel.so'
        misc=cls.work/'misc.o'
        subprocess.run(['gcc','-c','-fPIC','-O2','-std=gnu89','-ffunction-sections',
                        '-fvisibility=hidden','-I',str(SOURCE/'include'),
                        str(SOURCE/'geom/misc.c'),'-o',str(misc)],check=True)
        subprocess.run(['gcc','-shared','-fPIC','-O2','-std=gnu89','-ffunction-sections',
                        '-I',str(SOURCE/'include'),str(SOURCE/'dbase/mesh_patch.c'),
                        str(misc),'-Wl,--gc-sections','-lm','-o',str(so)],check=True)
        cls.lib=C.CDLL(str(so));cls.lib.mesh_patch_propose.argtypes=[C.c_int,C.c_int,C.POINTER(Point),C.c_int,C.POINTER(Triangle),C.POINTER(Result)]
        cls.lib.mesh_patch_release.argtypes=[C.POINTER(Result)]
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def run_kernel(self,p,t,regions=None,center=0):
        pts=(Point*len(p))(*(Point((C.c_float*2)(*xy)) for xy in p))
        rounded=[tuple(point.xy) for point in pts]
        tris=(Triangle*len(t))(*(Triangle((C.c_int*3)(*v),r) for v,r in zip(t,regions or [0]*len(t))))
        old=bytes(pts)+bytes(tris);result=Result()
        code=self.lib.mesh_patch_propose(center,len(p),pts,len(t),tris,C.byref(result))
        output=[tuple(result.triangle[i].vertex) for i in range(result.count)]
        reason=result.reason.decode();q=(result.oldq,result.newq)
        self.assertEqual(bytes(pts)+bytes(tris),old)
        self.lib.mesh_patch_release(C.byref(result))
        self.assertEqual(result.count,0)
        return code,output,reason,q,rounded
    def verify_output(self,p,old,new):
        self.assertEqual(sum(orient(p,t) for t in old),sum(orient(p,t) for t in new))
        from collections import Counter
        def boundary(ts):
            c=Counter((t[i],t[(i+1)%3]) for t in ts for i in range(3))
            return {e for e in c if c[e]-c[e[::-1]]==1}
        self.assertEqual(boundary(old),boundary(new))
        file=self.work/'output.str'
        file.write_text('\n'.join([f'c {i+1} {x!r} {y!r}' for i,(x,y) in enumerate(p)] +
            [f't {i+1} 1 '+' '.join(str(v+1) for v in t)+' 0 0 0' for i,t in enumerate(new)]))
        audit_partition(file)
    def test_convex_concave_and_collinear(self):
        for p in [[(0,0),(-1,-1),(1,-1),(1,1),(-1,1)],
                  [(.5,.5),(0,0),(3,0),(3,1),(1,1),(1,3),(0,3)],
                  [(0,0),(-2,-1),(2,-1),(2,1),(0,1),(-2,1)]]:
            t=[(0,i, i+1 if i<len(p)-1 else 1) for i in range(1,len(p))]
            code,new,reason,q,p=self.run_kernel(p,t)
            self.assertIn(code,(0,1),reason)
            if code==0:self.verify_output(p,t,new);self.assertGreaterEqual(q[1],q[0]*(1-1e-12))
    def test_generated_binary32_stars(self):
        rng=random.Random(15112026);accepted=0;refused=0
        for n in range(4,16):
            for rep in range(3):
                p=[(0,0)]+[(rng.uniform(.8,1.2)*math.cos(2*math.pi*i/n),
                             rng.uniform(.8,1.2)*math.sin(2*math.pi*i/n)) for i in range(n)]
                t=[(0,i+1,(i+1)%n+1) for i in range(n)];rng.shuffle(t)
                code,new,reason,q,p=self.run_kernel(p,t)
                self.assertIn(code,(0,1),reason)
                reference=propose_interior_star(0,dict(enumerate(p)),t)
                self.assertEqual(reference.status,'OK')
                if code==0:
                    accepted+=1;self.verify_output(p,t,new);self.assertGreaterEqual(q[1],q[0]*(1-1e-12))
                else:refused+=1;self.assertEqual(reason,'quality would deteriorate')
        self.assertGreater(accepted,0)
    def test_high_valence_quality_refusal(self):
        n=64
        p=[(0,0)]+[(math.cos(2*math.pi*i/n),math.sin(2*math.pi*i/n)) for i in range(n)]
        old=[(0,i+1,(i+1)%n+1) for i in range(n)]
        code,new,reason,q,_=self.run_kernel(p,old)
        self.assertEqual(code,1,reason)
        self.assertEqual(new,[])
        self.assertEqual(reason,'quality would deteriorate')
        self.assertTrue(math.isfinite(q[1]) and 0<q[1]<q[0])

    def test_captured_greedy_quality_deadlock(self):
        d=json.loads((ROOT/'tests/fixtures/rit_pmos_quality_deadlock/cavity.json').read_text())
        ids=[d['center']]+[int(k) for k in d['points_cm_hex'] if int(k)!=d['center']]
        p=[tuple(float.fromhex(x) for x in d['points_cm_hex'][str(k)]) for k in ids]
        old=[tuple(ids.index(v) for v in t) for t in d['old_triangles']]
        rng=random.Random(182026)
        for scale in (1.,1024.,1./1024):
            permutation=list(range(len(p)));rng.shuffle(permutation)
            points=[tuple(scale*x for x in p[v]) for v in permutation]
            triangles=[tuple(permutation.index(v) for v in t) for t in old];rng.shuffle(triangles)
            code,new,reason,q,rounded=self.run_kernel(points,triangles,center=permutation.index(0))
            self.assertEqual(code,0,reason)
            self.assertAlmostEqual(q[0],d['old_min_quality'],places=12)
            self.assertAlmostEqual(q[1],d['optimal_min_quality'],places=12)
            self.assertGreater(q[1],.8)
            self.verify_output(rounded,triangles,new)

    def test_optimal_score_against_exhaustive_ear_search(self):
        from mesh_patch_oracle import inside_or_on
        rng=random.Random(192026)
        for n in range(4,9):
            for rep in range(3):
                points=[(0,0)]+[(rng.uniform(.5,1.5)*math.cos(2*math.pi*i/n),
                                 rng.uniform(.5,1.5)*math.sin(2*math.pi*i/n)) for i in range(n)]
                old=[(0,i+1,(i+1)%n+1) for i in range(n)]
                code,new,reason,q,p=self.run_kernel(points,old)
                exact={i:tuple(F(x) for x in xy) for i,xy in enumerate(p)}
                def quality(t):
                    v=[exact[i] for i in t]
                    denominator=sum((a[0]-b[0])**2+(a[1]-b[1])**2 for a,b in zip(v,v[1:]+v[:1]))
                    return 2*math.sqrt(3)*float(orient(p,t)/denominator)
                @lru_cache(None)
                def exhaustive(poly):
                    if len(poly)==3:return quality(poly) if orient(p,poly)>0 else -1
                    best=-1
                    for i in range(len(poly)):
                        t=(poly[i-1],poly[i],poly[(i+1)%len(poly)])
                        if orient(p,t)<=0 or any(inside_or_on(exact,v,*t) for v in poly if v not in t):continue
                        rest=exhaustive(poly[:i]+poly[i+1:])
                        if rest>=0:best=max(best,min(quality(t),rest))
                    return best
                optimum=exhaustive(tuple(range(1,n+1)))
                self.assertIn(code,(0,1),reason)
                self.assertAlmostEqual(q[1],optimum,places=12)
                if code==0:self.verify_output(p,old,new)
                else:self.assertLess(optimum,q[0]*(1-1e-12))

    def test_captured_boundary_and_region_constraint(self):
        d=json.loads((ROOT/'tests/fixtures/rit_pmos_first_bad_removal/patch.json').read_text())
        ids=[d['center']['id']]+[int(k) for k in d['points'] if int(k)!=d['center']['id']]
        p=[tuple(float.fromhex(x) for x in d['points'][str(k)]) for k in ids]
        t=[tuple(ids.index(k) for k in tr['vertices']) for tr in d['old_triangles']]
        code,new,reason,_,_=self.run_kernel(p,t)
        self.assertEqual((code,new,reason),(1,[],'boundary constraint'))
        p=[(0,0),(-1,-1),(1,-1),(1,1),(-1,1)];t=[(0,1,2),(0,2,3),(0,3,4),(0,4,1)]
        self.assertEqual(self.run_kernel(p,t,[1,1,2,2])[2],'region constraint')
    def test_invalid_schema_and_nan(self):
        p=[(0,0),(-1,-1),(1,-1),(1,1),(-1,1)];t=[(0,1,2),(0,2,3),(0,3,4),(0,4,1)]
        for bad in [[(0,2,1)]+t[1:],t+[t[0]],[(0,1,99)]+t[1:]]:
            self.assertEqual(self.run_kernel(p,bad)[0],2)
        p[0]=(float('nan'),0)
        self.assertEqual(self.run_kernel(p,t)[0],2)
    def test_straight_boundary_subdivision(self):
        from collections import Counter
        p=[(-.99,0),(-1,0),(-1,-1),(1,-1),(1,0)]
        t=[(0,1,2),(0,2,3),(0,3,4)]
        code,new,reason,q,p=self.run_kernel(p,t)
        self.assertEqual(code,0,reason)
        self.assertEqual(len(new),len(t)-1)
        self.assertEqual(sum(orient(p,tr) for tr in t),sum(orient(p,tr) for tr in new))
        def boundary(ts):
            counts=Counter((tr[i],tr[(i+1)%3]) for tr in ts for i in range(3))
            return {e for e in counts if counts[e]-counts[e[::-1]]==1}
        expected=boundary(t)-{(4,0),(0,1)}|{(4,1)}
        self.assertEqual(boundary(new),expected)
        self.assertGreater(q[1],q[0])
        for tr in new:self.assertGreater(orient(p,tr),0)
        p[0]=(-.99,.01)
        self.assertEqual(self.run_kernel(p,t)[0],1)
        p=[(.01,0),(0,0),(0,-1),(2,0)]
        t=[(0,1,2),(0,2,3)]
        code,new,reason,_,p=self.run_kernel(p,t)
        self.assertEqual(code,0,reason)
        self.assertEqual(len(new),1)
        self.assertEqual(sum(orient(p,tr) for tr in t),orient(p,new[0]))

if __name__=='__main__':unittest.main()
