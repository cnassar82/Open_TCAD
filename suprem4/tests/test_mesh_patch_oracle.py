import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from mesh_patch_oracle import propose_interior_star, orient


class PatchOracleTest(unittest.TestCase):
    def test_convex_square(self):
        p={0:(0,0),1:(-1,-1),2:(1,-1),3:(1,1),4:(-1,1)}
        old=[(0,1,2),(0,2,3),(0,3,4),(0,4,1)]
        r=propose_interior_star(0,p,old)
        self.assertEqual(r.status,'OK')
        self.assertEqual(len(r.triangles),2)
        self.assertTrue(all(orient(p,*t)>0 for t in r.triangles))
        self.assertNotIn(0,{v for t in r.triangles for v in t})
        self.assertEqual(propose_interior_star(0,p,old[2:]+old[:2]).status,'OK')

    def test_concave_star(self):
        p={0:(Fraction(1,2),Fraction(1,2)),1:(0,0),2:(3,0),3:(3,1),
           4:(1,1),5:(1,3),6:(0,3)}
        old=[(0,a,b) for a,b in zip((1,2,3,4,5,6),(2,3,4,5,6,1))]
        r=propose_interior_star(0,p,old)
        self.assertEqual(r.status,'OK',r.reason)
        self.assertEqual(len(r.triangles),4)

    def test_high_valence_and_near_collinearity(self):
        p={0:(0,0)}
        for i in range(64):
            angle=2*math.pi*i/64
            p[i+1]=(math.cos(angle),math.sin(angle))
        old=[(0,i+1,(i+1)%64+1) for i in range(64)]
        r=propose_interior_star(0,p,old)
        self.assertEqual(r.status,'OK',r.reason)
        self.assertEqual(len(r.triangles),62)
        q={0:(0,0),1:(-2,-1),2:(2,-1),3:(2,1),4:(1,1+2**-20),5:(-2,1)}
        star=[(0,i,j) for i,j in ((1,2),(2,3),(3,4),(4,5),(5,1))]
        r=propose_interior_star(0,q,star)
        self.assertEqual(r.status,'OK',r.reason)

    def test_first_rit_failure_is_open_and_unsupported(self):
        f=Path(__file__).resolve().parent/'fixtures/rit_pmos_first_bad_removal/patch.json'
        data=json.loads(f.read_text())
        p={int(k):tuple(float.fromhex(x) for x in xy) for k,xy in data['points'].items()}
        r=propose_interior_star(data['center']['id'],p,
                                [t['vertices'] for t in data['old_triangles']])
        self.assertEqual((r.status,r.reason),('UNSUPPORTED','open boundary fan'))

    def test_rejects_nonmanifold_and_nonpositive(self):
        p={0:(0,0),1:(-1,-1),2:(1,-1),3:(1,1),4:(-1,1)}
        old=[(0,1,2),(0,2,3),(0,3,4),(0,4,1)]
        self.assertEqual(propose_interior_star(0,p,old+[old[0]]).status,'INVALID')
        self.assertEqual(propose_interior_star(0,p,[(0,2,1)]+old[1:]).status,'INVALID')

if __name__=='__main__':unittest.main()
