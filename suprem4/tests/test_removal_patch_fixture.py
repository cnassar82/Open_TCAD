"""Exact-geometry oracle for the first captured invalid PMOS removal.

Normal invocation verifies the extracted patch. --legacy-control intentionally
exits nonzero if the legacy angle-order proposal crosses its own boundary.
"""
import json
import sys
import unittest
from collections import Counter
from fractions import Fraction
from pathlib import Path

FIXTURE = Path(__file__).resolve().parent / 'fixtures/rit_pmos_first_bad_removal/patch.json'
DATA = json.loads(FIXTURE.read_text())
POINTS = {int(k): tuple(Fraction.from_float(float.fromhex(x)) for x in v)
          for k, v in DATA['points'].items()}
TRIANGLES = [t['vertices'] for t in DATA['old_triangles']]


def orient(a, b, c):
    a, b, c = (POINTS[p] for p in (a, b, c))
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def proper_cross(a, b, c, d):
    return (orient(a,b,c)*orient(a,b,d) < 0 and
            orient(c,d,a)*orient(c,d,b) < 0)


def boundary_edges():
    edges = Counter(tuple(sorted((t[i],t[(i+1)%3]))) for t in TRIANGLES
                    for i in range(3))
    return [edge for edge, count in edges.items() if count == 1]


def proposed_crossing(order):
    edges = [(order[i],order[(i+1)%len(order)]) for i in range(len(order))]
    return [(a,b,c,d) for i,(a,b) in enumerate(edges)
            for j,(c,d) in enumerate(edges) if j > i+1 and
            not (i == 0 and j == len(order)-1) and proper_cross(a,b,c,d)]


class FirstBadRemoval(unittest.TestCase):
    def test_original_star_is_valid_and_open(self):
        center = DATA['center']['id']
        self.assertEqual(len(TRIANGLES), 3)
        self.assertTrue(all(orient(*t)>0 for t in TRIANGLES))
        boundary = boundary_edges()
        self.assertEqual(len(boundary),5)
        rays = [e for e in boundary if center in e]
        self.assertEqual(len(rays),2)
        ends = [next(p for p in e if p != center) for e in rays]
        self.assertEqual(set(ends), {3793,3797})
        self.assertNotEqual(orient(ends[0],center,ends[1]),0)

    def test_complete_point_ownership_is_part_of_removal(self):
        owner = DATA['entry_ownership']
        self.assertTrue(owner['surface'])
        self.assertEqual([n['material'] for n in owner['nodes']], [0,1])
        self.assertEqual(DATA['rem_pt_order'], [5101,5102])
        self.assertEqual(DATA['call']['edge_id'], 12550)
        self.assertEqual(DATA['call']['removed_point'], DATA['center']['id'])

    def test_legacy_proposal_crosses_original_outer_edge(self):
        crossings = proposed_crossing(DATA['legacy_angle_order'])
        self.assertEqual(len(crossings),1)
        edges = {frozenset(crossings[0][:2]),frozenset(crossings[0][2:])}
        self.assertEqual(edges,{frozenset((3567,3566)),frozenset((3793,3797))})
        self.assertIn((3567,3566), [tuple(e) for e in boundary_edges()] +
                      [tuple(reversed(e)) for e in boundary_edges()])
        # Closing a noncollinear boundary fan changes the exact material domain.
        self.assertNotEqual(orient(3793,3796,3797),0)


if __name__ == '__main__':
    if '--legacy-control' in sys.argv:
        crosses = proposed_crossing(DATA['legacy_angle_order'])
        print(f'legacy proposal crossings: {crosses}')
        sys.exit(bool(crosses))
    unittest.main()
