import tempfile
from pathlib import Path
import unittest

from partition_audit import audit_partition


class PartitionAuditTest(unittest.TestCase):
    def check_mesh(self, triangles):
        coords = list(dict.fromkeys(p for triangle in triangles for p in triangle))
        vertices = [[coords.index(p) + 1 for p in tri] for tri in triangles]
        owners = {}
        for i, tri in enumerate(vertices, 1):
            for j in range(3):
                edge = tuple(sorted((tri[(j + 1) % 3], tri[(j + 2) % 3])))
                owners.setdefault(edge, []).append(i)
        records = [f"c {i} {p[0]:.17g} {p[1]:.17g} 0" for i, p in enumerate(coords, 1)]
        for i, tri in enumerate(vertices, 1):
            neighbors = []
            for j in range(3):
                edge = tuple(sorted((tri[(j + 1) % 3], tri[(j + 2) % 3])))
                neighbors.append(next((t for t in owners[edge] if t != i), -1022))
            records.append("t " + " ".join(map(str, [i, 1, *tri, *neighbors])))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "mesh.str"
            path.write_text("\n".join(records) + "\n")
            return audit_partition(path)

    def test_conforming_diagonal_interface(self):
        result = self.check_mesh([[(0, 0), (1, 0), (0, 1)],
                                  [(1, 0), (1, 1), (0, 1)]])
        self.assertEqual(result["triangles"], 2)

    def test_contained_triangle(self):
        with self.assertRaisesRegex(ValueError, "overlapping triangle interiors"):
            self.check_mesh([[(0, 0), (4, 0), (0, 4)],
                             [(1, 1), (2, 1), (1, 2)]])

    def test_crossed_triangles(self):
        with self.assertRaisesRegex(ValueError, "overlapping triangle interiors"):
            self.check_mesh([[(0, 0), (2, 0), (1, 2)],
                             [(0, 1), (1, -1), (2, 1)]])

    def test_unnoded_collinear_interface(self):
        with self.assertRaisesRegex(ValueError, "boundary"):
            self.check_mesh([[(0, 0), (2, 0), (0, 2)],
                             [(1, 1), (2, 2), (0, 2)]])

    def test_positive_gap_remains_valid(self):
        self.check_mesh([[(0, 0), (1, 0), (0, 1)],
                         [(1 + 1e-12, 0), (2, 0), (2, 1)]])


if __name__ == "__main__":
    unittest.main()
