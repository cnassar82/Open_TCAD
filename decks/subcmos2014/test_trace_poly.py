import tempfile
from pathlib import Path
import unittest
from trace_poly import inventory


class InventoryTest(unittest.TestCase):
    def test_linear_field_integral_and_surface_profile(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'linear.str'
            path.write_text('''c 1 0 0
c 2 10 0
c 3 10 8
c 4 0 8
r 1 3
t 1 1 1 2 3
t 2 1 1 3 4
s 2 5 3
n 0 3 1e15 1e16
n 1 3 1e15 11e16
n 2 3 1e15 11e16
n 3 3 1e15 1e16
''')
            result = inventory(path)
            si = result['materials']['3']
            self.assertAlmostEqual(si['area_um2'], 80)
            self.assertAlmostEqual(si['chemical_atoms_per_cm']['phosphorus']/4.8e10, 1.)
            self.assertAlmostEqual(si['chemical_atoms_per_cm']['boron']/8e8, 1.)
            body = result['probes']['body']['profile']['0.05']
            self.assertAlmostEqual(body['phosphorus_cm3']/2.25e16, 1.)


if __name__ == '__main__':
    unittest.main()
