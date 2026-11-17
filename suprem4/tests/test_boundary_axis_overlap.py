"""Independent boundary ownership cases from the RIT spacer failure."""
import unittest

from run_oxidation_perturbations import boundary_axis_overlap


class BoundaryAxisOverlapTest(unittest.TestCase):
    def test_shared_endpoint_is_valid(self):
        coordinates = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (2.0, 0.0)}
        boundary_axis_overlap(coordinates, [(1, 2), (2, 3)])

    def test_overlapping_exposed_edges_are_rejected(self):
        # A thin residual spacer oxide edge and old gate oxide edge both claim
        # the same exposed interval, with distinct point identities.
        coordinates = {1: (4.570329911075532, 0.07342216122196987),
                       2: (4.570341552607715, 0.07342216122196987),
                       3: (4.570333403535187, 0.07342216122196987)}
        with self.assertRaisesRegex(ValueError, "overlapping horizontal boundary"):
            boundary_axis_overlap(coordinates, [(1, 2), (3, 2)])

    def test_nearby_parallel_edges_are_valid(self):
        coordinates = {1: (0.0, 0.0), 2: (1.0, 0.0),
                       3: (0.0, 1.0e-9), 4: (1.0, 1.0e-9)}
        boundary_axis_overlap(coordinates, [(1, 2), (3, 4)])


if __name__ == "__main__":
    unittest.main()
