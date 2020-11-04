import unittest
from math import inf as infinity

from source.constraints import Constraints


class TestConstraints(unittest.TestCase):
    def setUp(self):
        pass

    def test_from_string(self):
        test_cases = [
            (True, True, 50.7, -12.34, 10, 2, 50, 20),
            (True, True, -16.7, -40.001, 1, 1, 2, 2),
        ]
        for (optimize, sort, max_energy, min_energy, max_polymers, min_polymers, max_merges, min_merges) in test_cases:
            test_constraints = Constraints.from_string(
                f"""
                {'NO ' if not optimize else ''}OPTIMIZE
                {'NO ' if not sort else ''}SORT
                MAX ENERGY {max_energy}
                MIN ENERGY {min_energy}
                MAX POLYMERS {max_polymers}
                MIN POLYMERS {min_polymers}
                MAX MERGES {max_merges}
                MIN MERGES {min_merges}
                """
            )
            with self.subTest(optimize=optimize):
                self.assertEqual(optimize, test_constraints.optimize())
            with self.subTest(sort=sort):
                self.assertEqual(sort, test_constraints.sort())
            with self.subTest(max_energy=max_energy, min_energy=min_energy):
                self.assertEqual(max_energy, test_constraints.max_energy())
                self.assertEqual(min_energy, test_constraints.min_energy())
            with self.subTest(max_energy=max_polymers, min_energy=min_polymers):
                self.assertEqual(max_polymers, test_constraints.max_polymers())
                self.assertEqual(min_polymers, test_constraints.min_polymers())
            with self.subTest(max_energy=max_merges, min_energy=min_merges):
                self.assertEqual(max_merges, test_constraints.max_merges())
                self.assertEqual(min_merges, test_constraints.min_merges())

    def test_add_new_constraint_from_string(self):
        constraints = Constraints.from_string("OPTIMIZE")
        self.assertTrue(constraints.optimize())
        constraints.add_new_constraint_from_string("NO OPTIMIZE")
        self.assertFalse(constraints.optimize())

    def test_with_fixed_polymers(self):
        test_cases = [5, 7]
        old_constraints = Constraints()
        for fixed_number_of_polymers in test_cases:
            new_constraints = old_constraints.with_fixed_polymers(fixed_number_of_polymers)
            self.assertEqual(fixed_number_of_polymers, new_constraints.max_polymers())
            self.assertEqual(fixed_number_of_polymers, new_constraints.min_polymers())
            self.assertNotEqual(fixed_number_of_polymers, old_constraints.max_polymers())
            self.assertNotEqual(fixed_number_of_polymers, old_constraints.min_polymers())
            old_constraints = new_constraints

    def test_with_fixed_merges(self):
        test_cases = [4, 6]
        old_constraints = Constraints()
        for fixed_number_of_merges in test_cases:
            new_constraints = old_constraints.with_fixed_merges(fixed_number_of_merges)
            self.assertEqual(fixed_number_of_merges, new_constraints.max_merges())
            self.assertEqual(fixed_number_of_merges, new_constraints.min_merges())
            self.assertNotEqual(fixed_number_of_merges, old_constraints.max_merges())
            self.assertNotEqual(fixed_number_of_merges, old_constraints.min_merges())
            old_constraints = new_constraints

    def test_with_fixed_energy(self):
        test_cases = [5.34, -7.12]
        old_constraints = Constraints()
        for fixed_amount_of_energy in test_cases:
            new_constraints = old_constraints.with_fixed_energy(fixed_amount_of_energy)
            self.assertEqual(fixed_amount_of_energy, new_constraints.max_energy())
            self.assertEqual(fixed_amount_of_energy, new_constraints.min_energy())
            self.assertNotEqual(fixed_amount_of_energy, old_constraints.max_energy())
            self.assertNotEqual(fixed_amount_of_energy, old_constraints.min_energy())
            old_constraints = new_constraints

    def test_with_unset_optimization_flag(self):
        old_constraints = Constraints.from_string("OPTIMIZE")
        self.assertTrue(old_constraints.optimize())
        new_constraints = old_constraints.with_unset_optimization_flag()
        self.assertFalse(new_constraints.optimize())
        self.assertTrue(old_constraints.optimize())
