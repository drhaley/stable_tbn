import unittest
from math import inf as infinity
from source.positive_multiset import PositiveMultiset


class TestPositiveMultiset(unittest.TestCase):
    def test_init(self):
        with self.subTest("accepts valid multisets"):
            PositiveMultiset(str, {"a": 2, "b": 5, "c": 3})
            PositiveMultiset(int, {1: 2, 2: 5, 3: 3}),

        with self.subTest("accepts infinite multisets"):
            PositiveMultiset(str, {"a": infinity, "b": 5, "c": 3}, allow_infinity=True)

        with self.subTest("allow_infinity flag set to False can disallow infinite quantities"):
            with self.assertRaises(AssertionError):
                PositiveMultiset(str, {"a": infinity, "b": 5, "c": 3}, allow_infinity=False)

        with self.subTest("all elements must be of the correct type"):
            with self.assertRaises(AssertionError):
                PositiveMultiset(str, {"a": 1, "b": 2, 5: 3})
            with self.assertRaises(AssertionError):
                PositiveMultiset(str, {"a": 1, 5: 2, "b": 3})
            with self.assertRaises(AssertionError):
                PositiveMultiset(str, {5: 1, "a": 2, "b": 3})
            with self.assertRaises(AssertionError):
                PositiveMultiset(int, {"a": 1, "b": 2, 5: 3})
            with self.assertRaises(AssertionError):
                PositiveMultiset(int, {"a": 1, 5: 2, "b": 3})
            with self.assertRaises(AssertionError):
                PositiveMultiset(int, {5: 1, "a": 2, "b": 3})

        with self.subTest("all counts must be positive integers or 'inf'"):
            for quantity in [0, -1, -2, 'a', '^', -infinity]:
                with self.subTest("Do not allow nonpositive quantities", quantity=quantity):
                    with self.assertRaises(AssertionError):
                        PositiveMultiset(str, {"x": 2, "y": quantity})
