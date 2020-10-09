import unittest
from source.polymer import Polymer
from source.monomer import Monomer


class TestPolymer(unittest.TestCase):
    def setUp(self):
        self.x = Monomer.from_string("x0 x1", "X")
        self.y = Monomer.from_string("2(y0) 1(y1) 3(y2)", "Y")

        self.polymer_1x = Polymer({self.x: 1})
        self.polymer_1y = Polymer({self.y: 1})
        self.polymer_1x_1y = Polymer({self.y: 1, self.x: 1})
        self.polymer_2x_3y = Polymer({self.y: 3, self.x: 2})

    def test_init(self):
        with self.subTest("do not allow creation of an empty polymer"):
            with self.assertRaises(AssertionError):
                Polymer({})

        for quantity in [0, -1, -2, 'a', '^', float("inf")]:
            with self.subTest("Do not allow nonpositive or infinite monomer quantities", quantity=quantity):
                with self.assertRaises(AssertionError):
                    Polymer({self.x: 2, self.y: quantity})

    def test_size(self):
        tests = [
            (self.polymer_1x, 1),
            (self.polymer_1y, 1),
            (self.polymer_1x_1y, 2),
            (self.polymer_2x_3y, 5),
        ]
        for polymer, polymer_size in tests:
            with self.subTest(polymer_size=polymer_size):
                self.assertEqual(polymer_size, polymer.size())

    def test_str(self):
        tests = [
            (self.polymer_1x, "{X}"),
            (self.polymer_1y, "{Y}"),
            (self.polymer_1x_1y, "{X, Y}"),
            (self.polymer_2x_3y, "{2(X), 3(Y)}"),
        ]
        for polymer, polymer_as_string in tests:
            with self.subTest(polymer_as_string=polymer_as_string):
                self.assertEqual(polymer_as_string, str(polymer))

    def test_lt(self):
        self.assertTrue(self.polymer_2x_3y < self.polymer_1x_1y)
        self.assertTrue(self.polymer_1x_1y < self.polymer_1x)
        self.assertTrue(self.polymer_1x < self.polymer_1y)

        self.assertFalse(self.polymer_1x_1y < self.polymer_2x_3y)
        self.assertFalse(self.polymer_1x < self.polymer_1x_1y)
        self.assertFalse(self.polymer_1y < self.polymer_1x)
