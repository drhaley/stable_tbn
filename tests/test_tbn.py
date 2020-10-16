import unittest
from math import inf as infinity

from source.tbn import Tbn
from source.monomer import Monomer
from source.domain import Domain


class TestTbn(unittest.TestCase):
    def setUp(self):
        self.x = Monomer.from_string("x0 x1", "X")
        self.y = Monomer.from_string("2(y0) 1(y1) 3(y2)", "Y")

        self.Tbn_1x = Tbn({self.x: 1})
        self.Tbn_1y = Tbn({self.y: 1})
        self.Tbn_1x_1y = Tbn({self.y: 1, self.x: 1})
        self.Tbn_2x_3y = Tbn({self.y: 3, self.x: 2})
        self.Tbn_infx_2y = Tbn({self.x: infinity, self.y: 2})
        self.Tbn_2x_infy = Tbn({self.y: infinity, self.x: 2})
        self.Tbn_1x_infy = Tbn({self.y: infinity, self.x: 1})

    def test_init(self):
        for quantity in [0, -1, -2, 'a', '^', '-inf']:
            with self.subTest("Do not allow nonpositive monomer quantities", quantity=quantity):
                with self.assertRaises(AssertionError):
                    Tbn({self.x: 2, self.y: quantity})

        with self.subTest("allow infinite monomer quantities"):
            Tbn({self.x: infinity, self.y: 2})
            Tbn({self.x: 2, self.y: infinity})

    def test_str(self):
        tests = [
            (self.Tbn_1x, "{X}"),
            (self.Tbn_1y, "{Y}"),
            (self.Tbn_1x_1y, "{X, Y}"),
            (self.Tbn_2x_3y, "{2(X), 3(Y)}"),
            (self.Tbn_infx_2y, "{inf(X), 2(Y)}"),
            (self.Tbn_2x_infy, "{2(X), inf(Y)}"),
            (self.Tbn_1x_infy, "{X, inf(Y)}"),
        ]
        for tbn, Tbn_as_string in tests:
            with self.subTest(Tbn_as_string=Tbn_as_string):
                self.assertEqual(Tbn_as_string, str(tbn))

    def test_lt(self):
        self.assertTrue(self.Tbn_2x_3y < self.Tbn_1x_1y)
        self.assertTrue(self.Tbn_1x_1y < self.Tbn_1x)
        self.assertTrue(self.Tbn_1x < self.Tbn_1y)
        self.assertTrue(self.Tbn_2x_infy < self.Tbn_2x_3y)
        self.assertTrue(self.Tbn_2x_3y < self.Tbn_1x_infy)
        self.assertTrue(self.Tbn_infx_2y < self.Tbn_1x_infy)

        self.assertFalse(self.Tbn_1x_1y < self.Tbn_2x_3y)
        self.assertFalse(self.Tbn_1x < self.Tbn_1x_1y)
        self.assertFalse(self.Tbn_1y < self.Tbn_1x)
        self.assertFalse(self.Tbn_2x_3y < self.Tbn_2x_infy)
        self.assertFalse(self.Tbn_1x_infy < self.Tbn_2x_3y)
        self.assertFalse(self.Tbn_1x_infy < self.Tbn_infx_2y)

    def test_eq(self):
        self.assertEqual(Tbn({Monomer.from_string("a"): 1}), Tbn({Monomer.from_string("a"): 1}))
        self.assertEqual(Tbn({Monomer.from_string("a"): infinity}), Tbn({Monomer.from_string("a"): infinity}))
        self.assertNotEqual(Tbn({Monomer.from_string("a"): 1}), Tbn({Monomer.from_string("a"): 2}))
        self.assertNotEqual(Tbn({Monomer.from_string("a"): 1}), Tbn({Monomer.from_string("b"): 1}))
        self.assertNotEqual(Tbn({Monomer.from_string("a"): 1}), Tbn({Monomer.from_string("a"): infinity}))

    def test_from_string(self):
        with self.subTest("single monomer example"):
            self.assertEqual(Tbn({Monomer.from_string("a"): 1}), Tbn.from_string("a"))

        with self.subTest("single monomer type, multiple monomer example"):
            self.assertEqual(Tbn({Monomer.from_string("a"): 2}), Tbn.from_string("2[a]"))

        with self.subTest("testing from_string with example from stablegen.net/help"):
            example_text = \
                """
                a*:b1 b*
                a b:b2 >m1
                a* >m2
                b*
                """
            example_tbn = Tbn({
                Monomer.from_string("a* b*"): 1,
                Monomer.from_string("a b", "m1"): 1,
                Monomer.from_string("a*", "m2"): 1,
                Monomer.from_string("b*"): 1,
            })
            self.assertEqual(example_tbn, Tbn.from_string(example_text))

        with self.subTest("testing from_string with multisets"):
            multiset_text = \
                """
                2[ 3(a*) a b ]
                [ c ]
                5[ a a:favorite_a b >bob ]
                7[ b a b b ]
                b a b b b*
                """
            multiset_tbn = Tbn({
                Monomer.from_string("a* a* a* a b"): 2,
                Monomer.from_string("c"): 1,
                Monomer.from_string("2(a) b", "bob"): 5,
                Monomer.from_string("3(b) a"): 7,
                Monomer.from_string("b* 3(b) a"): 1,
            })
            self.assertEqual(multiset_tbn, Tbn.from_string(multiset_text))

        with self.subTest("testing from_string with excess monomers"):
            multiset_text = \
                """
                2[ 3(a*) a b ]
                inf[ c ]
                5[ a a b >bob ]
                inf[ b a b b ]
                b a b b b*
                """
            multiset_tbn = Tbn({
                Monomer.from_string("a* a* a* a b"): 2,
                Monomer.from_string("c"): infinity,
                Monomer.from_string("2(a) b", "bob"): 5,
                Monomer.from_string("3(b) a"): infinity,
                Monomer.from_string("b* 3(b) a"): 1,
            })
            self.assertEqual(multiset_tbn, Tbn.from_string(multiset_text))

    def test_monomer_types(self):
        tests = [
            ({}, []),
            ({self.x: 3}, [self.x]),
            ({self.y: 5}, [self.y]),
            ({self.x: 5, self.y: 2}, [self.x, self.y]),
            ({self.x: infinity, self.y: infinity}, [self.x, self.y]),
        ]
        for monomer_multiset, monomer_types in tests:
            tbn = Tbn(monomer_multiset)
            with self.subTest("ordinary monomer type iterator", tbn=tbn):
                self.assertEqual(monomer_types, list(tbn.monomer_types()))
        flatten_tests = [
            ({}, []),
            ({self.x: 3}, [self.x]),
            ({self.y: 5}, [self.y]),
            ({self.x: 5, self.y: 2}, [self.x, self.y]),
        ]
        for monomer_multiset, monomer_types in flatten_tests:
            tbn = Tbn(monomer_multiset)
            with self.subTest("monomer type iterator with flatten", tbn=tbn):
                flattened_list = list(tbn.monomer_types(flatten=True))
                for monomer_type in monomer_types:
                    self.assertEqual(tbn.count(monomer_type), flattened_list.count(monomer_type))


    def test_limiting_domain_types(self):
        tests = [
            ({}, []),
            ({self.x: 3}, [Domain("x0*"), Domain("x1*")]),
            ({self.y: 5}, [Domain("y0*"), Domain("y1*"), Domain("y2*")]),
            ({self.x: 5, self.y: 2}, [Domain("x0*"), Domain("x1*"), Domain("y0*"), Domain("y1*"), Domain("y2*")]),
            ({Monomer.from_string("a*"): 2, Monomer.from_string("a"): 1}, [Domain("a")]),
            ({Monomer.from_string("3(a*)"): 1, Monomer.from_string("a"): 2}, [Domain("a")]),
            ({Monomer.from_string("a*"): 1, Monomer.from_string("a"): 2}, [Domain("a*")]),
            ({Monomer.from_string("2(a*)"): 1, Monomer.from_string("a"): 2}, [Domain("a*")]),
            ({Monomer.from_string("2(a*)"): 1, Monomer.from_string("a"): infinity}, [Domain("a*")]),
            ({Monomer.from_string("2(a*)"): infinity, Monomer.from_string("a"): 2}, [Domain("a")]),
        ]
        for monomer_multiset, expected_limiting_domain_types in tests:
            with self.subTest("limiting domain types", tbn=str(Tbn(monomer_multiset))):
                limiting_domain_types = list(Tbn(monomer_multiset).limiting_domain_types())
                self.assertEqual(expected_limiting_domain_types, limiting_domain_types)

        with self.subTest("cannot have conflicting excess domain types"):
            conflicting_excess_tbn = Tbn(
                {Monomer.from_string("a"): infinity, Monomer.from_string("a*"): infinity}
            )
            with self.assertRaises(AssertionError):
                list(conflicting_excess_tbn.limiting_domain_types())

    def test_limiting_monomer_types(self):
        test_tbn = Tbn({
            Monomer.from_string("a b c e"): infinity,
            Monomer.from_string("a d*"): 1,
            Monomer.from_string("d*"): infinity,
            Monomer.from_string("a*"): 2,  # a* is limiting in this example
            Monomer.from_string("a d e*"): 3,  # both d and e* are limiting in this example
            Monomer.from_string("f"): 1,
            Monomer.from_string("f*"): 1,  # f* is chosen as limiting to break the tie (favors stars)
        })
        limiting_monomer_types = sorted([
            Monomer.from_string("a d e*"),
            Monomer.from_string("a*"),
            Monomer.from_string("f*")
        ])
        self.assertEqual(limiting_monomer_types, list(test_tbn.limiting_monomer_types()))

    def test_count(self):
        tests = [
            (self.Tbn_1x, self.x, 1),
            (self.Tbn_1x, self.y, 0),
            (self.Tbn_1y, self.x, 0),
            (self.Tbn_1y, self.y, 1),
            (self.Tbn_1x_1y, self.x, 1),
            (self.Tbn_1x_1y, self.y, 1),
            (self.Tbn_2x_3y, self.x, 2),
            (self.Tbn_2x_3y, self.y, 3),
            (self.Tbn_2x_infy, self.x, 2),
            (self.Tbn_2x_infy, self.y, infinity),
            (self.Tbn_infx_2y, self.x, infinity),
            (self.Tbn_infx_2y, self.y, 2),
        ]
        for tbn, monomer, count in tests:
            with self.subTest("correct monomer counts from tbn", tbn=str(tbn), monomer=str(monomer), count=count):
                self.assertEqual(count, tbn.count(monomer))

    def test_subtract(self):
        tests = [
            (self.Tbn_1x, self.Tbn_1x_1y - self.Tbn_1y),
            (self.Tbn_1y, self.Tbn_2x_3y - self.Tbn_1x - self.Tbn_1x - self.Tbn_1y - self.Tbn_1y),
            (self.Tbn_2x_infy, self.Tbn_2x_infy - self.Tbn_1y),
        ]
        for first, second in tests:
            with self.subTest("subtracton subtest", first=str(first), second=str(second)):
                self.assertEqual(first, second)

        with self.subTest("do not allow subtraction unless it is a subset of the first"):
            with self.assertRaises(AssertionError):
                self.Tbn_1x - self.Tbn_1y
            with self.assertRaises(AssertionError):
                self.Tbn_1x_1y - self.Tbn_2x_3y
