import unittest
from typing import Dict
from source.monomer import Monomer
from source.domain import Domain


class TestMonomer(unittest.TestCase):
    def setUp(self):
        self.x = Monomer.from_string("x0 x1 >X")
        self.y = Monomer.from_string("2(y0) 1(y1) 3(y2) >Y")
        self.abc_star = Monomer.from_string("a b c*")
        self.triple_a_star = Monomer.from_string("a* a* a*")
        self.semi_self_saturated_monomer = Monomer.from_string("a a* b b b*")

    def test_init(self):
        for quantity in [0, -1, -2, 'a', '^', float("inf")]:
            with self.subTest("Do not allow nonpositive or infinite domain quantities", quantity=quantity):
                with self.assertRaises(AssertionError):
                    Monomer({Domain("xx"): 2, Domain("yy"): quantity}, "bad_monomer")

        with self.subTest("Do not allow empty monomers"):
            with self.assertRaises(AssertionError):
                Monomer.from_string("")
            with self.assertRaises(AssertionError):
                Monomer.from_string(">empty_monomer")

        for name in ["X", "Y"]:
            with self.subTest("No duplicate names", name=name):
                with self.assertRaises(AssertionError):
                    Monomer.from_string(f"z0 z1 >{name}")

        with self.subTest("Accept duplicate monomer name if domains match"):
            Monomer.from_string("x0 x1 >X")
            Monomer.from_string("2(y0) 1(y1) 3(y2) >Y")
            Monomer.from_string("2(y0) y1 3(y2) >Y")

        with self.subTest("Accept duplicate monomer name if domains match but are switched"):
            Monomer.from_string("x1 x0 >X")
            Monomer.from_string("1(y1) 2(y0) 3(y2) >Y")
            Monomer.from_string("y1 2(y0) 3(y2) >Y")
            Monomer.from_string("c* a b")
            Monomer.from_string("a b a* b b*")

        with self.subTest("Do not collapse complements; e.g. do not annihilate a with a*"):
            Monomer.from_string("b")
            Monomer.from_string("a a* b")

    def test_str(self):
        tests = [
            (self.x, "X"),
            (self.y, "Y"),
            (self.abc_star, "[a b c*]"),
            (self.triple_a_star, "[3(a*)]"),
            (self.semi_self_saturated_monomer, "[a a* 2(b) b*]"),
        ]
        for monomer, monomer_as_string in tests:
            with self.subTest(monomer_as_string=monomer_as_string):
                self.assertEqual(monomer_as_string, str(monomer))

    def test_eq(self):
        with self.subTest("equal to self"):
            self.assertEqual(self.x, self.x)
            self.assertEqual(self.y, self.y)
            self.assertEqual(self.abc_star, self.abc_star)
            self.assertEqual(self.triple_a_star, self.triple_a_star)
            self.assertEqual(self.semi_self_saturated_monomer, self.semi_self_saturated_monomer)

        with self.subTest("equal to equivalent monomer"):
            self.assertEqual(self.x, Monomer.from_string("x0 x1 >X"))
            self.assertEqual(self.y, Monomer.from_string("2(y0) 1(y1) 3(y2) >Y"))
            self.assertEqual(self.abc_star, Monomer.from_string("a b c*"))
            self.assertEqual(self.triple_a_star, Monomer.from_string("a* a* a*"))
            self.assertEqual(self.triple_a_star, Monomer.from_string("3(a*)"))
            self.assertEqual(self.semi_self_saturated_monomer, Monomer.from_string("a a* b b b*"))

        with self.subTest("equal to equivalent monomer with order switched"):
            self.assertEqual(self.x, Monomer.from_string("x1 x0 >X"))
            self.assertEqual(self.abc_star, Monomer.from_string("c* a b"))
            self.assertEqual(self.y, Monomer.from_string("y1 2(y0) 3(y2) >Y"))
            self.assertEqual(self.semi_self_saturated_monomer, Monomer.from_string("a b a* b b*"))

        with self.subTest("not equal to monomer with different name"):
            self.assertNotEqual(self.x, Monomer.from_string("z0 z1 >Z_test1"))
            self.assertNotEqual(self.x, Monomer.from_string("x0 x1 >Z_test2"))
            self.assertNotEqual(self.y, Monomer.from_string("2(y0) 1(y1) 3(y2) >Z_test3"))
            self.assertNotEqual(self.triple_a_star, Monomer.from_string("a* a* a* >Z_test4"))

    def test_lt(self):
        self.assertTrue((self.x < self.y) or (self.y < self.x))
        self.assertFalse((self.x < self.y) and (self.y < self.x))

    def test_hash(self):
        self.assertNotEqual(hash(self.x), hash(self.y))

    def test_from_string(self):
        class SensingMonomer(Monomer):
            def __init__(self, sensed_domain_multiset: Dict[Domain, int], sensed_name: str):
                super().__init__(sensed_domain_multiset, sensed_name)
                self.sensed_domain_multiset = sensed_domain_multiset
                self.sensed_name = sensed_name

        a = Domain("a")
        a_star = a.complement()
        b = Domain("b")
        b_star = b.complement()
        tests = [
            ({a: 2, b: 3, b_star: 1}, "aabbb_star"),
            ({a: 1, b_star: 1}, "ab_star"),
            ({a: 3, b_star: 1}, "aaab_star"),
        ]
        for domain_multiset, name in tests:
            domain_multiset_as_string = " ".join([
                " ".join(count * [str(domain)])
                    for domain, count in domain_multiset.items()
            ])
            with self.subTest("from_string() correctly parses and sends to __init__()"):
                created_monomer = SensingMonomer.from_string(f"{domain_multiset_as_string} >{name}")
                self.assertEqual(domain_multiset, created_monomer.sensed_domain_multiset)
                self.assertEqual(name, created_monomer.sensed_name)

            with self.subTest("from_string() correctly defaults to 'None' when no name is provided"):
                unnamed_monomer = SensingMonomer.from_string(domain_multiset_as_string)
                self.assertEqual(None, unnamed_monomer.sensed_name)

        with self.subTest("correctly collapses to multiset"):
            new_monomer = SensingMonomer.from_string("a a a b a*")
            self.assertEqual({a: 3, a_star: 1, b: 1}, new_monomer.sensed_domain_multiset)

        with self.subTest("ignores legacy domain names appearing after colon"):
            legacy_monomer = SensingMonomer.from_string("a a:name a* b:name2")
            self.assertEqual({a: 2, a_star: 1, b: 1}, legacy_monomer.sensed_domain_multiset)

    def test_name(self):
        tests = [
            (self.x, "X"),
            (self.abc_star, "[a b c*]"),
        ]
        for monomer, monomer_name in tests:
            self.assertEqual(monomer_name, monomer.name())

    def test_unstarred_domain_types(self):
        tests = [
            (self.x, [Domain("x0"), Domain("x1")]),
            (self.y, [Domain("y0"), Domain("y1"), Domain("y2")]),
            (self.abc_star, [Domain("a"), Domain("b"), Domain("c")]),
            (self.triple_a_star, [Domain("a")]),
            (self.semi_self_saturated_monomer, [Domain("a"), Domain("b")]),
        ]
        for monomer, unstarred_domain_types in tests:
            with self.subTest("must return unstarred domain types", monomer=str(monomer)):
                self.assertEqual(unstarred_domain_types, list(monomer.unstarred_domain_types()))

    def test_net_count(self):
        tests = [
            (self.x, "x0", 1),
            (self.x, "x1", 1),
            (self.x, "x2", 0),
            (self.x, "x0*", -1),
            (self.x, "x1*", -1),
            (self.x, "x2*", 0),
            (self.y, "y0", 2),
            (self.y, "y1", 1),
            (self.y, "y2", 3),
            (self.y, "y0*", -2),
            (self.y, "y1*", -1),
            (self.y, "y2*", -3),
            (self.abc_star, "a", 1),
            (self.abc_star, "b", 1),
            (self.abc_star, "c", -1),
            (self.triple_a_star, "a", -3),
            (self.triple_a_star, "a*", 3),
            (self.semi_self_saturated_monomer, "a", 0),
            (self.semi_self_saturated_monomer, "b", 1),
            (self.semi_self_saturated_monomer, "a*", 0),
            (self.semi_self_saturated_monomer, "b*", -1),
        ]
        for monomer, domain_string, net_count in tests:
            with self.subTest("net count", monomer=str(monomer), domain=domain_string, net_count=net_count):
                self.assertEqual(net_count, monomer.net_count(Domain(domain_string)))

    def test_regex(self):
        multiple_domain_regex = f"(?:{Domain.regex()}|[1-9]\\d*\\(\\s*{Domain.regex()}\\s*\\))"
        name_regex = r"[A-Za-z0-9_]+"
        domain_list_regex = f"{multiple_domain_regex}(?: {multiple_domain_regex})*"
        self.assertEqual(Monomer.regex(), f"{domain_list_regex}(?:|\\s*>{name_regex})")

    def test_as_explicit_list(self):
        self.assertEqual(
            [Domain("x0"), Domain("x1")],
            self.x.as_explicit_list()
        )
        self.assertEqual(
            [Domain("y0"), Domain("y0"), Domain("y1"), Domain("y2"), Domain("y2"), Domain("y2")],
            self.y.as_explicit_list()
        )
        self.assertEqual(
            [Domain("a"), Domain("b"), Domain("c*")],
            self.abc_star.as_explicit_list()
        )
        self.assertEqual(
            [Domain("a*"), Domain("a*"), Domain("a*")],
            self.triple_a_star.as_explicit_list()
        )
        self.assertEqual(
            [Domain("a"), Domain("a*"), Domain("b"), Domain("b"), Domain("b*")],
            self.semi_self_saturated_monomer.as_explicit_list()
        )
