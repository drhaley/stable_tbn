import unittest
from source.domain import Domain


class TestDomain(unittest.TestCase):
    def setUp(self):
        self.a = Domain("a")
        self.a_star = Domain("a*")

    def test_init(self):
        with self.subTest("do not accept empty string as a name"):
            with self.assertRaises(AssertionError):
                Domain("")
            with self.assertRaises(AssertionError):
                Domain("*")

        bad_star_strings = [
            "a**",
            "a*a",
            "x0*0*",
            "*00*",
            "*a",
        ]
        for bad_star_string in bad_star_strings:
            with self.subTest("do not accept stars except at end of domain name", bad_star_string=bad_star_string):
                with self.assertRaises(AssertionError):
                    Domain(bad_star_string)

        bad_parenthesis_strings = [
            "a(",
            "a)",
            "a()",
            "a(a",
            "a)a",
            "x0(0)",
            "x0)0(",
            "(00)",
            "(a",
            ")a",
        ]
        for bad_parenthesis_string in bad_parenthesis_strings:
            with self.subTest("no parenthesis in domain names", bad_parenthesis_string=bad_parenthesis_string):
                with self.assertRaises(AssertionError):
                    Domain(bad_parenthesis_string)

        with self.subTest("cannot contain spaces"):
            with self.assertRaises(AssertionError):
                Domain("my domain")

    def test_str(self):
        tests = [
            (self.a, "a"),
            (self.a_star, "a*"),
        ]
        for domain, domain_as_string in tests:
            with self.subTest(domain_as_string=domain_as_string):
                self.assertEqual(domain_as_string, str(domain))

    def test_eq(self):
        with self.subTest("equal to self"):
            self.assertEqual(self.a, self.a)
            self.assertEqual(self.a_star, self.a_star)

        with self.subTest("equal to equivalent domain"):
            self.assertEqual(self.a, Domain("a"))
            self.assertEqual(self.a_star, Domain("a*"))

        with self.subTest("not equal to complement"):
            self.assertNotEqual(self.a, self.a_star)

        with self.subTest("not equal to unlike domain"):
            self.assertNotEqual(self.a, Domain("x"))
            self.assertNotEqual(self.a_star, Domain("x*"))

    def test_lt(self):
        self.assertTrue(Domain("b") < Domain("c"))
        self.assertFalse(Domain("c") < Domain("b"))

        self.assertTrue(Domain("b*") < Domain("c"))
        self.assertFalse(Domain("c") < Domain("b*"))

    def test_hash(self):
        self.assertNotEqual(hash(self.a), hash(self.a_star))

    def test_is_starred(self):
        tests = [
            (self.a, False),
            (self.a_star, True),
        ]
        for domain, output in tests:
            with self.subTest(domain=str(domain)):
                self.assertEqual(output, domain.is_starred())

    def test_complement(self):
        with self.subTest("starring an unstarred domain"):
            self.assertEqual(self.a_star, self.a.complement())

        with self.subTest("unstarring a starred domain"):
            self.assertEqual(self.a, self.a_star.complement())

        for domain in [self.a, self.a_star]:
            with self.subTest("involution test -- complementing twice", domain=str(domain)):
                self.assertEqual(domain, domain.complement().complement())

    def test_regex(self):
        self.assertEqual(Domain.regex(), r"[A-Za-z0-9_]+(?:\*|)(?:\:[A-Za-z0-9_]+|)")
