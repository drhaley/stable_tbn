import unittest
from math import inf as infinity

from source.tbn import Tbn
from source.configuration import Configuration
from source.polymer import Polymer
from source.monomer import Monomer
from source.domain import Domain


class TestConfiguration(unittest.TestCase):
    def setUp(self):
        self.x = Monomer.from_string("x0 x1", "X")
        self.y = Monomer.from_string("2(y0) 1(y1) 3(y2)", "Y")

        self.polymer_1x = Polymer({self.x: 1})
        self.polymer_1y = Polymer({self.y: 1})
        self.polymer_1x_1y = Polymer({self.y: 1, self.x: 1})
        self.polymer_2x_3y = Polymer({self.y: 3, self.x: 2})
        self.polymer_1y_duplicate = Polymer({self.y: 1})

    def test_init(self):
        with self.subTest("Do not allow monomers in the configuration -- only polymers"):
            with self.assertRaises(AssertionError):
                Configuration({self.polymer_1x: 1, self.y: 1})
            with self.assertRaises(AssertionError):
                Configuration({self.x: 1, self.polymer_1y: 1})

        with self.subTest("Do not allow domains in the configuration -- only polymers"):
            with self.assertRaises(AssertionError):
                Configuration({self.polymer_1x: 1, Domain("y"): 1})
            with self.assertRaises(AssertionError):
                Configuration({Domain("x"): 1, self.polymer_1y: 1})

        with self.subTest("Allow multiple of same type of polymer"):
            test_configuration = Configuration({self.polymer_1x: 1, self.polymer_1y: 2})
            self.assertTrue(3, test_configuration.number_of_polymers())

        with self.subTest("Allow infinite quantities of polymers"):
            test_configuration = Configuration({self.polymer_1x: infinity, self.polymer_1y: 2})
            self.assertTrue(infinity, test_configuration.number_of_polymers())

    def test_str(self):
        test_configuration = Configuration({
            self.polymer_1y: 1, self.polymer_1x_1y: 2, self.polymer_2x_3y: 3, self.polymer_1x: 4
        })

        with self.subTest("able to display singletons in a configuration string"):
            intended_output_with_singletons = "3{2(X), 3(Y)}; 2{X, Y}; 4{X}; {Y}"
            self.assertEqual(intended_output_with_singletons, test_configuration.full_str())

        with self.subTest("able to suppress singleton output in a configuration string"):
            intended_output = "3{2(X), 3(Y)}; 2{X, Y}"
            self.assertEqual(intended_output, str(test_configuration))

        inf_test_configuration = Configuration({
            self.polymer_1y: infinity, self.polymer_1x_1y: infinity, self.polymer_2x_3y: 3, self.polymer_1x: 4
        })

        with self.subTest("able to display singletons in a configuration string"):
            intended_output_with_singletons = "3{2(X), 3(Y)}; inf{X, Y}; 4{X}; inf{Y}"
            self.assertEqual(intended_output_with_singletons, inf_test_configuration.full_str())

        with self.subTest("able to suppress singleton output in a configuration string"):
            intended_output = "3{2(X), 3(Y)}; inf{X, Y}"
            self.assertEqual(intended_output, str(inf_test_configuration))

    def test_eq(self):
        test_configurations = []
        for i in range(2):
            test_configurations.append(
                Configuration({self.polymer_1y: 1, self.polymer_1x_1y: 2, self.polymer_2x_3y: 3, self.polymer_1x: 4})
            )
        self.assertEqual(test_configurations[0], test_configurations[0])
        self.assertEqual(test_configurations[1], test_configurations[1])
        self.assertEqual(test_configurations[0], test_configurations[1])
        self.assertNotEqual(test_configurations[0], Configuration({}))
        self.assertNotEqual(Configuration({}), test_configurations[0])

    def test_size(self):
        test_polymers = [
            self.polymer_1x,
            self.polymer_1y,
            self.polymer_1x_1y,
            self.polymer_2x_3y,
        ]
        with self.subTest("Sums over single copies of different polymers"):
            self.assertEqual(len(test_polymers), Configuration({polymer: 1 for polymer in test_polymers}).number_of_polymers())
        with self.subTest("Sums over multiple copies of the same polymer(s)"):
            for i in range(1, 4):
                for j in range(1, 4):
                    self.assertEqual(i+j, Configuration({test_polymers[0]: i, test_polymers[1]: j}).number_of_polymers())
        with self.subTest("Reports 'inf' if any polymer has infinite quantity"):
            self.assertEqual(
                infinity,
                Configuration({test_polymers[2]: 3, test_polymers[3]: infinity}).number_of_polymers()
            )
            self.assertEqual(
                infinity,
                Configuration({test_polymers[2]: infinity, test_polymers[3]: infinity}).number_of_polymers()
            )

    def test_number_of_merges(self):
        self.assertEqual(0, Configuration({}).number_of_merges())
        self.assertEqual(1, Configuration({self.polymer_1x_1y: 1}).number_of_merges())
        self.assertEqual(4, Configuration({self.polymer_2x_3y: 1}).number_of_merges())
        self.assertEqual(5, Configuration({self.polymer_1x_1y: 1, self.polymer_2x_3y: 1}).number_of_merges())
        self.assertEqual(9, Configuration({self.polymer_1x_1y: 1, self.polymer_2x_3y: 2}).number_of_merges())
        self.assertEqual(10, Configuration({self.polymer_1x_1y: 2, self.polymer_2x_3y: 2}).number_of_merges())

        self.assertEqual(0, Configuration({self.polymer_1x: 1, self.polymer_1y: infinity}).number_of_merges())
        self.assertEqual(1, Configuration({self.polymer_1x_1y: 1, self.polymer_1y: infinity}).number_of_merges())
        self.assertEqual(infinity, Configuration({self.polymer_1x_1y: infinity, self.polymer_1y: 1}).number_of_merges())

    def test_energy(self):
        # values determined by regression testing
        monomer_A = Monomer.from_string("a a b b c d")
        monomer_B = Monomer.from_string("a* b* c*")
        polymer_A = Polymer({monomer_A: 1})
        polymer_B = Polymer({monomer_B: 1})
        polymer_AB = Polymer({monomer_A: 1, monomer_B: 1})
        polymer_AAB = Polymer({monomer_A: 2, monomer_B: 1})
        polymer_ABB = Polymer({monomer_A: 1, monomer_B: 2})
        self.assertEqual(0.0, Configuration({}).energy(bond_weighting_factor=0.4))
        self.assertEqual(1.0, Configuration({polymer_A: 1, polymer_AB: 1}).energy(bond_weighting_factor=0.4))
        self.assertEqual(2.2, Configuration({polymer_A: infinity, polymer_B: 1, polymer_AB: 1}).energy(bond_weighting_factor=0.4))
        self.assertEqual(1.8, Configuration({polymer_B: infinity, polymer_AB: 1}).energy(bond_weighting_factor=0.4))
        self.assertEqual(4.8, Configuration({polymer_B: infinity, polymer_AAB: 1}).energy(bond_weighting_factor=0.4))
        self.assertEqual(2.4, Configuration({polymer_A: infinity, polymer_ABB: 1}).energy(bond_weighting_factor=0.4))

        # case to test a fix to a bug in a specific example in which energy was misreported as 0.8 (should be 1.6)
        #  two limiting monomers with two domains each, all as singletons => 0.4(4 domains) = 1.6
        monomer_X = Monomer.from_string("a b")
        monomer_Y = Monomer.from_string("a* b*")
        polymer_X = Polymer({monomer_X: 1})
        polymer_Y = Polymer({monomer_Y: 1})
        self.assertEqual(1.6, Configuration({polymer_X: 2, polymer_Y: infinity}).energy(bond_weighting_factor=0.4))

    def test_flatten(self):
        self.assertEqual(
            Tbn({}),
            Configuration({}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 1, self.y: 1}),
            Configuration({self.polymer_1x_1y: 1}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 2, self.y: 3}),
            Configuration({self.polymer_2x_3y: 1}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 3, self.y: 4}),
            Configuration({self.polymer_1x_1y: 1, self.polymer_2x_3y: 1}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 5, self.y: 7}),
            Configuration({self.polymer_1x_1y: 1, self.polymer_2x_3y: 2}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 6, self.y: 8}),
            Configuration({self.polymer_1x_1y: 2, self.polymer_2x_3y: 2}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 1, self.y: infinity}),
            Configuration({self.polymer_1x: 1, self.polymer_1y: infinity}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: 1, self.y: infinity}),
            Configuration({self.polymer_1x_1y: 1, self.polymer_1y: infinity}).flatten()
        )
        self.assertEqual(
            Tbn({self.x: infinity, self.y: infinity}),
            Configuration({self.polymer_1x_1y: infinity, self.polymer_1y: 1}).flatten()
        )
