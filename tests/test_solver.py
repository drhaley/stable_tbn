import unittest
from math import inf as infinity

from source.tbn import Tbn
from source.solver import Solver, SolverMethod, SolverFormulation


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.cp_solver = Solver(SolverMethod.CONSTRAINT_PROGRAMMING)
        self.ip_solver = Solver(SolverMethod.INTEGER_PROGRAMMING)

    def test_stable_config(self):
        test_cases = [
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.SET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.LOW_W_FORMULATION),

            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.SET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.LOW_W_FORMULATION),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.SET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[a* b*] \n 2[a b]", infinity, 2, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[a* b*] \n 2[a b]", infinity, 2, self.ip_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.ip_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.ip_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.ip_solver, SolverFormulation.LOW_W_FORMULATION),
        ]
        for tbn_string, number_of_polymers, number_of_merges, solver, formulation in test_cases:
            with self.subTest(tbn_string=tbn_string, solver=solver, formulation=formulation):
                test_tbn = Tbn.from_string(tbn_string)
                configuration = solver.stable_config(test_tbn, formulation=formulation, bond_weighting_factor=2.0)
                self.assertEqual(number_of_polymers, configuration.number_of_polymers())
                self.assertEqual(number_of_merges, configuration.number_of_merges())

        low_w_test_cases = [
            ("a* b* \n a b \n a* \n b*", 4, 0, self.cp_solver, 0.4),
            ("a* b* \n a b \n a* \n b*", 4, 0, self.ip_solver, 0.4),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, 0.6),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, 0.6),

            ("2[a* b*] \n a b", 3, 0, self.cp_solver, 0.4),
            ("2[a* b*] \n a b", 3, 0, self.ip_solver, 0.4),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, 0.6),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, 0.6),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 5, 2, self.cp_solver, 0.4),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 5, 2, self.ip_solver, 0.4),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 3, self.cp_solver, 0.6),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 3, self.ip_solver, 0.6),

            ("inf[a* b*] \n 2[a b]", infinity, 0, self.cp_solver, 0.4),
            ("inf[a* b*] \n 2[a b]", infinity, 0, self.ip_solver, 0.4),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.cp_solver, 0.6),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.ip_solver, 0.6),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 2, self.cp_solver, 0.4),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 2, self.ip_solver, 0.4),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.cp_solver, 0.6),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.ip_solver, 0.6),
        ]

        for tbn_string, number_of_polymers, number_of_merges, solver, weight in low_w_test_cases:
            with self.subTest("low w tests", tbn_string=tbn_string, solver=solver, weight=weight):
                test_tbn = Tbn.from_string(tbn_string)
                configuration = solver.stable_config(
                    test_tbn, formulation=SolverFormulation.LOW_W_FORMULATION, bond_weighting_factor=weight,
                )
                self.assertEqual(number_of_polymers, configuration.number_of_polymers())
                self.assertEqual(number_of_merges, configuration.number_of_merges())

    def test_stable_configs(self):
        test_cases = [
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("a a \n a* a*", 2, 1, self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 5, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 5, self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[a* b*] \n 2[a b]", 1, 2, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[a* b*] \n 2[a b]", 1, 2, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 2, 4, self.cp_solver, SolverFormulation.BEYOND_MULTISET_FORMULATION),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 2, 4, self.cp_solver, SolverFormulation.LOW_W_FORMULATION),
        ]
        for tbn_string, number_of_configs, number_of_merges, solver, formulation in test_cases:
            with self.subTest(tbn_string=tbn_string, solver=solver, formulation=formulation):
                test_tbn = Tbn.from_string(tbn_string)
                configurations = solver.stable_configs(test_tbn, formulation=formulation, bond_weighting_factor=2.0)
                self.assertEqual(number_of_configs, len(configurations))
                for configuration in configurations:
                    self.assertEqual(number_of_merges, configuration.number_of_merges())

        low_w_test_cases = [
            ("a* b* \n a b \n a* \n b*", 1, 0.8, self.cp_solver, 0.4),
            ("a* b* \n a b \n a* \n b*", 1, 1.0, self.cp_solver, 0.6),

            ("2[a* b*] \n a b", 1, 0.8, self.cp_solver, 0.4),
            ("2[a* b*] \n a b", 1, 1.0, self.cp_solver, 0.6),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 1, 3.6, self.cp_solver, 0.4),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 1, 4.2, self.cp_solver, 0.6),

            ("inf[a* b*] \n 2[a b]", 1, 1.6, self.cp_solver, 0.4),
            ("inf[a* b*] \n 2[a b]", 1, 2.0, self.cp_solver, 0.6),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 1, 3.6, self.cp_solver, 0.4),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 2, 4.0, self.cp_solver, 0.6),
        ]

        for tbn_string, number_of_configs, energy, solver, weight in low_w_test_cases:
            with self.subTest("low w tests", tbn_string=tbn_string, solver=solver, weight=weight):
                test_tbn = Tbn.from_string(tbn_string)
                configurations = solver.stable_configs(
                    test_tbn, formulation=SolverFormulation.LOW_W_FORMULATION, bond_weighting_factor=weight,
                )
                self.assertEqual(number_of_configs, len(configurations))
                for configuration in configurations:
                    self.assertEqual(energy, configuration.energy(weight))

    def test_configs_with_number_of_polymers(self):
        test_cases = [
            # second argument is a list of number of configurations expected for specific numbers of polymers:
            #  e.g. [a,b,c,d] => 'a' configurations with 1 polymer, 'b' configurations with 2 polymers, etc., up to 4

            #("a* b* \n a b \n a* \n b*", [12, 7, 1, 0], self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("a* b* \n a b \n a* \n b*", [12, 7, 1, 0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("a* b* \n a b \n a* \n b*", [ 1, 4, 1, 0], self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("a* b* \n a b \n a* \n b*", [ 1, 4, 1, 0], self.cp_solver, SolverFormulation.MULTISET_FORMULATION),

            #("2[a* b*] \n a b", [ 5, 2, 0, 0], self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("2[a* b*] \n a b", [ 5, 2, 0, 0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("2[a* b*] \n a b", [ 1, 2, 0, 0], self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("2[a* b*] \n a b", [ 1, 1, 0, 0], self.cp_solver, SolverFormulation.MULTISET_FORMULATION),

            #("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [9,4,0,0], self.cp_solver, SolverFormulation.STABLEGEN_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [9,4,0,0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [1,4,0,0], self.cp_solver, SolverFormulation.SET_FORMULATION),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [1,3,0,0], self.cp_solver, SolverFormulation.MULTISET_FORMULATION),
        ]
        for tbn_string, number_of_configs_with_polymer_count, solver, formulation in test_cases:
            with self.subTest(tbn_string=tbn_string, solver=solver, formulation=formulation):
                test_tbn = Tbn.from_string(tbn_string)
                for polymer_count_minus_one, number_of_configs in enumerate(number_of_configs_with_polymer_count):
                    polymer_count = polymer_count_minus_one + 1
                    configurations = solver.configs_with_number_of_polymers(test_tbn, number_of_polymers=polymer_count,
                                                                            formulation=formulation)
                    self.assertEqual(number_of_configs, len(configurations))
