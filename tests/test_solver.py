import unittest
from math import inf as infinity

from source.tbn import Tbn
from source.solver import Solver, SolverMethod, SolverFormulation
from source.constraints import Constraints


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.cp_solver = Solver(SolverMethod.CONSTRAINT_PROGRAMMING)
        self.ip_solver = Solver(SolverMethod.INTEGER_PROGRAMMING)

    def test_stable_config(self):
        test_cases = [
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("a* b* \n a b \n a* \n b*", 3, 1, self.ip_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("2[a* b*] \n a b", 2, 1, self.ip_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 2, 5, self.ip_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[a* b*] \n 2[a b]", infinity, 2, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[a* b*] \n 2[a b]", infinity, 2, self.ip_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[a* b*] \n 2[a b]", infinity, 2, self.ip_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.ip_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", infinity, 4, self.ip_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),
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
                    test_tbn, formulation=SolverFormulation.VARIABLE_BOND_WEIGHT, bond_weighting_factor=weight,
                )
                self.assertEqual(number_of_polymers, configuration.number_of_polymers())
                self.assertEqual(number_of_merges, configuration.number_of_merges())

    def test_stable_configs(self):
        test_cases = [
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("a* b* \n a b \n a* \n b*", 1, 1, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("a a \n a* a*", 2, 1, self.cp_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("a a \n a* a*", 1, 1, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_AWARE_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("2[a* b*] \n a b", 2, 1, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("2[a* b*] \n a b", 1, 1, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 5, self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 4, 5, self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", 3, 5, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[a* b*] \n 2[a b]", 1, 2, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[a* b*] \n 2[a b]", 1, 2, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),

            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 2, 4, self.cp_solver, SolverFormulation.POLYMER_UNBOUNDED_MATRIX),
            ("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]", 2, 4, self.cp_solver, SolverFormulation.VARIABLE_BOND_WEIGHT),
        ]
        for tbn_string, number_of_configs, number_of_merges, solver, formulation in test_cases:
            with self.subTest(tbn_string=tbn_string, solver=solver, formulation=formulation):
                test_tbn = Tbn.from_string(tbn_string)
                configurations = solver.stable_configs(test_tbn, formulation=formulation, bond_weighting_factor=2.0)
                self.assertEqual(number_of_configs, len(list(configurations)))
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
                configurations = list(solver.stable_configs(
                    test_tbn, formulation=SolverFormulation.VARIABLE_BOND_WEIGHT, bond_weighting_factor=weight,
                ))
                self.assertEqual(number_of_configs, len(configurations))
                for configuration in configurations:
                    self.assertEqual(energy, configuration.energy(weight))

    def test_configs_with_number_of_polymers(self):
        test_cases = [
            # second argument is a list of number of configurations expected for specific numbers of polymers:
            #  e.g. [a,b,c,d] => 'a' configurations with 1 polymer, 'b' configurations with 2 polymers, etc., up to 4

            ("a* b* \n a b \n a* \n b*", [ 1, 4, 1, 0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("a* b* \n a b \n a* \n b*", [ 1, 4, 1, 0], self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("a* b* \n a b \n a* \n b*", [ 1, 4, 1, 0], self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),

            ("2[a* b*] \n a b", [ 1, 2, 0, 0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("2[a* b*] \n a b", [ 1, 2, 0, 0], self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("2[a* b*] \n a b", [ 1, 1, 0, 0], self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),

            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [ 1, 4, 0, 0], self.cp_solver, SolverFormulation.BOND_OBLIVIOUS_NETWORK),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [ 1, 4, 0, 0], self.cp_solver, SolverFormulation.POLYMER_BINARY_MATRIX),
            ("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)", [ 1, 3, 0, 0], self.cp_solver, SolverFormulation.POLYMER_INTEGER_MATRIX),
        ]
        for tbn_string, number_of_configs_with_polymer_count, solver, formulation in test_cases:
            with self.subTest(tbn_string=tbn_string, solver=solver, formulation=formulation):
                test_tbn = Tbn.from_string(tbn_string)
                for polymer_count_minus_one, number_of_configs in enumerate(number_of_configs_with_polymer_count):
                    polymer_count = polymer_count_minus_one + 1
                    constraints = Constraints().with_fixed_polymers(polymer_count).with_unset_optimization_flag()
                    configurations = list(
                        solver.stable_configs(
                            test_tbn,
                            constraints,
                            formulation=formulation
                        )
                    )
                    self.assertEqual(number_of_configs, len(configurations))
