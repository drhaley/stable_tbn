import unittest
from math import inf as infinity

from source.tbn import Tbn
from source.solver import Solver, SolverMethod, SolverFormulation

MultisetFormulation = {x for x in SolverFormulation
                            if x not in
                            {SolverFormulation.STABLE_GEN_FORMULATION, SolverFormulation.BOND_OBLIVIOUS_FORMULATION}
                        }
InfinityFormulation = {SolverFormulation.BEYOND_MULTISET_FORMULATION}


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.single_solvers = []
        self.complete_solvers = []
        for solver_method in SolverMethod:
            self.single_solvers.append((
                Solver(method=solver_method),
                f"SOLV:{solver_method}",
            ))
            if solver_method != SolverMethod.INTEGER_PROGRAMMING:
                self.complete_solvers.append((
                    Solver(method=solver_method),
                    f"SOLV:{solver_method}",
                ))

    def test_stable_config(self):
        for solver, parameter_string in self.single_solvers:
            self.__run_test_stable_config_with_specific_solver(solver, parameter_string)

    def test_stable_configs(self):
        for solver, parameter_string in self.complete_solvers:
            self.__run_test_stable_configs_with_specific_solver(solver, parameter_string)

    def test_configs_with_number_of_polymers(self):
        for solver, parameter_string in self.complete_solvers:
            self.__run_test_configs_with_number_of_polymers_with_specific_solver(solver, parameter_string)

    def __run_test_stable_config_with_specific_solver(self, solver: Solver, parameter_string: str):
        for formulation in SolverFormulation:
            with self.subTest("check for correctness in simple example", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("a* b* \n a b \n a* \n b*")
                config = solver.stable_config(test_tbn, formulation=formulation)
                self.assertEqual(3, config.number_of_polymers())
                self.assertEqual(1, config.number_of_merges())

        for formulation in MultisetFormulation:
            with self.subTest("check for no combinatorial duplicates", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("2[a* b*] \n a b")
                config = solver.stable_config(test_tbn, formulation=formulation)
                self.assertEqual(2, config.number_of_polymers())
                self.assertEqual(1, config.number_of_merges())

            with self.subTest("check for return of multiple stable configurations", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)")
                config = solver.stable_config(test_tbn, formulation=formulation)
                self.assertEqual(2, config.number_of_polymers())
                self.assertEqual(5, config.number_of_merges())

        for formulation in InfinityFormulation:
            with self.subTest("check for correctness with infinite monomers", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("inf[a* b*] \n 2[a b]")
                config = solver.stable_config(test_tbn, formulation=formulation)
                self.assertEqual(infinity, config.number_of_polymers())
                self.assertEqual(2, config.number_of_merges())

            with self.subTest("second check for correctness with infinite monomers", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]")
                config = solver.stable_config(test_tbn, formulation=formulation)
                self.assertEqual(infinity, config.number_of_polymers())
                self.assertEqual(4, config.number_of_merges())

    def __run_test_stable_configs_with_specific_solver(self, solver: Solver, parameter_string: str):
        for formulation in SolverFormulation:
            with self.subTest("check for correctness in simple example", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("a* b* \n a b \n a* \n b*")
                results = list(solver.stable_configs(test_tbn, formulation=formulation))
                self.assertEqual(1, len(results))
                self.assertEqual(1, results[0].number_of_merges())

        for formulation in MultisetFormulation:
            with self.subTest("check for no combinatorial duplicates", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("2[a* b*] \n a b")
                results = list(solver.stable_configs(test_tbn, formulation=formulation))
                self.assertEqual(1, len(results))
                self.assertEqual(1, results[0].number_of_merges())

            with self.subTest("check for return of multiple stable configurations", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)")
                results = list(solver.stable_configs(test_tbn, formulation=formulation))
                self.assertEqual(3, len(results))  # (6 + 3 - 5 - 4), (6 - 1 - 5), (6 - 2 - 4)
                self.assertEqual(5, results[0].number_of_merges())
                self.assertEqual(5, results[1].number_of_merges())
                self.assertEqual(5, results[2].number_of_merges())

        for formulation in InfinityFormulation:
            with self.subTest("check for correctness with infinite monomers", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("inf[a* b*] \n 2[a b]")
                results = list(solver.stable_configs(test_tbn, formulation=formulation))
                self.assertEqual(1, len(results))
                self.assertEqual(2, results[0].number_of_merges())

            with self.subTest("second check for correctness with infinite monomers", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("inf[2(a*) 2(b*)] \n 2[3(a) 3(b)]")
                results = list(solver.stable_configs(test_tbn, formulation=formulation))
                self.assertEqual(2, len(results))
                self.assertEqual(4, results[0].number_of_merges())

    def __run_test_configs_with_number_of_polymers_with_specific_solver(self, solver: Solver, parameter_string: str):
        for formulation in set(SolverFormulation).difference(InfinityFormulation):
            with self.subTest("check for correctness in simple example", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("a* b* \n a b \n a* \n b*")
                results = {
                    i: list(solver.configs_with_number_of_polymers(test_tbn, number_of_polymers=i, formulation=formulation))
                        for i in range(1, 5)
                }
                expected_number_of_configurations = {1: 1, 2: 4, 3: 1, 4: 0}
                for i in expected_number_of_configurations.keys():
                    self.assertEqual(expected_number_of_configurations[i], len(results[i]))
        for formulation in MultisetFormulation.difference(InfinityFormulation):
            with self.subTest("check for no combinatorial duplicates", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("2[a* b*] \n a b")
                results = {
                    i: list(solver.configs_with_number_of_polymers(test_tbn, number_of_polymers=i, formulation=formulation))
                        for i in range(1, 5)
                }
                expected_number_of_configurations = {1: 1, 2: 1, 3: 0, 4: 0}
                for i in expected_number_of_configurations.keys():
                    self.assertEqual(expected_number_of_configurations[i], len(results[i]))

            with self.subTest("check for return of multiple stable configurations", parameter_string=parameter_string, formulation=formulation):
                test_tbn = Tbn.from_string("6(a*) \n 2[3(a*)] \n a \n 5(a) \n 2(a) \n 4(a)")
                results = {
                    i: list(solver.configs_with_number_of_polymers(test_tbn, number_of_polymers=i, formulation=formulation))
                        for i in range(1, 5)
                }
                expected_number_of_configurations = {1: 1, 2: 3, 3: 0, 4: 0}
                for i in expected_number_of_configurations.keys():
                    self.assertEqual(expected_number_of_configurations[i], len(results[i]))
