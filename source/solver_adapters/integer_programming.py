from typing import Any, Iterator, Dict, List
from ortools.linear_solver import pywraplp
from source.solver_adapters import abstract


class IpModel(abstract.Model, pywraplp.Solver):
    def __init__(self):
        # calling superclasses explicitly here because of multiple inheritance
        abstract.Model.__init__(self)
        pywraplp.Solver.__init__(self, "stable_tbn-ip-model", pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
        self.__id_counter = 0

    def __get_id(self) -> int:
        self.__id_counter += 1
        return self.__id_counter

    def int_var(self, *args, **kargs) -> pywraplp.Variable:
        return self.IntVar(*args, **kargs)

    def bool_var(self, *args, **kargs) -> pywraplp.Variable:
        return self.BoolVar(*args, **kargs)

    def complement_var(self, var: pywraplp.Variable) -> Any:
        return -var + 1

    def __get_big_M(self) -> int:
        if self._big_M is None:
            raise AssertionError("Cannot make an implication without a big M value for integer programming")
        else:
            return int(self._big_M)

    def AddChainedImplication(self, *args) -> Any:
        # intended call: .AddChainedImplication(antecedent1, antecedent2, ..., consequent)
        big_M = self.__get_big_M()
        if len(args) < 2:
            raise AssertionError(
                "Call to AddChainedImplication with less than two arguments.  Need antecedent and consequent"
            )
        else:
            consequent = args[-1]
            antecedents = args[:-1]

        # delta1 and delta2 will be the indicators for the nonzeroness (i.e. truthfulness) of the consequent.
        # If consequent < 0, delta1 can be 1
        # If consequent > 0, delta2 can be 1
        # we force either of delta1 or delta2 to be the value 1
        delta1 = self.bool_var(f'indicator_chain1_{self.__get_id()}')
        delta2 = self.bool_var(f'indicator_chain2_{self.__get_id()}')
        self.Add(consequent <= -1 + big_M * (-delta1 + 1))
        self.Add(consequent >= 1 - big_M * (-delta2 + 1))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        constraint = self.Add(
            delta1 + delta2 + sum(self.complement_var(antecedent) for antecedent in antecedents) >= 1
        )
        return constraint

    def AddEqualToZeroImplication(self, *args) -> Any:
        # intended call: .AddZeroValueImplication(antecedent1, antecedent2, ..., consequent, big_M = very_large_value)
        big_M = self.__get_big_M()
        if len(args) < 2:
            raise AssertionError(
                "Call to AddEqualToZeroImplication with less than two arguments.  Need antecedent and consequent"
            )
        else:
            consequent = args[-1]
            antecedents = args[:-1]

        # delta will be the indicator for the zeroness of the consequent.
        # If consequent == 0, delta == True
        delta = self.bool_var(f'indicator_zero_{self.__get_id()}')
        self.Add(consequent <= big_M * (-delta + 1))
        self.Add(consequent >= -big_M * (-delta + 1))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        constraint = self.Add(
            delta + sum(self.complement_var(antecedent) for antecedent in antecedents) >= 1
        )
        return constraint

    def AddGreaterThanZeroImplication(self, *args) -> Any:
        # intended call: .AddGreaterThanImplication(antecedent1, antecedent2, ..., consequent, big_M = very_large_value)
        big_M = self.__get_big_M()
        if len(args) < 2:
            raise AssertionError(
                "Call to AddGreaterThanImplication with less than two arguments.  Need antecedent and consequent"
            )
        else:
            consequent = args[-1]
            antecedents = args[:-1]

        # delta will be the indicator for consequent greater than zero (e.g. if consequent > 0, delta == True)
        delta = self.bool_var(f'indicator_gt_zero_{self.__get_id()}')
        self.Add(consequent >= 1 - (big_M * (-delta + 1)))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        constraint = self.Add(
            delta + sum(self.complement_var(antecedent) for antecedent in antecedents) >= 1
        )
        return constraint


class Solver(abstract.SolverAdapter):
    @staticmethod
    def model() -> abstract.Model:
        return IpModel()

    def solve(self, model: abstract.Model, variables_with_values_to_keep: List[Any]) -> Any:
        status = model.Solve()
        return status

    def value(self, var: pywraplp.Variable) -> int:
        return int(var.solution_value())

    def solve_all(self, model: abstract.Model, variables_with_values_to_keep: List[Any]) -> Iterator[Dict[Any, int]]:
        raise NotImplementedError("Not implemented to query the complete solution set using IP")
