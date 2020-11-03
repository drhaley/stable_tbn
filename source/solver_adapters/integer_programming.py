from typing import Any, Iterator, Dict, List, Union
from ortools.linear_solver import pywraplp
from source.solver_adapters import abstract


class IpModel(abstract.Model, pywraplp.Solver):
    def __init__(self):
        # calling superclasses explicitly here because of multiple inheritance
        abstract.Model.__init__(self)
        pywraplp.Solver.__init__(self, "stable_tbn-ip-model", pywraplp.Solver.SCIP_MIXED_INTEGER_PROGRAMMING)
        self.__id_counter = 0
        self.OPTIMAL = pywraplp.Solver.OPTIMAL
        self.INFEASIBLE = pywraplp.Solver.INFEASIBLE

    def __get_id(self) -> int:
        self.__id_counter += 1
        return self.__id_counter

    def int_var(self, *args, **kargs) -> pywraplp.Variable:
        return self.IntVar(*args, **kargs)

    def bool_var(self, *args, **kargs) -> pywraplp.Variable:
        return self.BoolVar(*args, **kargs)

    def complement_var(self, var: pywraplp.Variable) -> Any:
        return -var + 1

    def __get_big_m(self) -> int:
        if self._big_M is None:
            raise AssertionError("Cannot make an implication without a big M value for integer programming")
        else:
            return int(self._big_M)

    def add_constraint(self, *args) -> Any:
        return self.Add(*args)

    def add_implication(self, *args) -> Any:
        # intended call: .AddChainedImplication(antecedent1, antecedent2, ..., consequent)
        big_m = self.__get_big_m()
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
        self.add_constraint(consequent <= -1 + big_m * (-delta1 + 1))
        self.add_constraint(consequent >= 1 - big_m * (-delta2 + 1))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        constraint = self.add_constraint(
            delta1 + delta2 + sum(self.complement_var(antecedent) for antecedent in antecedents) >= 1
        )
        return constraint

    def add_equal_to_zero_implication(self, *args) -> Any:
        # intended call:
        #   .add_equal_to_zero_implication(antecedent1, antecedent2, ..., consequent, big_m = very_large_value)
        big_m = self.__get_big_m()
        if len(args) < 2:
            raise AssertionError(
                "Call to add_equal_to_zero_implication with less than two arguments.  Need antecedent and consequent"
            )
        else:
            consequent = args[-1]
            antecedents = args[:-1]

        # delta will be the indicator for the zeroness of the consequent.
        # If consequent == 0, delta == True
        delta = self.bool_var(f'indicator_zero_{self.__get_id()}')
        self.add_constraint(consequent <= big_m * (-delta + 1))
        self.add_constraint(consequent >= -big_m * (-delta + 1))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        boolean_complement_sum: pywraplp.Variable = sum(self.complement_var(antecedent) for antecedent in antecedents)
        constraint = self.add_constraint(delta + boolean_complement_sum >= 1)
        return constraint

    def add_greater_than_zero_implication(self, *args) -> Any:
        # intended call:
        #   .add_greater_than_zero_implication(antecedent1, antecedent2, ..., consequent, big_m = very_large_value)
        big_m = self.__get_big_m()
        if len(args) < 2:
            raise AssertionError(
                "Call to add_greater_than_zero_implication with less than two arguments."
            )
        else:
            consequent = args[-1]
            antecedents = args[:-1]

        # delta will be the indicator for consequent greater than zero (e.g. if consequent > 0, delta == True)
        delta = self.bool_var(f'indicator_gt_zero_{self.__get_id()}')
        self.add_constraint(consequent >= 1 - (big_m * (-delta + 1)))

        # for p => (q => delta), change to delta + (1-p) + (1-q) >= 1.  for more antecedents, extrapolate
        boolean_complement_sum: pywraplp.Variable = sum(self.complement_var(antecedent) for antecedent in antecedents)
        constraint = self.add_constraint(delta + boolean_complement_sum >= 1)
        return constraint

    def minimize(self, *args, **kargs) -> None:
        self.Minimize(*args, **kargs)

    def maximize(self, *args, **kargs) -> None:
        self.Maximize(*args, **kargs)


class Solver(abstract.SolverAdapter):
    @staticmethod
    def model() -> abstract.Model:
        return IpModel()

    def solve(self, model: Union[abstract.Model, IpModel],
              variables_with_values_to_keep: List[Any], verbose: bool = False) -> Any:
        if verbose:
            model.EnableOutput()
        status = model.Solve()
        return status

    def value(self, var: Union[int, pywraplp.Variable]) -> int:
        if isinstance(var, int):
            return var
        else:
            return int(var.solution_value())

    def solve_all(self, model: abstract.Model, variables_with_values_to_keep: List[Any], verbose: bool = False)\
            -> Iterator[Dict[Any, int]]:
        raise NotImplementedError("Not implemented to query the complete solution set using IP")
