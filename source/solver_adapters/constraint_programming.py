from typing import Any, List, Iterator, Dict, Union
from ortools.sat.python import cp_model
from source.solver_adapters import abstract


class CpModel(abstract.Model, cp_model.CpModel):
    def __init__(self):
        # calling superclasses explicitly here because of multiple inheritance
        cp_model.CpModel.__init__(self)
        abstract.Model.__init__(self)
        self.OPTIMAL = cp_model.OPTIMAL
        self.INFEASIBLE = cp_model.INFEASIBLE

    def int_var(self, *args, **kargs) -> cp_model.IntVar:
        return self.NewIntVar(*args, **kargs)

    def bool_var(self, *args, **kargs) -> cp_model.IntVar:  # ortools internally stores this as type IntVar, not BoolVar
        return self.NewBoolVar(*args, **kargs)

    def complement_var(self, var: cp_model.IntVar) -> cp_model.IntVar:
        return var.Not()

    def add_constraint(self, *args) -> Any:
        return self.Add(*args)

    def add_implication(self, *args) -> Any:
        # intended call: .AddChainedImplication(antecedent1, antecedent2, ..., consequent)
        if len(args) < 2:
            raise AssertionError(
                "Called AddChainedImplication with too few arguments.  Need antecedent(s) & consequent"
            )
        else:
            consequent = args[-1]
            antecedents = list(args[:-1])

        for antecedent in antecedents:
            if type(antecedent) is int and antecedent == 0:  # have to do it this way because bool vars in CP == 0
                return None
        else:
            if type(consequent) is not cp_model.BoundedLinearExpression:
                constraint = self.add_constraint(consequent != int(False))
            else:
                constraint = self.add_constraint(consequent)

            constraint = constraint.OnlyEnforceIf([x for x in antecedents if x != 1 and x is not True])

            return constraint

    def add_equal_to_zero_implication(self, *args) -> Any:
        return self.add_implication(*args[:-1], args[-1] == 0)

    def add_greater_than_zero_implication(self, *args) -> Any:
        return self.add_implication(*args[:-1], args[-1] > 0)

    def minimize(self, *args, **kargs) -> None:
        self.Minimize(*args, **kargs)

    def maximize(self, *args, **kargs) -> None:
        self.Maximize(*args, **kargs)


class Solver(abstract.SolverAdapter):
    def __init__(self):
        super().__init__()
        self.__internal_solver = None

    @staticmethod
    def model() -> abstract.Model:
        return CpModel()

    def solve(self, model: abstract.Model, variables_with_values_to_keep: List[Any], verbose: bool = False) -> Any:
        self.__internal_solver = cp_model.CpSolver()
        self.__internal_solver.parameters.log_search_progress = verbose
        status = self.__internal_solver.Solve(model)
        return status

    def value(self, var: Union[int, cp_model.IntVar]) -> int:
        if isinstance(var, int):
            return var
        else:
            return self.__internal_solver.Value(var)

    def solve_all(self, model: abstract.Model, variables_with_values_to_keep, verbose: bool = False)\
            -> Iterator[Dict[Any, int]]:
        internal_solver = cp_model.CpSolver()
        internal_solver.parameters.log_search_progress = verbose

        found_solutions = []
        solution_accumulator = SolutionAccumulator(variables_with_values_to_keep, found_solutions)
        status = internal_solver.SearchForAllSolutions(model, solution_accumulator)

        if status == cp_model.INFEASIBLE:
            return iter(())
        elif status != cp_model.OPTIMAL:
            raise AssertionError(f"OR-Tools returned code {status}, but expected {cp_model.OPTIMAL}")

        return iter(found_solutions)


class SolutionAccumulator(cp_model.CpSolverSolutionCallback):
    def __init__(
            self,
            variables_with_values_to_keep: List[Any],
            found_solutions: List[Dict[Any, int]],  # solutions are passed out by reference
    ):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables_with_values_to_keep = variables_with_values_to_keep
        self.__found_solutions = found_solutions

    def on_solution_callback(self) -> None:
        this_solution = {v: self.Value(v) for v in self.__variables_with_values_to_keep}
        self.__found_solutions.append(this_solution)
