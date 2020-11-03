from abc import ABC, abstractmethod
from typing import Iterator, List, Any, Dict

from source.tbn import Tbn
from source.configuration import Configuration
from source.constraints import Constraints
from source.solver_adapters.abstract import SolverAdapter


class Formulation(ABC):
    def __init__(self, tbn: Tbn, solver: SolverAdapter, user_constraints: Constraints = Constraints()) -> None:
        self.tbn = tbn
        self.solver = solver
        self.user_constraints = user_constraints
        self.model = self.solver.model()
        self._populate_model()

    def get_configuration(self, verbose: bool = False) -> Configuration:
        solution_status = self.solver.solve(self.model, self._variables_to_keep(), verbose=verbose)
        self._assert_completed_status(solution_status)
        variable_to_value_dictionary = {
            var: self.solver.value(var) for var in self._variables_to_keep()
        }
        return self._interpret_solution(variable_to_value_dictionary)

    def get_all_configurations(self, verbose: bool = False) -> Iterator[Configuration]:
        for variable_to_value_dictionary in self.solver.solve_all(
                self.model, self._variables_to_keep(), verbose=verbose
        ):
            yield self._interpret_solution(variable_to_value_dictionary)

    def _assert_completed_status(self, status: int):
        if status == self.model.INFEASIBLE:
            raise AssertionError(f"Could not find solution to tbn, was reported infeasible")
        elif status != self.model.OPTIMAL:
            raise AssertionError(f"could not find solution to tbn, got code {status} instead")
        else:
            pass

    @abstractmethod
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        pass

    @abstractmethod
    def _variables_to_keep(self) -> List[Any]:
        """
        returns a list of the variables that are necessary to convert a solution back to a configuration
          (e.g. returns polymer composition variables but not 'internal' tie-breaker variables)
        """
        pass

    @abstractmethod
    def _interpret_solution(self, variable_to_value_dictionary: Dict[Any, int]) -> Configuration:
        """
        uses the provided dictionary to convert solution variables into solution values and from this,
          converts the solutions values into the corresponding configuration
        """
        pass
