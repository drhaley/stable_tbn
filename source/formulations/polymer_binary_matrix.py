from typing import Iterator, List, Any, Dict

from source.formulations.abstract import Formulation as AbstractFormulation
from source.tbn import Tbn
from source.configuration import Configuration
from source.constraints import Constraints


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        self._add_variables()
        self._add_constraints()

    def _add_variables(self) -> None:
        pass
        # self.model.bool_var()
        # self.model.int_var()

    def _add_constraints(self) -> None:
        pass
        # self.model.add_constraint()

    def _variables_to_keep(self) -> List[Any]:
        """
        returns a list of the variables that are necessary to convert a solution back to a configuration
          (e.g. returns polymer composition variables but not 'internal' tie-breaker variables)
        """
        pass

    def _interpret_solution(self, variable_to_value_dictionary: Dict[Any, int]) -> Configuration:
        """
        uses the provided dictionary to convert solution variables into solution values and from this,
          converts the solutions values into the corresponding configuration
        """
        pass
