import numpy as np
from typing import List, Tuple, Dict, Any, Set
from math import inf as infinity

from source.formulations.abstract import Formulation as AbstractFormulation
from source.tbn import Tbn
from source.solver_adapters.abstract import Model
from source.configuration import Configuration
from source.polymer import Polymer


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        monomers_and_slack_as_matrix = self._project_tbn_to_column_matrix(self.tbn)
        H, C = self._convert_matrix_to_hermite_normal_form(monomers_and_slack_as_matrix)
        kernel_basis = self._get_kernel_basis_from_hermite_normal_form(H, C)
        graver_basis = self._get_graver_basis_from_kernel_basis(kernel_basis)
        self.model = self._get_model_from_graver_basis(graver_basis)

    @staticmethod
    def _project_tbn_to_column_matrix(tbn: Tbn) -> np.array:
        pass

    @staticmethod
    def _convert_matrix_to_hermite_normal_form(monomers_and_slack_as_matrix: np.array) -> Tuple[np.array, np.array]:
        pass

    @staticmethod
    def _get_kernel_basis_from_hermite_normal_form(H: np.array, C: np.array) -> np.array:
        pass

    @staticmethod
    def _get_graver_basis_from_kernel_basis(kernel_basis: np.array) -> np.array:
        pass

    @staticmethod
    def _get_model_from_graver_basis(graver_basis: np.array) -> Model:
        pass

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
