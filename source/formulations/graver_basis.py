import os, sys
from random import randint
import subprocess
import numpy as np
from typing import List, Tuple, Dict, Any, Set
from math import inf as infinity

from source.formulations.abstract import Formulation as AbstractFormulation
from source.tbn import Tbn
from source.solver_adapters.abstract import Model
from source.configuration import Configuration
from source.polymer import Polymer

# TODO: Graver basis implementation does not respond or assert on some user constraints (merges, energy)
# TODO: implement verbosity flags


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        self._run_asserts()
        monomers_and_slack_as_matrix = self._project_tbn_to_column_matrix(self.tbn)
        polymer_basis = self._get_hilbert_basis_from_matrix(monomers_and_slack_as_matrix)
        self._display_polymer_basis(polymer_basis)
        self._populate_model_from_polymer_basis(polymer_basis)

    @staticmethod
    def _project_tbn_to_column_matrix(tbn: Tbn) -> np.array:
        monomer_types = list(tbn.monomer_types())
        limiting_domain_types = list(tbn.limiting_domain_types())
        filtered_limiting_domain_types = list(tbn.limiting_domain_types(filter_ties=True))
        monomer_matrix = np.array([
            [monomer.net_count(domain) for domain in limiting_domain_types]
                for monomer in monomer_types
        ], np.int64).T
        if filtered_limiting_domain_types:  # is not empty
            slack_matrix = np.array([
                [1 if filtered_domain == domain else 0 for domain in limiting_domain_types]
                    for filtered_domain in filtered_limiting_domain_types
            ], np.int64)
            full_matrix = np.hstack((monomer_matrix, slack_matrix))
        else:
            full_matrix = monomer_matrix
        return full_matrix

    @staticmethod
    def _get_hilbert_basis_from_matrix(matrix: np.array, quiet: bool = False) -> np.array:
        temporary_filename_prefix = os.path.join(os.getcwd(), f"TEMP_4ti2_CALLOUT_{randint(0,10000)}")
        matrix_filename = temporary_filename_prefix + '.mat'
        hilbert_basis_filename = temporary_filename_prefix + '.hil'

        # write a temporary matrix file
        with open(matrix_filename, 'w') as outFile:
            # shape is first line of the .mat file
            outFile.write(f"{matrix.shape[0]} {matrix.shape[1]}\n")
            for row in matrix:
                for entry in row:
                    outFile.write(f"{entry} ")
                outFile.write('\n')

        subprocess.call(["hilbert", temporary_filename_prefix, "-q" if quiet else ""])

        # read and interpret response
        with open(hilbert_basis_filename, 'r') as inFile:
            # shape is first line of the file
            shape_as_string = inFile.readline()
            shape = tuple(int(x) for x in shape_as_string.split())

            flat_array = np.array(inFile.read().split(), np.int64)
            basis = flat_array.reshape(*shape).T

        # clean up
        os.remove(matrix_filename)
        os.remove(hilbert_basis_filename)

        return basis

    def _populate_model_from_polymer_basis(self, basis: np.array) -> Model:
        self.polymer_basis = basis

        self.model.set_big_m(500)  # TODO: replace 500 with appropriate bound

        # declare variables
        size_of_basis_vectors, number_of_basis_vectors = basis.shape
        self.basis_coefficients = [
            self.model.int_var(0, 500, f"basis_coefficient_{i}")  # TODO: replace 500 with an appropriate upper bound
                for i in range(number_of_basis_vectors)
        ]

        # conservation constraints
        for i, monomer_type in enumerate(list(self.tbn.monomer_types())):
            self.model.Add(
                self.tbn.count(monomer_type) ==
                sum(
                    self.basis_coefficients[j] * basis_vector[i]
                        for j, basis_vector in enumerate(basis.T)
                )
            )

        # counting constraints
        number_of_polymers = sum(self.basis_coefficients)

        if self.user_constraints.max_polymers() != infinity:
            self.model.add_constraint(number_of_polymers <= self.user_constraints.max_polymers())
        if self.user_constraints.min_polymers() > 0:
            self.model.add_constraint(number_of_polymers >= self.user_constraints.min_polymers())

        # objective function
        if self.user_constraints.optimize():
            self.model.maximize(number_of_polymers)

    def _display_polymer_basis(self, basis: np.array) -> None:
        print("Found Polymer basis:")
        print("-------------------")
        for basis_vector in basis.T:
            polymer = self._make_polymer_from_vector(basis_vector)
            print(polymer)
        print("-------------------")

    def _make_polymer_from_vector(self, basis_vector: np.array) -> Polymer:
        this_polymer_dict = {}
        for i, monomer in enumerate(self.tbn.monomer_types()):
            monomer_count = int(basis_vector[i])
            if monomer_count > 0:
                this_polymer_dict[monomer] = monomer_count
        this_polymer = Polymer(this_polymer_dict)
        return this_polymer

    def _variables_to_keep(self) -> List[Any]:
        """
        returns a list of the variables that are necessary to convert a solution back to a configuration
          (e.g. returns polymer composition variables but not 'internal' tie-breaker variables)
        """
        return self.basis_coefficients

    def _interpret_solution(self, variable_to_value_dictionary: Dict[Any, int]) -> Configuration:
        """
        uses the provided dictionary to convert solution variables into solution values and from this,
          converts the solutions values into the corresponding configuration
        """
        this_configuration_dict = {}
        for i, basis_vector in enumerate(self.polymer_basis.T):
            this_polymer = self._make_polymer_from_vector(basis_vector)
            count = variable_to_value_dictionary[self.basis_coefficients[i]]
            if count > 0:
                this_configuration_dict[this_polymer] = count
        return Configuration(this_configuration_dict)

    def _run_asserts(self) -> None:
        for monomer_type in self.tbn.monomer_types():
            if self.tbn.count(monomer_type) == infinity:
                raise NotImplementedError("Not implemented to use graver basis formulation with infinite monomer counts")
