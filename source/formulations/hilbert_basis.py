import os
import sys
from random import randint
import subprocess
import numpy as np
from typing import List, Dict, Any
from math import inf as infinity

from source.formulations.abstract import Formulation as AbstractFormulation
from source.tbn import Tbn
from source.configuration import Configuration
from source.polymer import Polymer

# TODO: Hilbert basis implementation does not respond or assert on some user constraints (merges, energy)
# TODO: implement verbosity flags


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        # if constrain_quantities is enabled, truncates the Hilbert basis based upon available monomer counts
        #  sadly it does not seem to improve the runtime to enable this
        constrain_quantities = False

        self._run_asserts()
        monomers_and_slack_as_matrix = self._project_tbn_to_column_matrix(self.tbn)
        polymer_basis = self._get_hilbert_basis_from_matrix(monomers_and_slack_as_matrix,
                tbn = self.tbn if constrain_quantities else None)
        self._display_polymer_basis(polymer_basis)
        self._populate_model_from_tbn_and_polymer_basis(self.tbn, polymer_basis)

    @staticmethod
    def _project_tbn_to_column_matrix(tbn: Tbn) -> np.array:
        monomer_types = list(tbn.monomer_types())
        limiting_domain_types = list(tbn.limiting_domain_types())
        monomer_matrix = np.array([
            [-monomer.net_count(domain) for domain in limiting_domain_types]
                for monomer in monomer_types
        ], np.int64).T

        return monomer_matrix

    @staticmethod
    def _get_hilbert_basis_from_matrix(matrix: np.array, tbn: Tbn = None, quiet: bool = False) -> np.array:
        # if a TBN is specified, the Hilbert basis will be constrained to the quantities of the monomers in the TBN
        temporary_filename_prefix = os.path.join(os.getcwd(), f"TEMP_4ti2_CALLOUT_{randint(0,10000)}")
        matrix_filename = temporary_filename_prefix + '.mat'
        relations_filename = temporary_filename_prefix + '.rel'
        signs_filename = temporary_filename_prefix + '.sign'
        rhs_filename = temporary_filename_prefix + '.rhs'
        homogenous_basis_filename = temporary_filename_prefix + '.zhom'
        inhomogenous_basis_filename = temporary_filename_prefix + '.zinhom'

        number_of_domain_types = matrix.shape[0]
        number_of_monomer_types = matrix.shape[1]

        if tbn:
            identity_matrix = np.identity(number_of_monomer_types, np.int64)
            matrix = np.concatenate((matrix, identity_matrix))

        # write a temporary matrix file
        with open(matrix_filename, 'w') as outFile:
            # shape is first line of the .mat file
            number_of_rows = number_of_domain_types + (number_of_monomer_types if tbn else 0)
            outFile.write(f"{number_of_rows} {number_of_monomer_types}\n")
            for row in matrix:
                for entry in row:
                    outFile.write(f"{entry} ")
                outFile.write('\n')

        # write a temporary relations file
        with open(relations_filename, 'w') as outFile:
            entries = number_of_domain_types + (number_of_monomer_types if tbn else 0)
            outFile.write(f"1 {entries}\n")
            outFile.write(' '.join(['>']*number_of_domain_types))
            if tbn:
                outFile.write(' ')
                outFile.write(' '.join(['<']*number_of_monomer_types))
            outFile.write('\n')

        # write a temporary signs file
        with open(signs_filename, 'w') as outFile:
            outFile.write(f"1 {number_of_monomer_types}\n")
            outFile.write(' '.join(['1']*number_of_monomer_types))
            outFile.write('\n')

        # write a temporary right hand sides file
        with open(rhs_filename, 'w') as outFile:
            entries = number_of_domain_types + (number_of_monomer_types if tbn else 0)
            outFile.write(f"1 {entries}\n")
            outFile.write(' '.join(['0']*number_of_domain_types))
            if tbn:
                outFile.write(' ')
                outFile.write(' '.join([str(tbn.count(monomer)) for monomer in tbn.monomer_types()]))
            outFile.write('\n')

        try:
            subprocess.call(["4ti2-zsolve", temporary_filename_prefix] + (['-q'] if quiet else []))

            # read and interpret response
            with open(homogenous_basis_filename, 'r') as inFile:
                # shape is first line of the file
                shape_as_string = inFile.readline()
                shape = tuple(int(x) for x in shape_as_string.split())

                flat_array = np.array(inFile.read().split(), np.int64)
                homogenous_basis = flat_array.reshape(*shape)

            with open(inhomogenous_basis_filename, 'r') as inFile:
                # shape is first line of the file
                shape_as_string = inFile.readline()
                shape = tuple(int(x) for x in shape_as_string.split())

                flat_array = np.array(inFile.read().split(), np.int64)
                inhomogenous_basis = flat_array.reshape(*shape)

                # remove the row of all zeros
                inhomogenous_basis = inhomogenous_basis[~np.all(inhomogenous_basis == 0, axis=1)]

            if homogenous_basis.size == 0:
                basis = inhomogenous_basis.T
            elif inhomogenous_basis.size == 0:
                basis = homogenous_basis.T
            else:
                basis = np.concatenate(homogenous_basis, inhomogenous_basis).T

            print('='*80)
            print(homogenous_basis)
            print('-' * 80)
            print(inhomogenous_basis)
            print('-' * 80)
            print(basis)
            print('=' * 80)
        except FileNotFoundError:
            print("4ti2 was not able to complete")
            basis = None
        finally:
            # clean up
            for filename in [
                    matrix_filename,
                    relations_filename,
                    signs_filename,
                    rhs_filename,
                    inhomogenous_basis_filename,
                    homogenous_basis_filename,
                    ]:
                try:
                    os.remove(filename)
                except FileNotFoundError:
                    pass
        if basis.size > 0:
            return basis
        else:
            raise AssertionError("the callout to 4ti2 did not generate a Hilbert basis")

    def _populate_model_from_tbn_and_polymer_basis(self, tbn: Tbn, basis: np.array) -> None:
        self.polymer_basis = basis

        # upper bound on how many total monomers can be in non-singleton polymers
        upper_bound_on_total_monomers_in_complexes = sum(
            tbn.count(monomer_type) * (1 + abs(monomer_type.net_count(domain_type)))
            for monomer_type in tbn.limiting_monomer_types()
            for domain_type in tbn.limiting_domain_types()
        )
        monomer_counts = [tbn.count(monomer) for monomer in tbn.monomer_types()]
        upper_bound_on_total_monomers_in_complexes = min(
            upper_bound_on_total_monomers_in_complexes,
            sum(monomer_counts)   #total number of monomers
        )

        BIG_M = upper_bound_on_total_monomers_in_complexes

        self.model.set_big_m(BIG_M)

        # declare variables
        size_of_basis_vectors, number_of_basis_vectors = basis.shape
        self.basis_coefficients = [
            self.model.int_var(0, BIG_M, f"basis_coefficient_{i}")
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
                raise NotImplementedError("Not implemented to use Hilbert basis formulation with infinite monomer counts")
