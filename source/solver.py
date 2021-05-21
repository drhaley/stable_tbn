from typing import Iterator, Optional
from enum import Enum, auto

from source.tbn import Tbn
from source.configuration import Configuration
from source.solver_adapters import constraint_programming, integer_programming
from source.constraints import Constraints

from source.formulations.bond_aware_network import Formulation as BondAwareNetworkFormulation
from source.formulations.bond_oblivious_network import Formulation as BondObliviousNetworkFormulation
from source.formulations.polymer_binary_matrix import Formulation as PolymerBinaryMatrixFormulation
from source.formulations.polymer_integer_matrix import Formulation as PolymerIntegerMatrixFormulation
from source.formulations.polymer_unbounded_matrix import Formulation as PolymerUnboundedMatrixFormulation
from source.formulations.variable_bond_weight import Formulation as VariableBondWeightFormulation
from source.formulations.hilbert_basis import Formulation as HilbertBasisFormulation


class SolverMethod(Enum):
    CONSTRAINT_PROGRAMMING = auto()
    INTEGER_PROGRAMMING = auto()  # Only implemented for single queries


class SolverFormulation(Enum):
    BOND_AWARE_NETWORK = auto()
    BOND_OBLIVIOUS_NETWORK = auto()
    POLYMER_BINARY_MATRIX = auto()
    POLYMER_INTEGER_MATRIX = auto()
    POLYMER_UNBOUNDED_MATRIX = auto()
    VARIABLE_BOND_WEIGHT = auto()
    HILBERT_BASIS = auto()


class Solver:
    def __init__(
            self,
            method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
    ):
        if method == SolverMethod.CONSTRAINT_PROGRAMMING:
            self.__single_solve_adapter = constraint_programming.Solver()
            self.__sorted_polymers = True
        elif method == SolverMethod.INTEGER_PROGRAMMING:
            self.__single_solve_adapter = integer_programming.Solver()
            self.__sorted_polymers = False
        else:
            raise NotImplementedError(f"solver not implemented for method {method}")

        self.__multi_solve_adapter = constraint_programming.Solver()  # only implemented for CP to solve_all

    def stable_config(self,
                      tbn: Tbn,
                      user_constraints: Constraints = Constraints(),
                      formulation: SolverFormulation = SolverFormulation.POLYMER_UNBOUNDED_MATRIX,
                      bond_weighting_factor: Optional[float] = None,
                      verbose: bool = False,
                      ) -> Configuration:
        if bond_weighting_factor is not None:
            user_constraints = user_constraints.with_bond_weight(bond_weighting_factor)

        if formulation == SolverFormulation.BOND_AWARE_NETWORK:
            formulation = BondAwareNetworkFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.BOND_OBLIVIOUS_NETWORK:
            formulation = BondObliviousNetworkFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.POLYMER_BINARY_MATRIX:
            formulation = PolymerBinaryMatrixFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.POLYMER_INTEGER_MATRIX:
            formulation = PolymerIntegerMatrixFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.POLYMER_UNBOUNDED_MATRIX:
            formulation = PolymerUnboundedMatrixFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.VARIABLE_BOND_WEIGHT:
            formulation = VariableBondWeightFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        elif formulation == SolverFormulation.HILBERT_BASIS:
            formulation = HilbertBasisFormulation(tbn, self.__single_solve_adapter, user_constraints)
            return formulation.get_configuration(verbose=verbose)
        else:
            raise AssertionError(f"did not recognize formulation requested: {formulation}")

    def stable_configs(self,
                       tbn: Tbn,
                       user_constraints: Constraints = Constraints(),
                       formulation: SolverFormulation = SolverFormulation.POLYMER_UNBOUNDED_MATRIX,
                       bond_weighting_factor: Optional[float] = None,
                       verbose: bool = False,
                       ) -> Iterator[Configuration]:
        if bond_weighting_factor is not None:
            user_constraints = user_constraints.with_bond_weight(bond_weighting_factor)

        if user_constraints.optimize():  # do a first solve to find optimal objective value
            example_stable_configuration = self.stable_config(
                tbn,
                user_constraints=user_constraints,
                formulation=formulation,
                bond_weighting_factor=bond_weighting_factor,
                verbose=verbose,
            )
            user_constraints = user_constraints.with_unset_optimization_flag()
            number_of_polymers = example_stable_configuration.number_of_polymers()
            number_of_merges = example_stable_configuration.number_of_merges()
            amount_of_energy = example_stable_configuration.energy(user_constraints.bond_weight())
            fixed_polymer_user_constraints = user_constraints.with_fixed_polymers(number_of_polymers)
            fixed_merge_user_constraints = user_constraints.with_fixed_merges(number_of_merges)
            fixed_energy_user_constraints = user_constraints.with_fixed_energy(amount_of_energy)
        else:
            fixed_polymer_user_constraints = user_constraints
            fixed_merge_user_constraints = user_constraints
            fixed_energy_user_constraints = user_constraints

        if formulation == SolverFormulation.BOND_AWARE_NETWORK:
            formulation = BondAwareNetworkFormulation(
                tbn, self.__multi_solve_adapter, fixed_polymer_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.BOND_OBLIVIOUS_NETWORK:
            formulation = BondObliviousNetworkFormulation(
                tbn, self.__multi_solve_adapter, fixed_polymer_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.POLYMER_BINARY_MATRIX:
            formulation = PolymerBinaryMatrixFormulation(
                tbn, self.__multi_solve_adapter, fixed_polymer_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.POLYMER_INTEGER_MATRIX:
            formulation = PolymerIntegerMatrixFormulation(
                tbn, self.__multi_solve_adapter, fixed_polymer_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.POLYMER_UNBOUNDED_MATRIX:
            formulation = PolymerUnboundedMatrixFormulation(
                tbn, self.__multi_solve_adapter, fixed_merge_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.VARIABLE_BOND_WEIGHT:
            formulation = VariableBondWeightFormulation(
                tbn, self.__multi_solve_adapter, fixed_energy_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        elif formulation == SolverFormulation.HILBERT_BASIS:
            formulation = HilbertBasisFormulation(
                tbn, self.__multi_solve_adapter, fixed_polymer_user_constraints
            )
            return formulation.get_all_configurations(verbose=verbose)

        else:
            raise AssertionError(f"did not recognize formulation requested: {formulation}")
