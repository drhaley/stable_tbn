from typing import Iterator, Optional, Tuple, List, Dict, Any
from enum import Enum, auto
from math import inf as infinity

from source.tbn import Tbn
from source.configuration import Configuration
from source.polymer import Polymer
from source.monomer import Monomer
from source.solver_adapters import constraint_programming, integer_programming
from source.solver_adapters import abstract as adapter


class SolverMethod(Enum):
    CONSTRAINT_PROGRAMMING = auto()
    INTEGER_PROGRAMMING = auto()  # Only implemented for single queries


class Solver:
    def __init__(
                self,
                method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
    ):
        if method == SolverMethod.CONSTRAINT_PROGRAMMING:
            self.__adapter = constraint_programming.Solver()
        elif method == SolverMethod.INTEGER_PROGRAMMING:
            self.__adapter = integer_programming.Solver()
        else:
            raise NotImplementedError(f"solver not implemented for method {method}")

    def stable_config(self, tbn: Tbn) -> Configuration:
        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(tbn)
        status = self.__adapter.solve(model, polymer_composition_vars.values())
        if status is not model.OPTIMAL:
            raise AssertionError(f"could not find optimal solution to tbn, got {status} instead")
        else:
            value_dict = {var: self.__adapter.value(var) for var in polymer_composition_vars.values()}
            return self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, value_dict)

    def stable_configs(self, tbn: Tbn) -> Iterator[Configuration]:
        example_stable_configuration = self.stable_config(tbn)
        number_of_merges = example_stable_configuration.number_of_merges()
        return self.configs_with_number_of_merges(tbn, number_of_merges)

    def configs_with_number_of_polymers(self, tbn: Tbn, number_of_polymers: Optional[int] = None) -> Iterator[Configuration]:
        total_number_of_monomers = sum(
            tbn.count(monomer_type)
            for monomer_type in tbn.monomer_types()
        )
        if total_number_of_monomers == infinity:
            raise AssertionError(
                "cannot request a specific number of polymers when there are an infinite amount of monomers"
            )
        return self.configs_with_number_of_merges(tbn, number_of_merges=total_number_of_monomers - number_of_polymers)

    def configs_with_number_of_merges(self,
                                      tbn: Tbn,
                                      number_of_merges:Optional[int] = None,
                                      non_singleton_polymer_estimate: Optional[int] = None,
                                      ) -> Iterator[Configuration]:
        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(
            tbn, fixed_number_of_merges=number_of_merges, non_singleton_polymer_estimate=non_singleton_polymer_estimate
        )
        solutions = self.__adapter.solve_all(model, polymer_composition_vars.values())
        return [
            self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, solution)
                for solution in solutions
        ]

    def __build_model(self,
                      tbn: Tbn,
                      fixed_number_of_merges: Optional[int] = None,
                      non_singleton_polymer_estimate: Optional[int] = None,
                      ) -> Tuple[adapter.Model, List[Monomer], Dict[Tuple[int, int], Any]]:
        """
        :param tbn:
        :param fixed_number_of_polymers:
        if a number of polymers is fixed, can supply all the configurations with this number of polymers,
         otherwise supplies a single configuration with the largest number of polymers possible (i.e. fewest merges)
        :return:
            (model, ordered_monomer_list, polymer_inclusion_matrix_entry_dict)
        """
        model = self.__adapter.model()

        limiting_monomer_types = list(tbn.limiting_monomer_types())
        ordered_monomer_types = list(tbn.monomer_types())
        limiting_domain_types = list(tbn.limiting_domain_types())

        total_number_of_limiting_monomers = sum(
            tbn.count(monomer_type)
            for monomer_type in limiting_monomer_types
        )

        upper_bound_on_monomers_in_complexes = sum(
            tbn.count(monomer_type) * (1 + abs(monomer_type.net_count(domain_type)))
            for monomer_type in limiting_monomer_types
            for domain_type in monomer_type.unstarred_domain_types()
        )

        model.set_big_M(upper_bound_on_monomers_in_complexes)  # needed for some formulations, like IP

        order_the_polymers = True  # setting to this to False is a minor improvement for IP, but awful for CP

        if non_singleton_polymer_estimate is None:
            number_of_possible_polymers = total_number_of_limiting_monomers
        else:
            number_of_possible_polymers = non_singleton_polymer_estimate

        ############ Variables ############
        # polymer_composition_vars[i, j] = count of monomer i in polymer j
        polymer_composition_vars = {}
        for i, monomer_type in enumerate(ordered_monomer_types):
            for j in range(number_of_possible_polymers):
                polymer_composition_vars[i, j] = model.int_var(0, upper_bound_on_monomers_in_complexes, f'polymer_composition_{i}_{j}')

        if order_the_polymers:
            # order-enforcing variables; enforces lexicographical ordering of polymers based upon monomer counts within
            tiebreaker_vars = {}
            for i in range(-1, len(ordered_monomer_types)):
                for j in range(number_of_possible_polymers - 1):
                    tiebreaker_vars[i, j] = model.bool_var(f'tiebreaker_{i}_{j}')

        # indicator variables for polymers
        indicator_vars = {}
        for j in range(number_of_possible_polymers):
            indicator_vars[j] = model.bool_var(f'indicator_{j}')

        ############ Constraints ############
        # monomer conservation; must use all limiting monomers, and cannot exceed the count of other monomers
        for i, monomer in enumerate(ordered_monomer_types):
            if monomer in limiting_monomer_types:
                model.Add(
                    sum(
                        polymer_composition_vars[i, j]
                        for j in range(number_of_possible_polymers)
                    ) == tbn.count(monomer)
                )
            elif tbn.count(monomer) < infinity:
                model.Add(
                    sum(
                        polymer_composition_vars[i, j]
                        for j in range(number_of_possible_polymers)
                    ) <= tbn.count(monomer)
                )

        # must saturate the limiting domains in each polymer
        for domain in limiting_domain_types:
            for j in range(number_of_possible_polymers):
                model.Add(
                    sum(
                        monomer.net_count(domain) * polymer_composition_vars[i, j]
                        for i, monomer in enumerate(ordered_monomer_types)
                    ) <= 0
                )

        # this constraint encodes the logic for the polymer indicator variables: 0 if polymer is empty else 0/1
        for j in range(number_of_possible_polymers):
            model.Add(indicator_vars[j] <=
                    sum(
                        polymer_composition_vars[i, j]
                        for i, monomer in enumerate(ordered_monomer_types) if monomer in limiting_monomer_types
                    )
            )


        total_number_of_monomers = sum(
            polymer_composition_vars[i, j]
                for i in range(len(ordered_monomer_types))
                for j in range(number_of_possible_polymers)
        )
        number_of_non_singleton_polymers = sum(indicator_vars[j] for j in range(number_of_possible_polymers))
        number_of_merges = total_number_of_monomers - number_of_non_singleton_polymers

        if fixed_number_of_merges is None:
            model.Minimize(number_of_merges)
        else:
            model.Add(number_of_merges == fixed_number_of_merges)

        if order_the_polymers:
            # tiebreaker variables that enforce an ordering on the polymers
            # tiebreaker_vars[i,j] = (Is the count of all monomers <= i the same in both polymers j and j-1?)
            # boundary conditions
            for j in range(number_of_possible_polymers - 1):
                model.Add(
                    tiebreaker_vars[-1, j] == int(True))  # start out with "tied" condition and then test first entry
            # general case
            for i in range(len(ordered_monomer_types)):
                for j in range(number_of_possible_polymers - 1):
                    x = polymer_composition_vars[i, j]
                    y = polymer_composition_vars[i, j + 1]

                    # case 1: ties only make sense if a tie was not already broken above
                    model.AddChainedImplication(tiebreaker_vars[i, j], tiebreaker_vars[i - 1, j])

                    # case 2: try to resolve a tie, but still tied
                    model.AddEqualToZeroImplication(tiebreaker_vars[i, j], x - y)  # x == y

                    # case 3: try to resolve a tie, and succeed
                    model.AddGreaterThanZeroImplication(
                        model.complement_var(tiebreaker_vars[i, j]),
                        tiebreaker_vars[i - 1, j],
                        x - y,  # enforce that the above conditions imply x > y
                    )

        return model, ordered_monomer_types, polymer_composition_vars

    def interpret_solution(
            self,
            tbn: Tbn,
            ordered_monomer_types: List[Monomer],
            polymer_composition_vars: Dict[Tuple[int, int], Any],
            value_dict: Dict[Any, int],
    ) -> Configuration:

        number_of_polymers = 1 + max(j for _, j in polymer_composition_vars)
        this_configuration_dict = {}

        for j in range(number_of_polymers):
            this_polymer_dict = {}
            for i, monomer in enumerate(ordered_monomer_types):
                monomer_count = value_dict[polymer_composition_vars[i, j]]
                if monomer_count > 0:
                    this_polymer_dict[monomer] = monomer_count
            if this_polymer_dict:  # not an empty polymer
                this_polymer = Polymer(this_polymer_dict)
                this_configuration_dict[this_polymer] = 1 + this_configuration_dict.get(this_polymer, 0)

        partial_configuration = Configuration(this_configuration_dict)

        difference_tbn = tbn - partial_configuration.flatten()

        for monomer_type in difference_tbn.monomer_types():
            singleton_polymer = Polymer({monomer_type: 1})
            this_configuration_dict[singleton_polymer] = \
                difference_tbn.count(monomer_type) + this_configuration_dict.get(singleton_polymer, 0)

        return Configuration(this_configuration_dict)
