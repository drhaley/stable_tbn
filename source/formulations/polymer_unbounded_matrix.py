from typing import List, Tuple, Dict, Any
from math import inf as infinity

from source.formulations.abstract import Formulation as AbstractFormulation
from source.monomer import Monomer
from source.polymer import Polymer
from source.configuration import Configuration


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        self._construct_lists_and_calculate_constants()
        self._run_asserts()
        self._add_variables()
        self._add_constraints()
        self._add_sorting_constraints()
        self._apply_counting_constraints()
        if self.user_constraints.optimize():
            self._apply_objective_function()

    def _construct_lists_and_calculate_constants(self) -> None:
        self.ordered_monomer_types, self.monomer_counts = self._get_monomer_types_and_counts()
        self.total_number_of_monomers = sum(
            self.tbn.count(monomer_type)
            for monomer_type in self.ordered_monomer_types
        )
        self.limiting_domain_types = list(self.tbn.limiting_domain_types())
        self.limiting_monomer_types = self._get_limiting_monomer_types()
        self.total_number_of_limiting_monomers = sum(
            self.tbn.count(monomer_type)
            for monomer_type in self.limiting_monomer_types
        )
        # upper bound on how many total monomers can be in non-singleton polymers
        self.upper_bound_on_total_monomers_in_complexes = sum(
            self.tbn.count(monomer_type) * (1 + abs(monomer_type.net_count(domain_type)))
            for monomer_type in self.limiting_monomer_types
            for domain_type in self.limiting_domain_types
        )
        self.upper_bound_on_total_monomers_in_complexes = min(
            self.upper_bound_on_total_monomers_in_complexes,
            self.total_number_of_monomers
        )
        self.model.set_big_m(self.upper_bound_on_total_monomers_in_complexes)
        if self.user_constraints.max_polymers() == infinity:
            self.max_polymers = self.total_number_of_limiting_monomers
        else:
            self.max_polymers = self.user_constraints.max_polymers()

    def _add_variables(self) -> None:
        # polymer_composition_vars[i, j] = count of monomer i in polymer j
        self.polymer_composition_vars = {}
        for i, monomer_type in enumerate(self.ordered_monomer_types):
            for j in range(self.max_polymers):
                self.polymer_composition_vars[i, j] = self.model.int_var(
                    0,
                    min(self.monomer_counts[i], self.upper_bound_on_total_monomers_in_complexes),
                    f'polymer_composition_{i}_{j}'
                )

        # indicator variables for polymers
        self.indicator_vars = {}
        for j in range(self.max_polymers):
            self.indicator_vars[j] = self.model.bool_var(f'indicator_{j}')

    def _add_constraints(self) -> None:
        self._add_conservation_constraints()
        self._add_saturation_constraints()
        self._add_polymer_indicator_constraints()

    def _add_conservation_constraints(self) -> None:
        # monomer conservation; must use all limiting monomers, and cannot exceed the count of other monomers
        for i, monomer in enumerate(self.ordered_monomer_types):
            if monomer in self.limiting_monomer_types:
                self.model.add_constraint(
                    sum(
                        self.polymer_composition_vars[i, j]
                        for j in range(self.max_polymers)
                    ) == self.monomer_counts[i]
                )
            elif self.monomer_counts[i] < infinity:
                self.model.add_constraint(
                    sum(
                        self.polymer_composition_vars[i, j]
                        for j in range(self.max_polymers)
                    ) <= self.monomer_counts[i]
                )

    def _add_saturation_constraints(self) -> None:
        # must saturate the limiting domains in each polymer
        for domain in self.limiting_domain_types:
            for j in range(self.max_polymers):
                self.model.add_constraint(
                    sum(
                        monomer.net_count(domain) * self.polymer_composition_vars[i, j]
                        for i, monomer in enumerate(self.ordered_monomer_types)
                    ) <= 0
                )

    def _add_polymer_indicator_constraints(self) -> None:
        # this constraint encodes the logic for the polymer indicator variables: 0 if polymer is empty else 0/1
        #  can only be 1 if the polymer has a limiting monomer inside it,
        #  to prevent the algorithm for producing "optimal" configurations in which
        #  superfluous singletons are made explicit, which results in many isomorphic solutions
        for j in range(self.max_polymers):
            self.model.add_constraint(
                self.indicator_vars[j] <=
                sum(self.polymer_composition_vars[i, j]
                    for i, monomer in enumerate(self.ordered_monomer_types)
                    if monomer in self.limiting_monomer_types)
            )

    def _add_sorting_constraints(self) -> None:
        if self.user_constraints.sort():
            # tiebreaker variables that enforce lexicographical ordering of polymers based upon monomer counts within
            # tiebreaker_vars[i,j] = (Is the count of all monomers <= i the same in both polymers j and j-1?)
            tiebreaker_vars = {}
            for i in range(-1, len(self.ordered_monomer_types)):
                for j in range(self.max_polymers - 1):
                    tiebreaker_vars[i, j] = self.model.bool_var(f'tiebreaker_{i}_{j}')

            # boundary conditions
            for j in range(self.max_polymers - 1):
                self.model.add_constraint(
                    tiebreaker_vars[-1, j] == int(
                        True))  # start out with "tied" condition and then test first entry
            # general case
            for i in range(len(self.ordered_monomer_types)):
                for j in range(self.max_polymers - 1):
                    x = self.polymer_composition_vars[i, j]
                    y = self.polymer_composition_vars[i, j + 1]

                    # case 1: ties only make sense if a tie was not already broken above
                    self.model.add_implication(tiebreaker_vars[i, j], tiebreaker_vars[i - 1, j])

                    # case 2: try to resolve a tie, but still tied
                    self.model.add_equal_to_zero_implication(tiebreaker_vars[i, j], x - y)  # x == y

                    # case 3: try to resolve a tie, and succeed
                    self.model.add_greater_than_zero_implication(
                        self.model.complement_var(tiebreaker_vars[i, j]),
                        tiebreaker_vars[i - 1, j],
                        x - y,  # enforce that the above conditions imply x > y
                    )

    def _apply_counting_constraints(self) -> None:
        total_number_of_monomers_used = sum(
            self.polymer_composition_vars[i, j]
                for i in range(len(self.ordered_monomer_types))
                for j in range(self.max_polymers)
        )
        number_of_polymers = sum(self.indicator_vars[j] for j in range(self.max_polymers))
        self.number_of_merges = total_number_of_monomers_used - number_of_polymers

        if self.user_constraints.max_polymers() != infinity:
            self.model.add_constraint(number_of_polymers <= self.user_constraints.max_polymers())
        if self.user_constraints.min_polymers() > 0:
            self.model.add_constraint(number_of_polymers >= self.user_constraints.min_polymers())
        if self.user_constraints.max_merges() != infinity:
            self.model.add_constraint(self.number_of_merges <= self.user_constraints.max_merges())
        if self.user_constraints.min_merges() > 0:
            self.model.add_constraint(self.number_of_merges >= self.user_constraints.min_merges())

    def _apply_objective_function(self) -> None:
        self.model.minimize(self.number_of_merges)

    def _variables_to_keep(self) -> List[Any]:
        """
        returns a list of the variables that are necessary to convert a solution back to a configuration
          (e.g. returns polymer composition variables but not 'internal' tie-breaker variables)
        """
        return list(self.polymer_composition_vars.values())

    def _interpret_solution(self, variable_to_value_dictionary: Dict[Any, int]) -> Configuration:
        """
        uses the provided dictionary to convert solution variables into solution values and from this,
          converts the solutions values into the corresponding configuration
        """
        number_of_polymers = 1 + max(j for _, j in self.polymer_composition_vars)
        this_configuration_dict = {}

        for j in range(number_of_polymers):
            this_polymer_dict = {}
            for i, monomer in enumerate(self.ordered_monomer_types):
                monomer_count = variable_to_value_dictionary[self.polymer_composition_vars[i, j]]
                if monomer_count > 0:
                    this_polymer_dict[monomer] = monomer_count + this_polymer_dict.get(monomer, 0)
            if this_polymer_dict:  # not an empty polymer
                this_polymer = Polymer(this_polymer_dict)
                this_configuration_dict[this_polymer] = 1 + this_configuration_dict.get(this_polymer, 0)

        partial_configuration = Configuration(this_configuration_dict)

        difference_tbn = self.tbn - partial_configuration.flatten()

        for monomer_type in difference_tbn.monomer_types():
            singleton_polymer = Polymer({monomer_type: 1})
            this_configuration_dict[singleton_polymer] = \
                difference_tbn.count(monomer_type) + this_configuration_dict.get(singleton_polymer, 0)

        return Configuration(this_configuration_dict)

    def _run_asserts(self) -> None:
        pass

    def _get_monomer_types_and_counts(self) -> Tuple[List[Monomer], List[int]]:
        ordered_monomer_types = list(self.tbn.monomer_types())
        monomer_counts = [self.tbn.count(monomer) for monomer in ordered_monomer_types]
        return ordered_monomer_types, monomer_counts

    def _get_limiting_monomer_types(self) -> List[Monomer]:
        return list(self.tbn.limiting_monomer_types())
