import itertools
from typing import List, Any, Dict
from math import inf as infinity

from source.formulations.abstract import Formulation as AbstractFormulation
from source.configuration import Configuration
from source.polymer import Polymer


class Formulation(AbstractFormulation):
    def _populate_model(self) -> None:
        """
        populates self.model with the variables and constraints needed to solve a formulation
        """
        self._construct_lists_and_calculate_constants()
        self._run_asserts()
        self._add_variables()
        self._add_constraints()
        self._apply_counting_constraints()
        if self.user_constraints.optimize():
            self._apply_objective_function()

    def _construct_lists_and_calculate_constants(self):
        self.ordered_monomers = list(self.tbn.monomer_types(flatten=True))
        self.total_number_of_monomers = len(self.ordered_monomers)
        self.limiting_domain_types = list(self.tbn.limiting_domain_types())

        self.model.set_big_m(self.total_number_of_monomers)

    def _add_variables(self) -> None:
        # monomer_grouping_vars[i, j] = 1 if monomers i and j are grouped into the same polymer, 0 else
        self.grouping_vars = {}
        for i in range(self.total_number_of_monomers):
            self.grouping_vars[i, i] = 1  # enforce reflexivity
            for j in range(i + 1, self.total_number_of_monomers):
                self.grouping_vars[i, j] = self.model.bool_var(f'monomer_grouping_{i}_{j}')
                self.grouping_vars[j, i] = self.grouping_vars[i, j]  # enforce symmetry

        # rep_vars[i] = 1 if monomer i is the "leader" of its polymer, for counting purpose.  Else 0.
        self.rep_vars = {}
        for i in range(self.total_number_of_monomers):
            self.rep_vars[i] = self.model.bool_var(f'rep_{i}')

    def _add_constraints(self) -> None:
        # transitivity of grouping
        for m1, m2, m3 in itertools.permutations(range(self.total_number_of_monomers), 3):
            self.model.add_implication(
                self.grouping_vars[m1, m2],
                self.grouping_vars[m1, m3],
                self.grouping_vars[m2, m3]
            )

        # cannot be two representatives in the same polymer
        for i in range(self.total_number_of_monomers):
            for j in range(i + 1, self.total_number_of_monomers):
                self.model.add_implication(self.grouping_vars[i, j], self.model.complement_var(self.rep_vars[j]))

        self._add_saturation_constraints()

    def _add_saturation_constraints(self) -> None:
        # saturation constraint: limiting sites must be in the minority in any polymer
        for i in range(self.total_number_of_monomers):
            for domain in self.limiting_domain_types:
                self.model.add_constraint(
                    sum(self.grouping_vars[i, j] * self.ordered_monomers[j].net_count(domain)
                        for j in range(self.total_number_of_monomers)) <= 0
                )

    def _apply_counting_constraints(self) -> None:
        self.number_of_polymers = sum(self.rep_vars[i] for i in range(self.total_number_of_monomers))
        if self.user_constraints.max_polymers() != infinity:
            self.model.add_constraint(self.number_of_polymers <= self.user_constraints.max_polymers())
        if self.user_constraints.min_polymers() > 0:
            self.model.add_constraint(self.number_of_polymers >= self.user_constraints.min_polymers())

    def _apply_objective_function(self) -> None:
        self.model.maximize(self.number_of_polymers)

    def _variables_to_keep(self) -> List[Any]:
        return list(self.grouping_vars.values())

    def _interpret_solution(self, variable_to_value_dictionary: Dict[Any, int]) -> Configuration:
        """
        uses the provided dictionary to convert solution variables into solution values and from this,
          converts the solutions values into the corresponding configuration
        """
        this_configuration_dict = {}
        n = len(self.ordered_monomers)
        discovered = [False for _ in range(n)]
        for i in range(n):
            if not discovered[i]:
                discovered[i] = True
                this_polymer_dict = {self.ordered_monomers[i]: 1}
                for j in range(i+1, n):
                    if variable_to_value_dictionary[self.grouping_vars[i, j]] > 0:  # Monomers i and j are bound
                        discovered[j] = True
                        this_polymer_dict[self.ordered_monomers[j]] =\
                            1 + this_polymer_dict.get(self.ordered_monomers[j], 0)

                this_polymer = Polymer(this_polymer_dict)
                this_configuration_dict[this_polymer] = 1 + this_configuration_dict.get(this_polymer, 0)

        return Configuration(this_configuration_dict)

    def _run_asserts(self) -> None:
        pass
