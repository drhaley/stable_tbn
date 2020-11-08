from math import inf as infinity
from math import ceil, floor

from source.formulations.polymer_unbounded_matrix import Formulation as UnboundedFormulation


class Formulation(UnboundedFormulation):
    def _add_variables(self) -> None:
        super()._add_variables()

        # additional variables for low-W formulations which tracks number of unbound limiting sites
        self.bond_deficit_exists_vars = {}
        self.bond_deficit_amount_vars = {}
        for domain in self.limiting_domain_types:
            total_count_of_limiting_site = sum(
                self.tbn.count(monomer_type) * abs(monomer_type.net_count(domain))
                for monomer_type in self.limiting_monomer_types
            )
            for j in range(self.max_polymers):
                # must know boolean information on whether any unbound limiting sites exist
                self.bond_deficit_exists_vars[domain, j] = self.model.bool_var(
                    f'deficit_exists_{domain}_{j}'
                )
                # and if so, the number of these limiting sites
                self.bond_deficit_amount_vars[domain, j] = self.model.int_var(
                    0,
                    total_count_of_limiting_site,
                    f'deficit_amount_{domain}_{j}'
                )

    def _add_saturation_constraints(self) -> None:
        # must saturate the limiting domains in each polymer OR pay a deficit in bond weight
        for domain in self.limiting_domain_types:
            for j in range(self.max_polymers):
                # net_count <= deficit_amount
                net_count_of_unbound_limiting_domains = sum(
                    monomer.net_count(domain) * self.polymer_composition_vars[i, j]
                    for i, monomer in enumerate(self.ordered_monomer_types)
                )
                self.model.add_constraint(
                     net_count_of_unbound_limiting_domains <= self.bond_deficit_amount_vars[domain, j]
                )
                # but with only the above, the deficit variables are not tightly constrained,
                #  and the algorithm could "take a hit" to increase energy to meet energy/polymer goals
                #  specified by the user, which in practice creates several isomorphic solutions
                #  for non-optimal energy configurations.
                #  Thus, we add the following constraint to keep the deficit variables tight

                # if any deficit exists, then net_count == deficit_amount
                self.model.add_equal_to_zero_implication(
                    self.bond_deficit_exists_vars[domain, j],
                    net_count_of_unbound_limiting_domains - self.bond_deficit_amount_vars[domain, j]
                )
                # if no deficit exists, then deficit_amount is zero
                self.model.add_equal_to_zero_implication(
                    self.model.complement_var(self.bond_deficit_exists_vars[domain, j]),
                    self.bond_deficit_amount_vars[domain, j]
                )
                # if deficit exists, then deficit_amount is greater than zero
                self.model.add_greater_than_zero_implication(
                    self.bond_deficit_exists_vars[domain, j],
                    self.bond_deficit_amount_vars[domain, j]
                )
                # TODO: find a way to do the above with fewer constraints

    def _apply_counting_constraints(self) -> None:
        super()._apply_counting_constraints()
        self.total_bond_deficit = sum(
            self.bond_deficit_amount_vars[domain, j]
            for domain in self.limiting_domain_types
            for j in range(self.max_polymers)
        )
        scaling_factor = 100
        self.scaled_energy = (
                round(scaling_factor * self.user_constraints.bond_weight()) * self.total_bond_deficit
                + scaling_factor * self.number_of_merges
        )
        if self.user_constraints.max_energy() != infinity:
            self.model.add_constraint(self.scaled_energy <= ceil(scaling_factor * self.user_constraints.max_energy()))
        if self.user_constraints.min_energy() != -infinity:
            self.model.add_constraint(self.scaled_energy >= floor(scaling_factor * self.user_constraints.min_energy()))

    def _apply_objective_function(self) -> None:
        self.model.minimize(self.scaled_energy)

    def _run_asserts(self) -> None:
        if self.user_constraints.bond_weight() is None or self.user_constraints.bond_weight() <= 0.0:
            raise AssertionError("For low-W formulation, must supply positive bond weighting factor.")
