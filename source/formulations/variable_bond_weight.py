from math import inf as infinity
from math import ceil, floor

from source.formulations.polymer_unbounded_matrix import Formulation as UnboundedFormulation


class Formulation(UnboundedFormulation):
    def _add_variables(self) -> None:
        super()._add_variables()

        # additional variables for low-W formulations which tracks number of unbound sites
        self.bond_deficit_vars = {}
        for domain in self.limiting_domain_types:
            total_count_of_limiting_site = sum(
                self.tbn.count(monomer_type) * abs(monomer_type.net_count(domain))
                for monomer_type in self.limiting_monomer_types
            )
            for j in range(self.max_polymers):
                self.bond_deficit_vars[domain, j] = self.model.int_var(
                    0,
                    total_count_of_limiting_site,
                    f'deficit_{domain}_{j}'
                )

    def _add_saturation_constraints(self) -> None:
        # must saturate the limiting domains in each polymer OR pay a deficit in bond weight
        for domain in self.limiting_domain_types:
            for j in range(self.max_polymers):
                self.model.add_constraint(
                    sum(
                        monomer.net_count(domain) * self.polymer_composition_vars[i, j]
                        for i, monomer in enumerate(self.ordered_monomer_types)
                    ) <= self.bond_deficit_vars[domain, j]
                )

    def _apply_counting_constraints(self) -> None:
        super()._apply_counting_constraints()
        self.total_bond_deficit = sum(
            self.bond_deficit_vars[domain, j]
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
        super()._run_asserts()

        if self.user_constraints.bond_weight() is None or self.user_constraints.bond_weight() <= 0.0:
            raise AssertionError("For low-W formulation, must supply positive bond weighting factor.")
