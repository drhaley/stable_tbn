from source.formulations.bond_oblivious_network import Formulation as BondObliviousFormulation


class Formulation(BondObliviousFormulation):
    def _add_variables(self) -> None:
        super()._add_variables()
        self.site_bonding_vars = {}
        for monomer_number, monomer in enumerate(self.ordered_monomers):
            for domain_number, domain_type in enumerate(monomer.as_explicit_list()):
                for second_monomer_number, second_monomer in enumerate(self.ordered_monomers):
                    for second_domain_number, second_domain_type in enumerate(second_monomer.as_explicit_list()):
                        first_domain_tuple = (monomer_number, domain_number)
                        second_domain_tuple = (second_monomer_number, second_domain_number)
                        if domain_type == second_domain_type.complement():
                            if (monomer > second_monomer) or \
                                    (monomer == second_monomer and domain_number > second_domain_number):
                                # symmetric condition
                                new_variable_or_value = self.site_bonding_vars[second_domain_tuple, first_domain_tuple]
                            else:
                                new_variable_or_value = self.model.bool_var(
                                    f"domain_bind_{monomer_number}_{domain_number}" +
                                    f"_{second_monomer_number}_{second_domain_number}"
                                )
                        else:
                            new_variable_or_value = 0

                        self.site_bonding_vars[first_domain_tuple, second_domain_tuple] = new_variable_or_value

    def _add_saturation_constraints(self) -> None:
        # saturation constraint: all limiting sites are bound
        for this_monomer_number, monomer in enumerate(self.ordered_monomers):
            for this_domain_number, domain in enumerate(monomer.as_explicit_list()):
                if domain in self.limiting_domain_types:
                    self.model.add_constraint(
                        1 == sum(val for (((monomer_number, domain_number), _), val) in self.site_bonding_vars.items()
                                 if this_monomer_number == monomer_number and this_domain_number == domain_number)
                    )

        # cannot bind more than once
        for this_monomer_number, monomer in enumerate(self.ordered_monomers):
            for this_domain_number, domain in enumerate(monomer.as_explicit_list()):
                self.model.add_constraint(
                    1 >= sum(val for (((monomer_number, domain_number), _), val) in self.site_bonding_vars.items()
                             if this_monomer_number == monomer_number and this_domain_number == domain_number)
                )

        # bond implies grouping
        for ((monomer_number, domain_number), (second_monomer_number, second_domain_number)), var \
                in self.site_bonding_vars.items():
            self.model.add_implication(var, self.grouping_vars[monomer_number, second_monomer_number])
