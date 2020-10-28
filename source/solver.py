import itertools

from typing import Iterator, Optional, Tuple, List, Dict, Any
from enum import Enum, auto
from math import inf as infinity
from math import ceil

from source.tbn import Tbn
from source.configuration import Configuration
from source.polymer import Polymer
from source.monomer import Monomer
from source.solver_adapters import constraint_programming, integer_programming
from source.solver_adapters import abstract as adapter


class SolverMethod(Enum):
    CONSTRAINT_PROGRAMMING = auto()
    INTEGER_PROGRAMMING = auto()  # Only implemented for single queries


class SolverFormulation(Enum):
    STABLEGEN_FORMULATION = auto()  # TODO: NYI
    BOND_OBLIVIOUS_FORMULATION = auto()
    SET_FORMULATION = auto()
    MULTISET_FORMULATION = auto()
    BEYOND_MULTISET_FORMULATION = auto()
    LOW_W_FORMULATION = auto()


class Solver:
    def __init__(
            self,
            method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
    ):
        if method == SolverMethod.CONSTRAINT_PROGRAMMING:
            self.__adapter = constraint_programming.Solver()
            self.__sorted_polymers = True
        elif method == SolverMethod.INTEGER_PROGRAMMING:
            self.__adapter = integer_programming.Solver()
            self.__sorted_polymers = False
        else:
            raise NotImplementedError(f"solver not implemented for method {method}")

    def stable_config(self,
                      tbn: Tbn,
                      formulation: SolverFormulation = SolverFormulation.BEYOND_MULTISET_FORMULATION,
                      bond_weighting_factor: Optional[float] = None,
                      ) -> Configuration:
        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(
            tbn, formulation=formulation, bond_weighting_factor=bond_weighting_factor
        )
        status = self.__adapter.solve(model, list(polymer_composition_vars.values()))
        if status == model.INFEASIBLE:
            raise AssertionError(f"Could not find optimal solution to tbn, was reported infeasible")
        elif status != model.OPTIMAL:
            raise AssertionError(f"could not find optimal solution to tbn, got code {status} instead")
        else:
            value_dict = {var: self.__adapter.value(var) for var in polymer_composition_vars.values()}
            return self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, value_dict, formulation=formulation)

    def stable_configs(self,
                       tbn: Tbn,
                       formulation: SolverFormulation = SolverFormulation.BEYOND_MULTISET_FORMULATION,
                       bond_weighting_factor: Optional[float] = None,
                       ) -> Iterator[Configuration]:
        example_stable_configuration = self.stable_config(
            tbn, formulation=formulation, bond_weighting_factor=bond_weighting_factor
        )
        if formulation == SolverFormulation.LOW_W_FORMULATION:
            max_energy = example_stable_configuration.energy(bond_weighting_factor)
            return self.configs_with_energy(tbn, formulation=formulation, max_energy=max_energy, bond_weighting_factor=bond_weighting_factor)
        elif formulation == SolverFormulation.BEYOND_MULTISET_FORMULATION:
            number_of_merges = example_stable_configuration.number_of_merges()
            return self.configs_with_number_of_merges(tbn, formulation=formulation, number_of_merges=number_of_merges)
        else:
            number_of_polymers = example_stable_configuration.number_of_polymers()
            return self.configs_with_number_of_polymers(tbn, formulation=formulation, number_of_polymers=number_of_polymers)

    def configs_with_number_of_polymers(self,
                                        tbn: Tbn,
                                        formulation: SolverFormulation,
                                        number_of_polymers: Optional[int] = None,
                                        ) -> Iterator[Configuration]:
        total_number_of_monomers = sum(
            tbn.count(monomer_type)
            for monomer_type in tbn.monomer_types()
        )
        if total_number_of_monomers == infinity:
            raise AssertionError(
                "cannot request a specific number of polymers when there are an infinite amount of monomers"
            )
        if formulation in [SolverFormulation.BEYOND_MULTISET_FORMULATION, SolverFormulation.LOW_W_FORMULATION]:
            raise AssertionError("to specify a specific number of polymers with this formulation, " +
                                 "must input an appropriate constraint file")

        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(
            tbn, formulation=formulation, a_priori_number_of_polymers=number_of_polymers,
        )
        solutions = self.__adapter.solve_all(model, list(polymer_composition_vars.values()))
        return [
            self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, solution, formulation=formulation)
            for solution in solutions
        ]

    def configs_with_number_of_merges(self,
                                      tbn: Tbn,
                                      formulation: SolverFormulation,
                                      number_of_merges: Optional[int] = None,
                                      max_number_of_polymers: Optional[int] = None,
                                      ) -> Iterator[Configuration]:
        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(
            tbn,
            formulation=formulation,
            max_number_of_merges=number_of_merges,
            a_priori_number_of_polymers=max_number_of_polymers,
        )
        solutions = self.__adapter.solve_all(model, list(polymer_composition_vars.values()))
        return [
            self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, solution, formulation=formulation)
            for solution in solutions
        ]

    def configs_with_energy(self,
                            tbn: Tbn,
                            formulation: SolverFormulation,
                            max_energy: Optional[float] = None,
                            bond_weighting_factor: Optional[float] = None,
                            max_number_of_polymers: Optional[int] = None,
    ) -> Iterator[Configuration]:
        model, ordered_monomer_types, polymer_composition_vars = self.__build_model(
            tbn,
            formulation=formulation,
            max_energy=max_energy,
            bond_weighting_factor=bond_weighting_factor,
            a_priori_number_of_polymers=max_number_of_polymers,
        )
        solutions = self.__adapter.solve_all(model, list(polymer_composition_vars.values()))
        return [
            self.interpret_solution(tbn, ordered_monomer_types, polymer_composition_vars, solution, formulation=formulation)
            for solution in solutions
        ]

    def __build_model(self,
                      tbn: Tbn,
                      formulation: SolverFormulation,
                      a_priori_number_of_polymers: Optional[int] = None,
                      max_number_of_merges: Optional[int] = None,
                      max_energy: Optional[float] = None,
                      bond_weighting_factor: Optional[float] = None,
                      ) -> Tuple[adapter.Model, List[Monomer], Dict[Tuple[int, int], Any]]:
        """
        :param tbn:
        :param formulation:
        Specifies the type of programming formulation to be used (from the SolverFormulation enum)
        :param a_priori_number_of_polymers:
        For all formulations except BEYOND_MULTISET_FORMULATION and LOW_W_FORMULATION,
         will supply all the configurations with this number of polymers, otherwise supplies a single configuration
          with the most number of polymers possible
        If using STABLEGEN_FORMULATION or BOND_OBLIVIOUS FORMULATION, this instead acts as the lower bound
         on the polymer count, since those formulations are not constrained by a polymer composition matrix,
         instead they are constrained by the counting method of Rep() variables, which can undercount
        If using BEYOND_MULTISET_FORMULATION or LOW_W_FORMULATION, optionally specifies a priori knowledge
         of the upper bound of the number of non-singleton polymers (to improve efficiency)
        :param max_number_of_merges:
        Only used for BEYOND_MULTISET_FORMULATION.
        If specified, the model will supply all the configurations with this number of merges or fewer, otherwise
         supplies a single configuration with the fewest number of merges possible
        :param max_energy:
        Only used for LOW_W_FORMULATION.
        If specified, the model will supply all the configurations with this energy or less, otherwise
         supplies a single configuration with the lowest energy possible
        :return:
            (model, ordered_monomer_list, polymer_inclusion_matrix_entry_dict)
        """

        # StableGen formulation is qualitatively different, so factored that code out.  Some code is shared, but
        #  seems more straightforward at this point to isolate them
        if formulation in [SolverFormulation.STABLEGEN_FORMULATION, SolverFormulation.BOND_OBLIVIOUS_FORMULATION]:
            return self.__build_stablegen_model(tbn, formulation, a_priori_number_of_polymers)

        if formulation == SolverFormulation.LOW_W_FORMULATION:
            if bond_weighting_factor is None or bond_weighting_factor <= 0.0:
                raise AssertionError("For low-W formulation, must supply positive bond weighting factor.")
            if max_energy is None:
                optimize = True
            else:
                optimize = False
        elif formulation == SolverFormulation.BEYOND_MULTISET_FORMULATION:
            if max_number_of_merges is None:
                optimize = True
            else:
                optimize = False
        else:
            if a_priori_number_of_polymers is None:
                optimize = True
            else:
                optimize = False

        model = self.__adapter.model()

        if formulation == SolverFormulation.SET_FORMULATION:
            ordered_monomer_types = list(tbn.monomer_types(flatten=True))
            monomer_counts = [1 for _ in ordered_monomer_types]
        else:
            ordered_monomer_types = list(tbn.monomer_types())
            monomer_counts = [tbn.count(monomer) for monomer in ordered_monomer_types]

        total_number_of_monomers = sum(
            tbn.count(monomer_type)
            for monomer_type in ordered_monomer_types
        )

        limiting_domain_types = list(tbn.limiting_domain_types())

        limiting_monomer_types = list(tbn.limiting_monomer_types())
        total_number_of_limiting_monomers = sum(
            tbn.count(monomer_type)
            for monomer_type in limiting_monomer_types
        )

        # upper bound on how many total monomers can be in non-singleton polymers
        upper_bound_on_total_monomers_in_complexes = sum(
            tbn.count(monomer_type) * (1 + abs(monomer_type.net_count(domain_type)))
            for monomer_type in limiting_monomer_types
            for domain_type in limiting_domain_types
        )
        upper_bound_on_total_monomers_in_complexes = min(
            upper_bound_on_total_monomers_in_complexes,
            total_number_of_monomers
        )

        # set big M, in case integer programming algorithm is called for instead of constraint programming
        #  this value must be at least as large as the maximum number of (tracked) monomers that can be in a polymer
        if formulation in {SolverFormulation.BEYOND_MULTISET_FORMULATION, SolverFormulation.LOW_W_FORMULATION}:
            model.set_big_M(upper_bound_on_total_monomers_in_complexes)
        else:
            model.set_big_M(total_number_of_monomers)

        if a_priori_number_of_polymers is None:
            if formulation in {SolverFormulation.BEYOND_MULTISET_FORMULATION, SolverFormulation.LOW_W_FORMULATION}:
                a_priori_number_of_polymers = total_number_of_limiting_monomers
            else:
                a_priori_number_of_polymers = total_number_of_monomers

        ############ Variables ############
        # polymer_composition_vars[i, j] = count of monomer i in polymer j
        polymer_composition_vars = {}
        for i, monomer_type in enumerate(ordered_monomer_types):
            for j in range(a_priori_number_of_polymers):
                polymer_composition_vars[i, j] = model.int_var(
                    0,
                    min(monomer_counts[i], upper_bound_on_total_monomers_in_complexes),  # latter used if count is inf
                    f'polymer_composition_{i}_{j}'
                )

        # indicator variables for polymers
        indicator_vars = {}
        for j in range(a_priori_number_of_polymers):
            indicator_vars[j] = model.bool_var(f'indicator_{j}')

        # additional variables for low-W formulations which tracks number of unbound sites
        if formulation == SolverFormulation.LOW_W_FORMULATION:
            bond_deficit_vars = {}
            for domain in limiting_domain_types:
                total_count_of_limiting_site = sum(
                    tbn.count(monomer_type) * abs(monomer_type.net_count(domain))
                    for monomer_type in limiting_monomer_types
                )
                for j in range(a_priori_number_of_polymers):
                    bond_deficit_vars[domain, j] = model.int_var(
                        0,
                        total_count_of_limiting_site,
                        f'deficit_{domain}_{j}'
                    )

        ############ Constraints ############
        # monomer conservation; must use all limiting monomers, and cannot exceed the count of other monomers
        # if not using BEYOND_MULTISET_FORMULATION or LOW_W_FORMULATION, all monomers must be used in exact counts
        for i, monomer in enumerate(ordered_monomer_types):
            if (formulation not in [SolverFormulation.BEYOND_MULTISET_FORMULATION, SolverFormulation.LOW_W_FORMULATION])\
                    or (monomer in limiting_monomer_types):
                model.Add(
                    sum(
                        polymer_composition_vars[i, j]
                        for j in range(a_priori_number_of_polymers)
                    ) == monomer_counts[i]
                )
            elif monomer_counts[i] < infinity:
                model.Add(
                    sum(
                        polymer_composition_vars[i, j]
                        for j in range(a_priori_number_of_polymers)
                    ) <= monomer_counts[i]
                )

        # must saturate the limiting domains in each polymer
        for domain in limiting_domain_types:
            for j in range(a_priori_number_of_polymers):
                if formulation != SolverFormulation.LOW_W_FORMULATION:
                    deficit = 0
                else:
                    deficit = bond_deficit_vars[domain, j]
                model.Add(
                    sum(
                        monomer.net_count(domain) * polymer_composition_vars[i, j]
                        for i, monomer in enumerate(ordered_monomer_types)
                    ) <= deficit
                )

        # this constraint encodes the logic for the polymer indicator variables: 0 if polymer is empty else 0/1
        # if using BEYOND_MULTISET_FORMULATION or LOW_W_FORMULATION, can only be 1 if the polymer has
        #  a limiting monomer inside it,
        #  to prevent the algorithm for producing "optimal" configurations in which
        #  superfluous singletons are made explicit, which results in many isomorphic solutions
        for j in range(a_priori_number_of_polymers):
            if formulation in {SolverFormulation.BEYOND_MULTISET_FORMULATION, SolverFormulation.LOW_W_FORMULATION}:
                model.Add(indicator_vars[j] <=
                          sum(polymer_composition_vars[i, j]
                              for i, monomer in enumerate(ordered_monomer_types) if monomer in limiting_monomer_types))
            else:
                model.Add(indicator_vars[j] <=
                          sum(polymer_composition_vars[i, j]
                              for i, _ in enumerate(ordered_monomer_types)))

        # these aggregate expressions aid in readability of the code, but are not explicit variables or constraints
        total_number_of_monomers_used = sum(
            polymer_composition_vars[i, j]
            for i in range(len(ordered_monomer_types))
            for j in range(a_priori_number_of_polymers)
        )
        number_of_polymers = sum(indicator_vars[j] for j in range(a_priori_number_of_polymers))
        number_of_merges = total_number_of_monomers_used - number_of_polymers

        if formulation == SolverFormulation.LOW_W_FORMULATION:
            total_bond_deficit = sum(
                bond_deficit_vars[domain, j]
                    for domain in limiting_domain_types
                    for j in range(a_priori_number_of_polymers)
            )
            SCALING_FACTOR = 100
            energy = (
                    round(SCALING_FACTOR * bond_weighting_factor) * total_bond_deficit \
                    + SCALING_FACTOR * number_of_merges
            )
            if optimize:
                model.Minimize(energy)
            else:
                model.Add(energy <= ceil(SCALING_FACTOR * max_energy))
        elif formulation == SolverFormulation.BEYOND_MULTISET_FORMULATION:
            if optimize:
                model.Minimize(number_of_merges)
            else:
                model.Add(number_of_merges <= max_number_of_merges)
        else:
            if optimize:
                model.Maximize(number_of_polymers)
            else:
                model.Add(number_of_polymers == a_priori_number_of_polymers)

        if self.__sorted_polymers:
            # tiebreaker variables that enforce lexicographical ordering of polymers based upon monomer counts within
            # tiebreaker_vars[i,j] = (Is the count of all monomers <= i the same in both polymers j and j-1?)
            tiebreaker_vars = {}
            for i in range(-1, len(ordered_monomer_types)):
                for j in range(a_priori_number_of_polymers - 1):
                    tiebreaker_vars[i, j] = model.bool_var(f'tiebreaker_{i}_{j}')

            # boundary conditions
            for j in range(a_priori_number_of_polymers - 1):
                model.Add(
                    tiebreaker_vars[-1, j] == int(True))  # start out with "tied" condition and then test first entry
            # general case
            for i in range(len(ordered_monomer_types)):
                for j in range(a_priori_number_of_polymers - 1):
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

    def __build_stablegen_model(self,
                                tbn: Tbn,
                                formulation: SolverFormulation,
                                a_priori_number_of_polymers: Optional[int] = None,
                                ) -> Tuple[adapter.Model, List[Monomer], Dict[Tuple[int, int], Any]]:

        if a_priori_number_of_polymers is None:
            optimize = True
        else:
            optimize = False

        model = self.__adapter.model()

        ordered_monomers = list(tbn.monomer_types(flatten=True))
        total_number_of_monomers = len(ordered_monomers)
        limiting_domain_types = list(tbn.limiting_domain_types())

        # set big M, in case integer programming algorithm is called for instead of constraint programming
        #  this value must be at least as large as the maximum number of (tracked) monomers that can be in a polymer
        model.set_big_M(total_number_of_monomers)

        ############ Variables ############
        if formulation == SolverFormulation.STABLEGEN_FORMULATION:
            pair_vars = {}
            for monomer_number, monomer in enumerate(ordered_monomers):
                for domain_number, domain_type in enumerate(monomer.as_explicit_list()):
                    for second_monomer_number, second_monomer in enumerate(ordered_monomers):
                        for second_domain_number, second_domain_type in enumerate(second_monomer.as_explicit_list()):
                            if domain_type == second_domain_type.complement():
                                if (monomer > second_monomer) or \
                                        (monomer == second_monomer and domain_number > second_domain_number):
                                    # symmetric condition
                                    new_variable_or_value = pair_vars[(second_monomer_number, second_domain_number), (monomer_number, domain_number)]
                                else:
                                    new_variable_or_value = model.bool_var(f"pair_{monomer_number}_{domain_number}_{second_monomer_number}_{second_domain_number}")
                            else:
                                new_variable_or_value = 0
                            pair_vars[(monomer_number, domain_number), (second_monomer_number, second_domain_number)] = new_variable_or_value

        # bind_vars[i, j] = 1 if monomers i and j are bound into the same polymer, 0 else
        bind_vars = {}
        for i in range(total_number_of_monomers):
            bind_vars[i, i] = 1  # enforce reflexivity
            for j in range(i+1, total_number_of_monomers):
                bind_vars[i, j] = model.bool_var(f'bind_{i}_{j}')
                bind_vars[j, i] = bind_vars[i, j]  # enforce symmetry

        # rep_vars[i] = 1 if monomer i is the "leader" of its polymer, for counting purpose.  Else 0.
        rep_vars = {}
        for i in range(total_number_of_monomers):
            rep_vars[i] = model.bool_var(f'rep_{i}')

        ############ Constraints ############
        # transitivity of binding
        for m1, m2, m3 in itertools.permutations(range(total_number_of_monomers), 3):
            model.AddChainedImplication(bind_vars[m1, m2], bind_vars[m1, m3], bind_vars[m2, m3])

        # cannot be two representatives in the same polymer
        for i in range(total_number_of_monomers):
            for j in range(i+1, total_number_of_monomers):
                model.AddChainedImplication(bind_vars[i, j], model.complement_var(rep_vars[j]))

        if formulation == SolverFormulation.STABLEGEN_FORMULATION:
            # saturation constraint: all limiting sites are paired
            for this_monomer_number, monomer in enumerate(ordered_monomers):
                for this_domain_number, domain in enumerate(monomer.as_explicit_list()):
                    if domain in limiting_domain_types:
                        model.Add(
                            1 == sum(val for (((monomer_number, domain_number), (_,_)), val) in pair_vars.items()
                                     if this_monomer_number == monomer_number and this_domain_number == domain_number)
                        )

            # cannot pair more than once
            for this_monomer_number, monomer in enumerate(ordered_monomers):
                for this_domain_number, domain in enumerate(monomer.as_explicit_list()):
                    model.Add(
                        1 >= sum(val for (((monomer_number, domain_number), (_,_)), val) in pair_vars.items()
                                 if this_monomer_number == monomer_number and this_domain_number == domain_number)
                    )

            # pair implies bond
            for ((monomer_number, domain_number), (second_monomer_number, second_domain_number)), var in pair_vars.items():
                model.AddChainedImplication(var, bind_vars[monomer_number, second_monomer_number])
        else:
            # saturation constraint: limiting sites must be in the minority in any polymer
            for i in range(total_number_of_monomers):
                for domain in limiting_domain_types:
                    model.Add(
                        sum(bind_vars[i, j] * ordered_monomers[j].net_count(domain)
                            for j in range(total_number_of_monomers)) <= 0
                    )

        # this aggregate expression aids in readability of the code, but is not an explicit variable or constraints
        number_of_polymers = sum(rep_vars[i] for i in range(total_number_of_monomers))

        if optimize:
            model.Maximize(number_of_polymers)
        else:
            model.Add(number_of_polymers == a_priori_number_of_polymers)

        return model, ordered_monomers, bind_vars

    def interpret_solution(
                self,
                tbn: Tbn,
                ordered_monomer_types: List[Monomer],
                polymer_composition_vars: Dict[Tuple[int, int], Any],
                value_dict: Dict[Any, int],
                formulation: SolverFormulation,
            ) -> Configuration:

        if formulation in [SolverFormulation.STABLEGEN_FORMULATION, SolverFormulation.BOND_OBLIVIOUS_FORMULATION]:
            return self.interpret_bind_map_solution(tbn, ordered_monomer_types, polymer_composition_vars, value_dict)
        else:
            return self.interpret_polymer_dict_solution(tbn, ordered_monomer_types, polymer_composition_vars, value_dict)

    @staticmethod
    def interpret_bind_map_solution(
                _: Tbn,
                ordered_monomers: List[Monomer],
                bind_vars: Dict[Tuple[int, int], Any],
                value_dict: Dict[Any, int],
            ) -> Configuration:

        this_configuration_dict = {}
        n = len(ordered_monomers)
        discovered = [False for _ in range(n)]
        for i in range(n):
            if not discovered[i]:
                discovered[i] = True
                this_polymer_dict = {ordered_monomers[i]: 1}
                for j in range(i+1, n):
                    if value_dict[bind_vars[i, j]] > 0:  # Monomers i and j are bound
                        discovered[j] = True
                        this_polymer_dict[ordered_monomers[j]] = 1 + this_polymer_dict.get(ordered_monomers[j], 0)

                this_polymer = Polymer(this_polymer_dict)
                this_configuration_dict[this_polymer] = 1 + this_configuration_dict.get(this_polymer, 0)

        return Configuration(this_configuration_dict)

    def interpret_polymer_dict_solution(
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
                    this_polymer_dict[monomer] = monomer_count + this_polymer_dict.get(monomer, 0)
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
