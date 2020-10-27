from typing import Dict, Union
from math import inf as infinity

from source.tbn import Tbn
from source.polymer import Polymer
from source.positive_multiset import PositiveMultiset


class Configuration:
    def __init__(self, polymer_counts: Dict[Polymer, Union[int, float]]):
        for polymer in polymer_counts:
            if type(polymer) is not Polymer:
                raise AssertionError(f"created configuration with object of non-polymer type: {type(polymer)}")
        self.__polymer_counts = PositiveMultiset(Polymer, polymer_counts, allow_infinity=True)

    def number_of_polymers(self) -> int:
        return sum(self.__polymer_counts.values())

    def number_of_merges(self):
        return \
            sum(
                self.__polymer_counts[polymer] * (polymer.size() - 1)
                    for polymer in self.__polymer_counts if polymer.size() > 1
            )

    def __str__(self) -> str:
        return self.full_str(singletons=False)

    def full_str(self, singletons=True):
        polymer_strings_as_list = []
        for polymer in sorted(self.__polymer_counts.keys()):
            if singletons or polymer.size() > 1:
                count = self.__polymer_counts[polymer]
                if count > 1:
                    polymer_as_string = f"{count}{polymer}"
                else:
                    polymer_as_string = str(polymer)

                polymer_strings_as_list.append(polymer_as_string)
        polymers_as_string = "; ".join(polymer_strings_as_list)
        return f"{polymers_as_string}"

    def __hash__(self) -> int:
        return hash(self.full_str())

    def __eq__(self, other: "Configuration") -> bool:
        return self.full_str() == other.full_str()

    def flatten(self) -> Tbn:
        monomer_counts = {}
        for polymer, polymer_count in self.__polymer_counts.items():
            for monomer, monomer_count in polymer.items():
                monomer_counts[monomer] = (polymer_count * monomer_count) + monomer_counts.get(monomer, 0)
        return Tbn(monomer_counts)

    def energy(self, bond_weighting_factor: float) -> float:
        if bond_weighting_factor <= 0.0:
            raise AssertionError("For low-W formulation, must supply positive bond weighting factor.")
        limiting_domain_types = self.flatten().limiting_domain_types()
        bond_deficit = 0
        for limiting_domain_type in limiting_domain_types:
            for polymer, polymer_count in self.__polymer_counts.items():
                local_deficit = polymer_count * sum(
                    monomer.net_count(limiting_domain_type) * monomer_count
                        for monomer, monomer_count in polymer.items()
                )
                if local_deficit > 0:
                    bond_deficit += local_deficit

        return round(bond_weighting_factor * bond_deficit + self.number_of_merges(), 8)
