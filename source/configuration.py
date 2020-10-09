from typing import Dict, Union
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
                    for polymer in self.__polymer_counts
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
