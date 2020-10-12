from typing import Dict
from source.monomer import Monomer
from source.positive_multiset import PositiveMultiset


class Polymer:
    def __init__(self, monomer_counts: Dict[Monomer, int]):
        if not monomer_counts:
            raise AssertionError("received request to create empty polymer")

        self.__monomer_counts = PositiveMultiset(Monomer, monomer_counts)

    def size(self) -> int:
        return sum(self.__monomer_counts.values())

    def __str__(self) -> str:
        monomer_strings_as_list = []
        for monomer in sorted(self.__monomer_counts.keys()):
            count = self.__monomer_counts[monomer]
            if count > 1:
                monomer_as_string = f"{count}({monomer})"
            else:
                monomer_as_string = str(monomer)

            monomer_strings_as_list.append(monomer_as_string)
        monomers_as_string = ", ".join(monomer_strings_as_list)
        return f"{{{monomers_as_string}}}"

    def __hash__(self) -> int:
        return hash(str(self))

    def __lt__(self, other: "Polymer") -> bool:
        all_keys = set(self.__monomer_counts.keys()).union(other.__monomer_counts.keys())
        for monomer in sorted(all_keys):
            if self.__monomer_counts.get(monomer, 0) > other.__monomer_counts.get(monomer, 0):
                return True
            elif other.__monomer_counts.get(monomer, 0) > self.__monomer_counts.get(monomer, 0):
                return False
        else:
            return False

    def __eq__(self, other: "Polymer") -> bool:
        all_keys = set(self.__monomer_counts.keys()).union(other.__monomer_counts.keys())
        for monomer in sorted(all_keys):
            if self.__monomer_counts.get(monomer, 0) != other.__monomer_counts.get(monomer, 0):
                return False
        else:
            return True

    def items(self):
        return self.__monomer_counts.items()
