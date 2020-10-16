import re
from math import inf as infinity
from math import isnan as isNotANumber
from typing import Dict, Iterator, Union

from source.monomer import Monomer
from source.domain import Domain
from source.positive_multiset import PositiveMultiset


class Tbn:
    def __init__(self, monomer_counts: Dict[Monomer, Union[int, float]]):
        self.__monomer_counts = PositiveMultiset(Monomer, monomer_counts, allow_infinity=True)

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

    def __lt__(self, other: "Tbn") -> bool:
        all_keys = set(self.__monomer_counts.keys()).union(other.__monomer_counts.keys())
        for monomer in sorted(all_keys):
            if self.__monomer_counts.get(monomer, 0) > other.__monomer_counts.get(monomer, 0):
                return True
            elif other.__monomer_counts.get(monomer, 0) > self.__monomer_counts.get(monomer, 0):
                return False
        else:
            return False

    def __eq__(self, other: "Tbn") -> bool:
        return self.__monomer_counts == other.__monomer_counts

    def __sub__(self, other: "Tbn") -> "Tbn":
        difference_tbn = {}
        for monomer_type in set(self.monomer_types()).union(set(other.monomer_types())):
            difference = self.count(monomer_type) - other.count(monomer_type)
            if difference > 0:
                difference_tbn[monomer_type] = difference
            elif difference < 0:
                raise AssertionError("Can only subtract Tbns when second is a subset of the first")

        return Tbn(difference_tbn)

    @classmethod
    def from_string(cls, text) -> "Tbn":
        monomer_counts = {}

        for raw_line in text.split('\n'):
            line = raw_line.strip()
            if line:  # not just whitespace
                quantity_search_result = re.match(
                    f"^(|inf|[1-9]\\d*)\\[\\s*({Monomer.regex()})\\s*\\]$", line
                )
                if quantity_search_result:
                    raw_count, monomer_string = quantity_search_result.groups()
                    if raw_count == "inf":
                        count = infinity
                    elif raw_count == "":
                        count = 1
                    else:
                        count = int(raw_count)
                else:
                    count = 1
                    monomer_string = line

                monomer = Monomer.from_string(monomer_string)
                monomer_counts[monomer] = monomer_counts.get(monomer, 0) + count

        return Tbn(monomer_counts)

    def monomer_types(self, flatten:bool = False) -> Iterator[Monomer]:
        for monomer in sorted(self.__monomer_counts):
            if flatten:
                if self.count(monomer) == infinity:
                    raise AssertionError("cannot flatten a tbn with infinite quantities of monomers")
                else:
                    for _ in range(self.count(monomer)):
                        yield monomer
            else:
                yield monomer

    def limiting_domain_types(self) -> Iterator[Domain]:
        domain_tally = {}
        for monomer in self.__monomer_counts:
            count_of_monomer = self.__monomer_counts[monomer]
            for domain in monomer.unstarred_domain_types():
                domain_tally[domain] = monomer.net_count(domain) * count_of_monomer + domain_tally.get(domain, 0)
                if isNotANumber(domain_tally[domain]):  # result of infinity - infinity
                    raise AssertionError(f"domain {domain} exists in opposing infinite quantities")

        for unstarred_domain in sorted(domain_tally.keys()):
            if domain_tally[unstarred_domain] >= 0:
                yield unstarred_domain.complement()
            else:
                yield unstarred_domain  # Note: if tied, also yields starred_domain

    def limiting_monomer_types(self) -> Iterator[Monomer]:
        limiting_domain_types = list(self.limiting_domain_types())
        for monomer_type in sorted(self.__monomer_counts.keys()):
            for domain_type in limiting_domain_types:
                if monomer_type.net_count(domain_type) > 0:
                    yield monomer_type
                    break  # consider next monomer now; do not need to compare other domains on this monomer

    def count(self, monomer: Monomer) -> Union[int, float]:
        return self.__monomer_counts.get(monomer, 0)

    def number_of_monomers(self) -> Union[int, float]:
        return sum(self.__monomer_counts.values())
