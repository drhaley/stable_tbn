import re
from typing import Dict, List
from source.domain import Domain
from source.positive_multiset import PositiveMultiset


class Monomer:
    name_regex = r"[A-Za-z0-9_]+"
    multiple_domain_regex = f"(?:{Domain.regex()}|[1-9]\\d*\\(\\s*{Domain.regex()}\\s*\\))"

    __known_monomers = {}  # used to catalogue monomers; used to makes sure that all monomer names are unique

    def __init__(self, domain_counts: Dict[Domain, int], name: str):
        if not domain_counts:
            raise AssertionError("attempted to create an empty monomer")

        self.__domain_counts = PositiveMultiset(Domain, domain_counts)

        domain_strings_as_list = []
        for domain in sorted(domain_counts.keys()):
            count = domain_counts[domain]
            if count > 1:
                domain_as_string = f"{count}({domain})"
            else:
                domain_as_string = str(domain)

            domain_strings_as_list.append(domain_as_string)
        domains_as_string = " ".join(domain_strings_as_list)

        if name is None:
            final_name = f"[{domains_as_string}]"
        else:
            stripped_name = name.strip()
            if stripped_name == "":
                raise AssertionError("Cannot give whitespace string as a name for a monomer")
            else:
                final_name = stripped_name

        if final_name in self.__known_monomers:
            previously_known_monomer = self.__known_monomers[final_name]
            if domain_counts != previously_known_monomer.__domain_counts:
                raise AssertionError(f"Cannot have two distinct monomers with the same name: {final_name}")

        self.__name = final_name
        self.__known_monomers[final_name] = self

    def __str__(self) -> str:
        return self.__name

    def __eq__(self, other: "Monomer") -> bool:
        return self.__name == other.__name  # relies on __init__() to assert that the names are unique

    def __lt__(self, other: "Monomer") -> bool:
        return str(self) < str(other)

    def __hash__(self) -> int:
        return hash(str(self))

    @classmethod
    def from_string(cls, monomer_as_string: str, name: str = None) -> "Monomer":
        # name extraction first
        domain_list_regex = f"{cls.multiple_domain_regex}(?: {cls.multiple_domain_regex})*"
        
        name_search_pattern = f"^({domain_list_regex})\\s*(|>{cls.name_regex})$"
        name_search_result = re.match(name_search_pattern, monomer_as_string)
        if not name_search_result:
            raise AssertionError(f"could not parse monomer from string '{monomer_as_string}'")

        composition_string, raw_name_string = name_search_result.groups()
        
        if raw_name_string:
            if name:
                raise AssertionError(
                    "received call to Monomer.from_string() specifying a name in the string and in the passed argument"
                )
            else:
                name = raw_name_string[1:]  # remove the '>'

        # now parse the composition of the monomer
        domain_counts = {}
        domain_strings_as_list = composition_string.split()
        for domain_string_with_optional_quantity in domain_strings_as_list:
            quantity_search_result = re.match(
                f"^([1-9]\\d*)\\(\\s*({Domain.regex()})\\s*\\)$",
                domain_string_with_optional_quantity
            )
            if quantity_search_result:
                count = int(quantity_search_result.groups()[0])
                domain_string = quantity_search_result.groups()[1]
            else:
                count = 1
                domain_string = domain_string_with_optional_quantity

            domain_search_result = re.match(f"^({Domain.regex()})$", domain_string)
            if not domain_search_result:
                raise AssertionError(f"Could not parse domain name from {domain_search_result}")

            domain = Domain(domain_search_result.groups()[0])
            domain_counts[domain] = domain_counts.get(domain, 0) + count

        return cls(domain_counts, name)

    def name(self) -> str:
        return self.__name

    def unstarred_domain_types(self) -> List[Domain]:
        unstarred_domain_types = set()
        for domain_type in self.__domain_counts:
            if domain_type.is_starred():
                unstarred_domain_type = domain_type.complement()
            else:
                unstarred_domain_type = domain_type
            unstarred_domain_types.add(unstarred_domain_type)
        return sorted(unstarred_domain_types)

    def net_count(self, domain: Domain) -> int:
        return self.__domain_counts.get(domain, 0) - self.__domain_counts.get(domain.complement(), 0)

    @classmethod
    def regex(cls):
        domain_list_regex = f"{cls.multiple_domain_regex}(?: {cls.multiple_domain_regex})*"
        return f"{domain_list_regex}(?:|\\s*>{cls.name_regex})"
