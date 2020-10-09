import re
from copy import deepcopy


class Domain:
    name_regex = r"[A-Za-z0-9_]+"
    optional_star_regex = r"(?:\*|)"
    optional_assigned_name_regex = r"(?:\:[A-Za-z0-9_]+|)"

    def __init__(self, domain_as_string: str):
        # alphanumerics with underscores and an optional star.
        # for backwards compatibility with StableGen, also allows a colon and an assigned name (which are ignored)

        search_pattern = f"^({self.name_regex})({self.optional_star_regex}){self.optional_assigned_name_regex}$"
        name_search_result = re.match(search_pattern, domain_as_string)
        if not name_search_result:
            parsing_error_message = f"could not parse domain: '{domain_as_string}', format must be '{self.regex()}'"
            raise AssertionError(parsing_error_message)

        self.__name, optional_star = name_search_result.groups()
        self.__starred = True if optional_star == "*" else False

    def __str__(self) -> str:
        if self.__starred:
            return self.__name + "*"
        else:
            return self.__name

    def __eq__(self, other: "Domain") -> bool:
        return (self.__starred == other.__starred) and (self.__name == other.__name)

    def __lt__(self, other: "Domain") -> bool:
        if self.__name < other.__name:
            return True
        elif (self.__name == other.__name) and not self.is_starred() and other.is_starred():
            return True
        else:
            return False

    def __hash__(self) -> int:
        return hash(str(self))

    def is_starred(self) -> bool:
        return self.__starred

    def complement(self) -> "Domain":
        new_domain = deepcopy(self)
        new_domain.__starred = not self.__starred
        return new_domain

    @classmethod
    def regex(cls):
        return f"{cls.name_regex}{cls.optional_star_regex}{cls.optional_assigned_name_regex}"
