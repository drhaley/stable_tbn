import re
from math import inf as infinity
from copy import copy


class Constraints:
    """
    User-defined constraints are placed into this container class
    """
    def __init__(self):
        self._max_polymers = infinity
        self._min_polymers = 0
        self._max_merges = infinity
        self._min_merges = 0
        self._max_energy = infinity
        self._min_energy = -infinity
        self._sort = True
        self._optimize = True
        self._bond_weight = 2.0

    @classmethod
    def from_string(cls, text: str) -> "Constraints":
        this_constraints = Constraints()
        for raw_line in text.split('\n'):
            line = raw_line.strip()
            if line:  # not just whitespace
                this_constraints.add_new_constraint_from_string(line)
        return this_constraints

    def add_new_constraint_from_string(self, line: str) -> None:
        # TODO: refactor the boiler plate switch structure
        line = line.upper()

        unsigned_floating_point_regex = r"(?:[0-9]*[.])?[0-9]+"
        floating_point_regex = r"[+-]?(?:[0-9]*[.])?[0-9]+"
        non_negative_integer_regex = r"[0-9]+"

        if re.match("OPTIMIZE", line):
            self._optimize = True
        elif re.match("NO\\s+OPTIMIZE", line):
            self._optimize = False
        elif re.match("SORT", line):
            self._sort = True
        elif re.match("NO\\s+SORT", line):
            self._sort = False
        else:
            still_searching = True
            search_results = re.match(f"MAX ENERGY ({floating_point_regex})", line)
            if still_searching and search_results:
                self._max_energy = float(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"MIN ENERGY ({floating_point_regex})", line)
            if still_searching and search_results:
                self._min_energy = float(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"MAX MERGES ({non_negative_integer_regex})", line)
            if still_searching and search_results:
                self._max_merges = int(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"MIN MERGES ({non_negative_integer_regex})", line)
            if still_searching and search_results:
                self._min_merges = int(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"MAX POLYMERS ({non_negative_integer_regex})", line)
            if still_searching and search_results:
                self._max_polymers = int(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"MIN POLYMERS ({non_negative_integer_regex})", line)
            if still_searching and search_results:
                self._min_polymers = int(search_results.groups()[0])
                still_searching = False
            search_results = re.match(f"BOND WEIGHT ({unsigned_floating_point_regex})", line)
            if still_searching and search_results:
                self._bond_weight = float(search_results.groups()[0])
                still_searching = False
            if still_searching:
                raise AssertionError(f"Cannot parse line '{line}' in constraints file")

    def with_fixed_polymers(self, number_of_polymers: int) -> "Constraints":
        this = copy(self)
        this._max_polymers = number_of_polymers
        this._min_polymers = number_of_polymers
        return this

    def with_fixed_merges(self, number_of_merges: int) -> "Constraints":
        this = copy(self)
        this._max_merges = number_of_merges
        this._min_merges = number_of_merges
        return this

    def with_fixed_energy(self, amount_of_energy: float) -> "Constraints":
        this = copy(self)
        this._max_energy = amount_of_energy
        this._min_energy = amount_of_energy
        return this

    def with_bond_weight(self, bond_weight: float) -> "Constraints":
        this = copy(self)
        this._bond_weight = bond_weight
        return this

    def with_unset_optimization_flag(self) -> "Constraints":
        this = copy(self)
        this._optimize = False
        return this

    def max_polymers(self):
        return self._max_polymers

    def min_polymers(self):
        return self._min_polymers

    def max_merges(self):
        return self._max_merges

    def min_merges(self):
        return self._min_merges

    def max_energy(self):
        return self._max_energy

    def min_energy(self):
        return self._min_energy

    def sort(self):
        return self._sort

    def optimize(self):
        return self._optimize

    def bond_weight(self):
        return self._bond_weight
