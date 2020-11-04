from math import inf as infinity

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

    @classmethod
    def from_string(cls, text: str) -> "Constraints":
        this_constraints = Constraints()
        for raw_line in text.split('\n'):
            line = raw_line.strip()
            if line:  # not just whitespace
                this_constraints.add_new_constraint_from_string(line)
        return this_constraints

    def add_new_constraint_from_string(self, line: str) -> None:
        pass  # TODO: make a syntax

    def set_fixed_polymers(self, number_of_polymers: int):
        self._max_polymers = number_of_polymers
        self._min_polymers = number_of_polymers

    def set_fixed_merges(self, number_of_merges: int):
        self._max_merges = number_of_merges
        self._min_merges = number_of_merges

    def set_fixed_energy(self, amount_of_energy: float):
        self._max_energy = amount_of_energy
        self._min_energy = amount_of_energy

    def unset_optimization_flag(self):
        self._optimize = False

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
