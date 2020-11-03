from abc import ABC, abstractmethod
from typing import Any, Iterator, Dict, List


class Model(ABC):
    def __init__(self):
        ABC.__init__(self)
        self._big_M = None
        self.OPTIMAL = "ABSTRACT CLASS OPTIMAL"
        self.INFEASIBLE = "ABSTRACT CLASS INFEASIBLE"

    @abstractmethod
    def int_var(self, *args, **kargs) -> Any:
        pass

    @abstractmethod
    def bool_var(self, *args, **kargs) -> Any:
        pass

    @abstractmethod
    def complement_var(self, var: Any) -> Any:
        pass

    @abstractmethod
    def add_constraint(self, *args) -> Any:
        pass

    @abstractmethod
    def add_implication(self, *args) -> Any:
        # intended call: .add_implication(antecedent1, antecedent2, ..., consequent)
        pass

    @abstractmethod
    def add_equal_to_zero_implication(self, *args) -> Any:
        # same as add_implication, except last expression will be enforced to be equal to zero
        pass

    @abstractmethod
    def add_greater_than_zero_implication(self, *args) -> Any:
        # same as add_implication, except last expression will be enforced to be equal to zero
        pass

    @abstractmethod
    def minimize(self, *args, **kargs) -> None:
        pass

    @abstractmethod
    def maximize(self, *args, **kargs) -> None:
        pass

    def set_big_m(self, big_M: int) -> None:
        # not used by all solvers.  this should be a large value (i.e. for big M formulations for integer programming)
        self._big_M = big_M


class SolverAdapter(ABC):
    def __init__(self):
        super().__init__()

    @staticmethod
    @abstractmethod
    def model() -> Model:
        pass

    @abstractmethod
    def solve(self, model: Model, variables_with_values_to_keep: List[Any], verbose: bool = False) -> Any:
        pass

    @abstractmethod
    def value(self, var: Any) -> int:
        pass

    @abstractmethod
    def solve_all(self, model: Model, variables_with_values_to_keep: List[Any], verbose: bool = False)\
            -> Iterator[Dict[Any, int]]:
        pass
