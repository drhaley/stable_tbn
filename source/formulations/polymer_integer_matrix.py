from typing import List, Tuple
from math import inf as infinity

from source.monomer import Monomer
from source.formulations.polymer_unbounded_matrix import Formulation as UnboundedFormulation


class Formulation(UnboundedFormulation):
    def _get_limiting_monomer_types(self) -> List[Monomer]:
        return list(self.tbn.limiting_monomer_types())

    def _run_asserts(self) -> None:
        super()._run_asserts()
        if self.monomer_counts == infinity:
            raise AssertionError("Cannot run integer matrix formulation on a tbn with infinite monomer counts")
