from typing import List, Tuple

from source.monomer import Monomer
from source.formulations.polymer_integer_matrix import Formulation as IntegerFormulation


class Formulation(IntegerFormulation):
    def _get_monomer_types_and_counts(self) -> Tuple[List[Monomer], List[int]]:
        # all monomers are labelled in this formulation (no multisets)
        ordered_monomer_types = list(self.tbn.monomer_types(flatten=True))
        monomer_counts = [1 for _ in ordered_monomer_types]
        return ordered_monomer_types, monomer_counts

    def _run_asserts(self) -> None:
        super()._run_asserts()
