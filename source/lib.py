import os
from typing import Iterator, Optional
from source.solver import Solver, SolverMethod, SolverFormulation
from source.configuration import Configuration
from source.tbn import Tbn


def get_stable_configs(
        tbn_filename: str,
        constraint_filename: Optional[str] = None,
        solver_method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
        formulation: SolverFormulation = SolverFormulation.BEYOND_MULTISET_FORMULATION,
        bond_weighting_factor: Optional[float] = None,
        verbose: bool = False,
    ) -> Iterator[Configuration]:

    tbn = get_tbn_from_filename(tbn_filename)
    _ = get_constraints_from_filename(constraint_filename)
    solver = Solver(method=solver_method)
    stable_configurations = solver.stable_configs(
        tbn, formulation=formulation, bond_weighting_factor=bond_weighting_factor, verbose=verbose
    )

    return stable_configurations


def get_stable_config(
        tbn_filename: str,
        constraint_filename: Optional[str] = None,
        solver_method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
        formulation: SolverFormulation = SolverFormulation.BEYOND_MULTISET_FORMULATION,
        bond_weighting_factor: Optional[float] = None,
        verbose: bool = False,
    ) -> Configuration:

    tbn = get_tbn_from_filename(tbn_filename)
    _ = get_constraints_from_filename(constraint_filename)
    solver = Solver(method=solver_method)
    stable_configuration = solver.stable_config(
        tbn, formulation=formulation, bond_weighting_factor=bond_weighting_factor, verbose=verbose
    )

    return stable_configuration


def get_tbn_from_filename(tbn_filename) -> Tbn:
    with open(tbn_filename) as tbnFile:
        tbn_as_string = tbnFile.read()

    return Tbn.from_string(tbn_as_string)


def get_constraints_from_filename(constraint_filename) -> str:
    if constraint_filename:
        raise NotImplementedError("constraint file is not yet implemented")
        # with open(constraint_filename) as constraintFile:
        #     constraints_as_string = constraintFile.read()
    return ""
