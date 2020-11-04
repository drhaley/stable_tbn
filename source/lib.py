from typing import Iterator, Optional
from source.solver import Solver, SolverMethod, SolverFormulation
from source.configuration import Configuration
from source.tbn import Tbn
from source.constraints import Constraints


def get_stable_configs(
        tbn_filename: str,
        constraints_filename: Optional[str] = None,
        solver_method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
        formulation: SolverFormulation = SolverFormulation.POLYMER_UNBOUNDED_MATRIX,
        bond_weighting_factor: Optional[float] = None,
        verbose: bool = False,
) -> Iterator[Configuration]:
    tbn = get_tbn_from_filename(tbn_filename)
    user_constraints = get_constraints_from_filename(constraints_filename)
    solver = Solver(method=solver_method)
    stable_configurations = solver.stable_configs(
        tbn,
        user_constraints=user_constraints,
        formulation=formulation,
        bond_weighting_factor=bond_weighting_factor,
        verbose=verbose,
    )

    return stable_configurations


def get_stable_config(
        tbn_filename: str,
        constraints_filename: Optional[str] = None,
        solver_method: SolverMethod = SolverMethod.CONSTRAINT_PROGRAMMING,
        formulation: SolverFormulation = SolverFormulation.POLYMER_UNBOUNDED_MATRIX,
        bond_weighting_factor: Optional[float] = None,
        verbose: bool = False,
) -> Configuration:
    tbn = get_tbn_from_filename(tbn_filename)
    user_constraints = get_constraints_from_filename(constraints_filename)
    solver = Solver(method=solver_method)
    stable_configuration = solver.stable_config(
        tbn,
        user_constraints=user_constraints,
        formulation=formulation,
        bond_weighting_factor=bond_weighting_factor,
        verbose=verbose,
    )

    return stable_configuration


def get_tbn_from_filename(tbn_filename) -> Tbn:
    with open(tbn_filename) as tbnFile:
        tbn_as_string = tbnFile.read()

    return Tbn.from_string(tbn_as_string)


def get_constraints_from_filename(constraints_filename) -> Constraints:
    if constraints_filename:
        with open(constraints_filename) as constraintsFile:
            constraints_as_string = constraintsFile.read()
        return Constraints.from_string(constraints_as_string)
    else:
        return Constraints()
