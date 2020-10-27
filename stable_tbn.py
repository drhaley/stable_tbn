import argparse
import timeit
from source.solver import SolverMethod, SolverFormulation
from source import lib

alternate_formulations = [
    "STABLEGEN_FORMULATION",
    "BOND_OBLIVIOUS_FORMULATION",
    "SET_FORMULATION",
    "MULTISET_FORMULATION",
]


def main() -> None:
    args = get_command_line_arguments()

    tic = timeit.default_timer()

    if args.weight is None:
        if args.formulation is not None:
            for formulation_as_string, formulation_as_enum in SolverFormulation.__members__.items():
                if formulation_as_string == args.formulation:
                    formulation = formulation_as_enum
                    break
            else:
                raise AssertionError(f"Did not recognize alternate formulation method requested '{args.formulation}'")
        else:
            formulation = SolverFormulation.BEYOND_MULTISET_FORMULATION
        bond_weighting_factor = None
    else:
        if args.method is not None:
            print("alternate method requested but bond weight was specified, falling back to LOW_W_FORMULATION")
        formulation = SolverFormulation.LOW_W_FORMULATION
        bond_weighting_factor = float(args.weight)

    if not args.single:
        stable_configurations = lib.get_stable_configs(
            tbn_filename=args.tbn_filename,
            constraint_filename=args.constraint_filename,
            solver_method=SolverMethod.INTEGER_PROGRAMMING if args.ip else SolverMethod.CONSTRAINT_PROGRAMMING,
            formulation=formulation,
            bond_weighting_factor=bond_weighting_factor,
        )

        toc = timeit.default_timer()
        for i, configuration in enumerate(stable_configurations):
            configuration_string = configuration.full_str() if args.full else str(configuration)
            print(f"Configuration {i+1}:\n{configuration_string}")
    else:
        stable_configuration = lib.get_stable_config(
            tbn_filename=args.tbn_filename,
            constraint_filename=args.constraint_filename,
            solver_method=SolverMethod.INTEGER_PROGRAMMING if args.ip else SolverMethod.CONSTRAINT_PROGRAMMING,
            formulation=formulation,
            bond_weighting_factor=bond_weighting_factor,
        )

        toc = timeit.default_timer()
        configuration_string = stable_configuration.full_str() if args.full else str(stable_configuration)
        print(f"Configuration: {configuration_string}")

    if args.timed:
        print(f"seconds elapsed: {toc-tic}")


def get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "tbn_filename",
        metavar="tbn_filename",
        type=str,
        help="filename for tbn text file",
    )
    parser.add_argument(
        '-i',
        "--ip",
        action="store_true",
        help="use the integer programming formulation (instead of constraint programming)",
    )
    parser.add_argument(
        '-1',
        dest="single",
        action="store_true",
        help="only report one stable configuration",
    )
    parser.add_argument(
        "-t",
        "--timed",
        action="store_true",
        help="print elapsed time",
    )
    parser.add_argument(
        "-f",
        "--full",
        action="store_true",
        help="print full configuration, including singletons",
    )
    parser.add_argument(
        "-w",
        "--weight",
        metavar="bond_weight",
        type=str,
        help="energy weight/worth of bonds vs polymers formed, e.g. 0.5",
    )
    parser.add_argument(
        "-c",
        dest="constraint_filename",
        metavar="constraint_filename",
        type=str,
        help="filename for constraint text file",
    )
    parser.add_argument(
        "--formulation",
        type=str,
        help=f"specify an alternate solution formulation, one of:\n{', '.join(alternate_formulations)}",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
