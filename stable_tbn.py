import argparse
import timeit
from source.solver import SolverMethod, SolverFormulation
from source import lib


def main() -> None:
    args = get_command_line_arguments()

    tic = timeit.default_timer()

    formulation = SolverFormulation.BEYOND_MULTISET_FORMULATION  # TODO: create command-line arguments to choose

    if not args.single:
        stable_configurations = lib.get_stable_configs(
            tbn_filename=args.tbn_filename,
            constraint_filename=args.constraint_filename,
            solver_method=SolverMethod.INTEGER_PROGRAMMING if args.ip else SolverMethod.CONSTRAINT_PROGRAMMING,
            formulation=formulation,
        )

        toc = timeit.default_timer()
        for i, configuration in enumerate(stable_configurations):
            print(f"Configuration {i+1}:\n{configuration}")
    else:
        stable_configuration = lib.get_stable_config(
            tbn_filename=args.tbn_filename,
            constraint_filename=args.constraint_filename,
            solver_method=SolverMethod.INTEGER_PROGRAMMING if args.ip else SolverMethod.CONSTRAINT_PROGRAMMING,
            formulation=formulation,
        )

        toc = timeit.default_timer()
        print(stable_configuration)

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
            "-c",
            dest="constraint_filename",
            metavar="constraint_filename",
            type=str,
            help="filename for constraint text file",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
