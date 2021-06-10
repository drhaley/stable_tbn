import argparse
import re
from typing import Dict, List


def main() -> None:
    args = get_command_line_arguments()

    def print_verbose(output: str) -> None:
        if args.verbose:
            print(output)

    if not args.filename:
        print("Syntax: python3 crn_conversion.py <crn_text_filename>")
    else:
        crn = {}
        w = None
        with open(args.filename) as inFile:
            for i, line in enumerate(inFile):
                line = line.strip()
                find_w = re.findall(r"^w\s*=\s*(\d*)$", line)  # why do re.match and re.search return the whole string?
                if not line:
                    pass
                elif find_w:
                    # w = find_w.group()
                    w = int(find_w[0])
                    print_verbose(f"w = {w}")
                elif not re.match(r"^(\w+)\s*:(\s*(\w+\^[T|F|t|f]))*$", line):
                    raise AssertionError(f"Syntax error in line {i}:\n{line}")
                else:
                    print_verbose(f"\nLINE {i}")
                    preline, postline = line.split(':')
                    reactant = re.search(r"^(\w+)\s*", preline).group()
                    catalysts = re.findall(r"(\w+\^[T|F|t|f])", postline)
                    print_verbose(f"reactant: {reactant}")
                    print_verbose(f"catalysts: {catalysts}")
                    if reactant in crn:
                        raise AssertionError(f"duplicate declaration of reactant \"{reactant}\"")
                    else:
                        crn[reactant] = catalysts
        tbn = process_crn(crn, w, verbose=args.verbose)
        with open("processed_tbn.txt", "w") as outFile:
            outFile.write(tbn)


def get_command_line_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filename",
        metavar="filename",
        type=str,
        help="filename for crn text file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="verbose output",
    )
    return parser.parse_args()


def process_crn(crn: Dict[str, List[str]], w: int, verbose:bool=False) -> str:
    def print_verbose(output: str) -> None:
        if verbose:
            print(output)

    tbn = {}
    k = max(len(catalysts) for catalysts in crn.values())
    k = max(k, 2)
    n= k*w

    for reactant in crn:
        tbn[f"{reactant}_G"] = [f"{reactant}_domain_{i}{j}*" for i in range(n) for j in range(n)]
        print_verbose(f"{reactant}_G: " + str(tbn[f"{reactant}_G"]))
        for i in range(n):
            tbn[f"{reactant}_H{i}"] = [f"{reactant}_domain_{i}{j}" for j in range(n)]
            print_verbose(f"{reactant}_H{i}: " + str(tbn[f"{reactant}_H{i}"]))
        for j in range(n):
            tbn[f"{reactant}_V{j}"] = [f"{reactant}_domain_{i}{j}" for i in range(n)]
            print_verbose(f"{reactant}_V{j}: " + str(tbn[f"{reactant}_V{j}"]))

    for reactant, catalysts in crn.items():
        diag = 0
        for catalyst in catalysts:
            name, raw_version = catalyst.split("^")
            version = "H" if raw_version.upper() == "F" else "V"
            for i in range(w):
                new_domains = [f"{reactant}_domain_{j}{diag-j}" for j in range(diag+1)]
                tbn[f"{name}_{version}{i}"].extend(new_domains)
                diag += 1

        for missing_cat_number in range(k - len(catalysts)):
            monomer_name = f"{reactant}_C{len(catalysts) + missing_cat_number}"
            for i in range(w):
                new_domains = [f"{reactant}_domain_{j}{diag-j}" for j in range(diag+1)]
                tbn.setdefault(monomer_name, []).extend(new_domains)
                diag += 1

    tbn_as_string = '\n'.join([f"{' '.join(domains)} >{monomer}" for monomer, domains in tbn.items()])
    return tbn_as_string


if __name__ == "__main__":
    main()
