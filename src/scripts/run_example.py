#!/usr/bin/env python
import argparse

from aspect_py_project_template.example import get_massnumber


def main():
    parser = argparse.ArgumentParser(
        description="Retrun the mass number of the given molecule."
    )
    parser.add_argument(
        "--molecule", type=str, required=True, help="Valid molecule formula (E.g. H2O)."
    )
    args = parser.parse_args()
    mass_number = get_massnumber(args.molecule)
    print(f"The mass number of {args.molecule} is {mass_number}")


if __name__ == "__main__":
    main()
